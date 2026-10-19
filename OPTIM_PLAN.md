# fastKarOptim: optimization plan

Branch `optim` of fastKar, checked out as a worktree at `~/git/fastKar_optim`.
The package is renamed `fastKarOptim` (DESCRIPTION, NAMESPACE `useDynLib`,
Rcpp exports, `fastKar::` data references, `tests/testthat.R`) so it loads in the
same R session as the reference `fastKar`. Neither package exports anything
(no `export()` in NAMESPACE), so cross-package code uses `fastKar:::f` and
`fastKarOptim:::f`.

## 1. Baseline (profiled 2026-10-18, one walk set, `mc.cores = 1`)

| Event (CCLE) | Footprint | Walks / walk nodes | pix | Tiles | `forward_simulate` |
|---|---|---|---|---|---|
| pyrgo NCI-H1048 | 0.65 Mb | 52 / 100 | 10 kb | 69 | 0.46 s (0.2 s is `prep_for_sim`) |
| dm NCI-H69 | 2.9 Mb | 380 / 4623 | 100 kb / 10 kb | 70 / 334 | 4.6 s / 18 s |
| chromothripsis NCI-H2066 | 41 Mb | 122 / 1767 | 100 kb / 10 kb | 439 / 4131 | 0.6 s / 38 s |

`compmaps` on the 4131-tile map takes 3.0 s and `make_noisydat` takes 2.0 s.

Where the time goes (Rprof, 10 kb runs):
- **`sum_matrices` + `calculate_interchrom`: ~50%.** Each walk returns a *full copy* of
  the reference template (P = n(n+1)/2 rows; 8.5 M for 4131 tiles). These are
  `rbind`-ed (W x P rows), grouped by `id` and `unique`d. This dominates the
  time, and memory grows as W x P.
- **`simulate_map`: 15-25%.** It builds all local pairs as a data.table, then
  evaluates the spline in R. Pairs whose tiles fall outside the target
  (`orig.id = NA`) are computed and thrown away. That is most pairs in
  full-graph mode.
- **Sparse lift `t(T) %*% L %*% T`, `symmetrize`, `TsparseMatrix` conversion: ~10%.**
- **`data.table::copy` of the template per walk: ~12%.**
- **Fixed per-call overhead ~0.3 s**: `prep_for_sim` (GRanges merges, 0.2 s) and
  `gM()` construction (0.1 s). `bestfit_search` pays this for every candidate,
  e.g. 1000 candidates x 0.3 s = 5 min of pure overhead per event.

Model facts that make an exact C++ port easy:
- `distance_spline` is an R `splinefun` with method fmm and 24,862 knots on
  log(d) in [9.21, 19.33]. Evaluating it is a binary search plus a cubic.
- `diag_spline` is linear: `-9.76352525944288 + 0.611138644669642 * log(widthprod)`.
- `interchrom_density` is a scalar (2.8e-13).

## 2. Ground rules

1. **Bit-faithful first.** Each optimization must reproduce the current
   `fastKar` outputs. Integer and structural outputs (`i`, `j`, `id`, tile tables,
   hashes, walks) must be `identical()`. Floating-point values must agree to a
   relative 1e-12, with absolute slack for values below 1e-300. Only summation
   order may change.
2. **Keep current behaviour, quirks included.** Bugs found along the way
   (list in section 7) are fixed in separate, labelled commits, each with its own
   test, and only after the equivalence suite passes on the faithful version. That
   way "matches the reference" always has a clear meaning.
3. **RNG streams are preserved by default.** With the same `set.seed`, a function
   must make the same draws in the same order as the reference
   (`sample.gwalks`, `local.sampling`, `make_noisydat`, `paste_loose_ends`).
   Where breaking the stream buys a large speedup, the stream-preserving version
   stays the default and the fast one goes behind an option. This needs your
   sign-off first.
4. **Public signatures and return shapes do not change.** New fast paths are
   either internal or new optional arguments with backwards-compatible
   defaults.
5. **One change per commit, with timings before and after**, recorded by the
   benchmark harness (section 6).

## 3. forward.R (main target)

### F1. Accumulate into one reference vector (pure R, exact)
- In `simulate_walks`, each walk's contribution goes into a single numeric
  vector `acc` of length P, indexed by pixel `id`. The walk's local sparse result
  is added to `acc` via `id = (i-1)*n - (i-1)*(i-2)/2 + (j-i+1)`. The template
  is no longer copied per walk.
- `sum_matrices` becomes "add vectors". `calculate_interchrom` becomes
  `acc + widthprod * rho * cn[i] * cn[j]` on the template, with no merge and
  no `allow.cartesian`.
- `if.sum = FALSE` still returns one data.table per walk. It is built only at
  the end, from per-walk sparse vectors.
- Expected effect: removes ~60% of the time and the W x P memory term.

### F2. C++ walk kernel (exact)
`sim_walk_ref_cpp(mid, width, orig_id, circular, total_len, n_ref, acc, model, cn)`
does everything `simulate_map` + lift + symmetrize + upper-triangle extraction
does now, in one pass:
- It loops only over local pairs where both tiles are in the target
  (`orig.id` not NA). Distances still come from the full walk's `cumsum`.
- Lift semantics are copied exactly from `t(T) L T` on the symmetrized local matrix.
  For local `a < b` mapped to reference `p, q`:
  - if `p != q`, add `L_ab` to `(min, max)`;
  - if `p == q` (a revisit), add `2 L_ab`;
  - for `a == a`, add `L_aa`.
- It evaluates the spline in C++: an fmm cubic with the same interval search and
  extrapolation as R's `spline_eval`; the below-10 kb power law with the knot
  slope; the 1/d tail above 1e8; the linear `diag_spline`. Long-read model:
  `exp(-d/model)`. Spline coefficients are read once at load from `splineobj`
  (so retrained splines are picked up automatically), never hard-coded.
- `make_template_dat_cpp` moves from a load-time `cppFunction()` into
  `src/`. Today it is compiled at load time, which slows every
  `load_all` and breaks on machines without a compiler.
- Expected effect: the remaining `simulate_map` time drops by 10-50x, and
  full-graph mode gains the most.

### F3. Pay per-event costs once
- `forward_simulate(..., prepped = NULL)`: accepts a precomputed
  `prep_for_sim()` result. `bestfit_search`, `hic_kl`/`get_dists`,
  `group_klcomp` and `neighbor_separable` prepare the tiling once per (graph,
  target, pix) and reuse it. The tiling does not depend on the walks.
- `gm.out = FALSE` paths in internal callers skip `gM()` construction.
- `get_dists`/`get_all_pair_distances`/`group_klcomp`: simulate each karyotype
  once, then compute all pairwise KLs. Today `hic_kl` re-simulates both maps
  for every pair, so n karyotypes cost n(n-1) sims; this cuts it to n.
- `prep_for_sim` itself: replace the GRanges merges with an interval join on
  integers where it is safe. This is low priority once F3 has removed the
  repeated calls.

### F4. Candidate sets that share walks (optional, exact)
Candidates from one graph often share most walks. A per-call cache, keyed by
canonical walk hash and walk CN, stores each walk's sparse reference
contribution, so `bestfit_search` and `identifiability_hic_nlls` simulate each
distinct walk once. The cache is opt-in and has a memory cap.

### F5. Parallelism
After F1-F2, `mclapply` over walks inside `simulate_walks` costs more in fork
overhead than it saves. Plan:
- Make the single-core path fast.
- Parallelize over candidates (callers already do this).
- Optionally add OpenMP inside the kernel over walks, with per-thread
  accumulators and a deterministic reduction order, so results do not depend on
  the thread count.

## 4. comparisons.R

- **C1 `compdats`.** Add a fast path for the common case where both inputs are
  on the same template. In that case it works on aligned vectors directly; today
  it does a `setkey` + `all.equal` + merge.
  - Keep `round(value.x)`, the `value.y == 0 -> 0` rule, the `widthprod`/`area0`
    weighting and `na.rm`.
  - Keep the `setkey` side effect on the caller's tables, because callers may
    rely on the sorted order (`identifiability_hic_nlls` passes `copy(obs_dat)`
    for this reason). It is removed only as a flagged fix.
- **C2 `make_noisydat`.** Keep the single vectorized `rnbinom`/`rpois` call, so
  the RNG stream is unchanged. Build the `nsamp` tables by slicing one value
  vector instead of `rbind`-ing `nsamp` copies and calling `split`.
- **C3 `compmaps`/`rebin_matrix`.** Skip `make_template_dat` when the input is
  already on the target bins. Cache templates by bin count.
- **C4 `liktest_separable(_lr)`, `longread_probdist`.** Precompute the reference
  side once. Replace string-named vector indexing (`setNames(...)[all_words]`)
  with `match()`.
- `estimate.depthratio` (I/O bound) and `callloops` (rarely used) are left alone.

## 5. graphstats.R

- **G1 `gg.to.wiring`.** Remove the string pasting, `grepl` and `substr` on node
  ids and use integer side codes. Callers can pass a precomputed `wiring`
  (`sample.gwalks`, `local.sampling`, `sample.neighborhood`, `markov.gwalk`).
- **G2 `sample.gwalks`.**
  - Precompute the shufflable groups once (non-frozen, size > 1), so the
    per-sample R loop only visits them. Draws stay `sample(out[g])` in the same
    group order, so the RNG stream is unchanged.
  - Replace `digest()` of each permutation with a C++ hash of the integer vector
    (it is only used for dedup).
  - `gW()` construction per walk is the remaining cost. Measure it before
    deciding anything.
- **G3 `local.sampling` and relatives.** Replace the per-step
  `data.table::copy` and `edges[n == pivot]` with integer vectors and
  precomputed group indices. Keep the exact `sample.int`/`sample` draws.
  `traverse_graph_cpp` is ~2.2 s per call on high-copy graphs: profile the C++
  separately (incremental re-traversal after a 2-swap is the candidate idea).
- **G4 `reads_fromwalk`.** Replace `foverlaps` per walk and `paste` of words
  with a C++ sweep that emits the same words and multiplicities. Words must be
  `identical()`. Row order is canonicalized before comparison, since callers
  key by word.
- **G5 `edit_dist_cpp`/`split_walk_fragments`.** `split_walk_fragments` grows
  `out = c(out, pieces)` inside a loop. Preallocate it, or move it to C++ next to
  `compute_cost_matrix_cpp`. `solve_LSAP` is left as is.

## 6. Tests

All tests live in `tests/testthat/` and run with
`setDTthreads(1); devtools::test('~/git/fastKar_optim')`.

### T0. Pinned reference oracle
- Install the reference fastKar **at commit `0a5497e`** into a private library
  (`~/git/fastKar_optim/tests/oracle_lib/`, git-ignored), so the oracle does not
  move when fastKar master changes.
- Install it with the untracked `R/identifiability.R` included, and record the
  commit and file hashes in `tests/oracle_lib/ORACLE_VERSION`.
- Helper `skip_if_no_oracle()` skips the test if the library is missing.

### T1. Golden fixtures (fast, no oracle needed)
`tests/testthat/fixtures/*.rds`: small inputs plus the reference outputs,
generated once from the oracle by `tests/make_golden.R`.

Fixture graphs:
- toy double translocation (`makewalk`);
- `make_threeway_repdup`;
- `makerepdup`;
- `gimme_bfb` (high CN, foldbacks);
- a circular-walk graph;
- a graph with 2-3 bp nodes revisited back to back (the spline-blowup case);
- local subgraphs (`lk_local_graph`-style, cut to the footprint, small enough to
  check in) of the three profiled CCLE events plus one tyfonas.

### T2. Differential tests against the oracle (`test-equiv-forward.R`, etc.)
For each fixture graph, a seeded set of karyotypes is run through both packages
over a grid of arguments:
- `pix.size` in {1e4, 1e5, 1e6};
- `if.sum` in {T, F} x `gm.out` in {T, F};
- `model` in {0, 1e5};
- `if.interchr` in {T, F};
- `purity` in {1, 0.6, 0};
- target equal to the footprint, a strict sub-window, a window extending past
  the graph (`node.id = NA` tiles), and several disjoint ranges on different
  chromosomes;
- walks with CN in `walks$dt$cn`, walks on the minus strand, and walks that
  miss the target;
- a walk whose local map is empty.

Comparison: `gr`/tile tables and `i`/`j`/`id` must be `identical()`; values
must agree within the tolerance in section 2.

Covered the same way:
- `prep_for_sim`, `simulate_diploid_map`, `make_template_dat`;
- `compdats`/`compmaps` (NB and Poisson, KL mode, `if.diag`, masks);
- `make_noisydat`, with the same seed giving identical draws;
- `kl_nb`, `rebin_matrix` (both paths), `hic_cov`;
- `sample.gwalks`, with the same seed giving identical walk lists and hashes,
  for frozen nodes, `remove.dups`, `onlyhash`, `keep.circular` and `chunksize`
  boundaries;
- `local.sampling`, `sample.neighborhood`, `gg.to.wiring`;
- `reads_fromwalk` (linear and circular; `minsize`; `use.nodes`; read length
  longer than a walk);
- `longread_probdist`, `edit_dist_cpp` (all `edit_mode`s);
- `paste_loose_ends`.

### T3. Kernel unit tests
- C++ spline vs `splinefun` on a dense log grid. The grid covers 1 bp-1e9 bp,
  every knot, and points just either side of each knot, the 10 kb join and the
  1e8 tail. Required agreement: 1e-13 relative.
- Lift semantics on hand-built walks with exact expected matrices: a revisit
  mapping to the same tile (the `2 L_ab` case), reverse-strand tile order, a
  circular wrap distance, and a mix of in-target and out-of-target tiles.
- The pixel `id` formula vs `make_template_dat_cpp` for n = 1..200.

### T4. Property tests (independent of the oracle; catch errors both share)
- Output maps are symmetric-consistent, non-negative and finite.
- Linear in `depth` and in walk CN. Summing the `if.sum = FALSE` maps plus the
  interchromosomal term gives the `if.sum = TRUE` map.
- A walk and its reverse complement `-rev(w)` give the same map. A circular walk
  gives the same map under rotation.
- `purity` mixing is linear: the `p` map equals `p` x the tumor map plus
  (1-`p`) x the diploid map.
- Splitting the target into disjoint pieces leaves each pixel's value unchanged
  wherever both maps define it.
- Results do not depend on `mc.cores` (1 vs 4) or on `setDTthreads`.

### T5. Downstream end-to-end regression (slower, `skip_on_ci()`-style flag)
On two CCLE fixture events, with fixed seeds, the following must be identical
to the oracle in candidate order, with NLLs to 1e-10:
- `bestfit_search(topk = -1)`;
- `identifiability_hic_nlls`;
- `local_hic_nll_matrix` / `local_longread_nll_matrix`;
- `run_graph_local_identifiability` (picks table identical).

### T6. Benchmarks (not testthat)
`bench/bench_forward.R` times every workhorse on the fixture events (the three
profiled events above plus one tyfonas) at 1e5 and 1e4. It records median of 5
runs and peak RSS per (function, event, commit), appended to
`bench/history_forward.tsv`, as the existing `bench/` does for traversal. Each
optimization commit adds a row, and a slowdown above 10% on any case blocks it.

## 7. Bugs found while reading (keep for now; fix later as flagged commits)
- `forward_simulate` purity branch: `combined.dat[is.na(value.normal), normal.value := 0]`
  writes to a new column instead of `value.normal`. It is harmless today, because
  both sides cover every pixel.
- `compdats` calls `setkey()` on its inputs, which reorders the caller's tables
  in place.
- `f1score_comparemaps` mask branch: `gM(..., dat = dat)` refers to an undefined
  `dat`.
- `hic_kl`/`longread_kl` call `error()`, which does not exist (it should be `stop()`).
- `neighbor_separable`: `frozen.nodes` is undefined when `use.nodes` is `NULL`.
- `sample.neighborhood`/`markov.gwalk`: `sample(edges[cn>1]$n, 1)` has the
  single-value `sample(x)` trap that `local.sampling` already fixed.
- `prep_for_sim` docs say `pix.size = 0` gives node-half tiles. Unverified;
  test it under T2.

## 8. Order of work
1. T0 + T1 + T3 scaffolding and the benchmark harness, with no optimization yet. The
   suite passes trivially on the unchanged code, which validates the harness.
2. F1 → F2 → F3 (forward.R). Run T2/T4/T5 and the benchmark after each step.
3. C1-C3, then G1-G3, then G4/G5 and F4/F5.
4. Bug fixes from section 7, each as its own commit.
5. Merge back into fastKar. Either the `optim` branch is merged into `master`
   (rename reverted), or the package name stays separate for a while. That is
   your call.
