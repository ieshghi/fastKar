```
  __          _   _   __
 / _|        | | | | / /
| |_ __ _ ___| |_| |/ /  __ _ _ __
|  _/ _` / __| __|    \ / _` | '__|
| || (_| \__ \ |_| |\  \ (_| | |
|_| \__,_|___/\__\_| \_/\__,_|_|
```

Simulate Hi-C data for any rearranged karyotype — fast.

`fastKar` is an R package for sampling and analysing **karyotype walks** through
structurally rearranged cancer genomes, plus forward-modeling Hi-C and long-read
data from those karyotypes. It builds on the [gGnome](https://github.com/mskilab-org/gGnome)
framework.

The 3D model is parameterised from the high-depth Hi-C library published by
[Rao et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/25497547/).

---

## Tutorial — the 5-step arc

The full walkthrough lives in [`fastKar.ipynb`](fastKar.ipynb). What follows is
the headline version, anchored on a synthetic **three-way reciprocal
translocation** between chromosomes 1, 2, and 3.

### 1 · Build a graph

```r
library(fastKar)
gw <- make_threeway_repdup(regsize = 1e6)
gg <- gw$graph
gg
#> gGraph with 9 nodes, 6 loose ends, and 12 edges (6 REF and 6 ALT)
```

Each of chromosomes 1, 2, 3 is tiled into segments `{A,B,C}`, `{D,E,F}`,
`{G,H,I}`. Three reciprocal translocations form a directed cycle
`B → E → H → B`, producing derivative chromosomes alongside the three
normal ones.

![graph](man/figures/fig_graph.png)

Middle segments (B, E, H) sit at CN = 3; outer segments at CN = 2; six loose
ends mark the chromosomal termini.

---

### 2 · Sample karyotypes

```r
walks_gw <- sample.gwalks(gg, N = 500, return.gw = TRUE,
                          remove.dups = TRUE, verbose = FALSE)
length(walks_gw)
#> [1] 8
```

The full karyotype space is enumerable here: 13,824 raw wirings collapse to
8 unique karyotypes after canonical-hashing under cyclic, RC, and multiset
symmetries. Each karyotype is a set of 6 walks covering all 21 node-copies.

![karyotypes](man/figures/fig_karyotypes.png)

Most karyotypes produce six linear chromosomes (three normal + three
derivatives); one of the eight closes the `B → E → H → B` cycle into a small
circular walk.

---

### 3 · Forward-simulate Hi-C

```r
sim <- forward_simulate(walks_gw[[1]], pix.size = 1e5, depth = 1, model = 0)
```

`forward_simulate` predicts Hi-C contacts under the assumption that contact
frequency scales with **linear distance along the karyotype walk**, not along
the reference. Junctions that bring distant reference positions into adjacency
on the walk produce strong off-diagonal Hi-C signal.

![hic](man/figures/fig_hic_single.png)

Different karyotypes look superficially similar — same CN, same intra-chromosome
diagonals — because the karyotype-discriminating signal lives in
**interchromosomal contacts**:

![hic-diff](man/figures/fig_hic_diff.png)

Per-pixel difference maps for `karyo 1 − 2` and `karyo 1 − 3` show that
intra-chromosomal blocks carry zero `Δ`; the differences are entirely between
chromosomes, exactly where the translocation junctions route contacts.

---

### 4 · Walk-to-walk distances

```r
d <- get_dists(walks_gw, graph = gg,
               pix.size = 1e5, readL = 3e6)
# returns a data.table of pairwise (edit, hic, longread) distances
```

Three metrics quantify "closeness" between karyotypes:

| function           | distance                                  | what it measures                                            |
|--------------------|-------------------------------------------|-------------------------------------------------------------|
| `edit_dist_cpp`    | walk-set edit distance (bp)               | Hungarian assignment over Needleman-Wunsch alignment costs  |
| `hic_kl`           | Hi-C KL                                   | NB divergence between simulated Hi-C maps                   |
| `longread_kl`      | long-read KL                              | KL between pseudo-read distributions at a given read length |

![distance matrices](man/figures/fig_distance_matrices.png)

For this graph the three metrics largely agree on the qualitative ordering;
edit distance picks up segment-by-segment differences in walk structure, while
Hi-C KL picks up quantitative differences in interchromosomal contact patterns.

---

### 5 · Fit karyotypes against observed Hi-C

```r
# Treat one karyotype as ground truth, simulate noisy observation
true_walk <- walks_gw[[1]]
true_hic  <- forward_simulate(true_walk, pix.size = 1e5, depth = 20)
obs_hic   <- make_noisymap(true_hic, theta = 2)[[1]]

# Search candidate karyotypes against the observation
fit <- bestfit_search(gg, hic = obs_hic, nsample = 50,
                      depth = 20, pix.size = 1e5, topk = 5)
```

`bestfit_search` samples candidate karyotypes, forward-simulates Hi-C for each,
and ranks by negative-binomial NLL via `compmaps` against the observation.

![NLL ranking](man/figures/fig_nll_ranking.png)

The data confidently localizes truth into the top tier of the 8-karyotype
space; depending on the noise realization, truth lands at rank 1 or 2, with a
large NLL gap to the bottom half. The rank-1 fit's predicted Hi-C is nearly
indistinguishable from the ground-truth's:

![best fit vs truth](man/figures/fig_bestfit_vs_truth.png)

For simple and long-range SV events (e.g. multi-way translocations) the karyotype space is small
and Hi-C is sufficient to localize truth tightly. For more complex and local alterations
that duplicate segments (e.g. BFB amplicons, ecDNA) the karyotype space will be exponentially larger and
`bestfit_search` may be becomes a localiser over a wider neighborhood of similar
karyotypes rather than a point estimator — the same pipeline scales to both regimes.

---

## Karyotype inference (kinf)

`bestfit_search` ranks randomly sampled karyotypes, which cannot reach the truth when the karyotype space is
huge. `infer_karyotype` instead **searches** the wirings of the graph for the best fit to the Hi-C: an exact
negative-binomial NLL with incremental updates per move, driven by parallel tempering, on a GPU (JAX). The Python
code ships in `inst/python/kinf`; it runs in its own Python environment.

### Install (once)

```r
library(fastKar)
# NVIDIA GPU nodes (needs only the driver; ~5 GB). Put envdir on a filesystem the compute nodes see.
install_kinf(envdir = '/shared/envs/fastkar_kinf', build = 'cuda12')
# or CPU only (~1 GB, ~20x slower):
install_kinf(build = 'cpu')
kinf_check()   # JAX devices + a self test that must recover a known karyotype
```

`install_kinf()` downloads a standalone Python with [uv](https://github.com/astral-sh/uv) (it runs without
environment modules, e.g. inside Slurm jobs), installs the pinned versions from
`inst/python/kinf/requirements-*.txt`, and records the interpreter in
`tools::R_user_dir('fastKar', 'config')/kinf.json`. Override with `options(fastKar.kinf.python = ...)` or the
`FASTKAR_KINF_PYTHON` environment variable.

### Run

```r
fit <- infer_karyotype(gg, hic = 'sample.mcool', pix.size = 1e5,   # or a gMatrix + depth =
                       target_region = region, chains = 16, steps = 20000,
                       topk = 5, rescale = TRUE)
fit$gw        # best karyotype, a gWalk on gg
fit$kars      # top-k distinct karyotypes visited (best first), fit$nlls their NLLs
fit$summary   # best NLL, device, timing, acceptance / swap rates, depth correction
fit$trace     # best NLL per block of proposals
```

Coarse-to-fine: `infer_karyotype(..., pix.size = 2.5e5, init = fit_1mb)` warm-starts the chains from a previous
run on the same graph and target region (its top-k karyotypes if it has them), or from a gWalk / `best_b.npy` /
`gwalks.json`; `init_random = n` keeps n chains on random starts as a control. `history_every = n` records every
chain's NLL every n proposals (`fit$history`) for convergence plots.

`gg` must have node and junction copy numbers (`loosefix()` it first). With `truth =` (a gWalk on `gg`, for
simulations) the summary also reports the NLL gap to the truth and whether it was recovered exactly.
`rescale = TRUE` corrects the depth like `bestfit_search(rescale = TRUE)`, using random wirings as candidates.

### On a Slurm cluster

```r
install_kinf(envdir = '/shared/envs/fastkar_kinf', build = 'cuda12',
             slurm = list(partition = 'gpu', gres = 'gpu:1', time = '2:00:00', account = 'mylab'))
job <- infer_karyotype(gg, hic = 'sample.mcool', pix.size = 1e5, backend = 'slurm', wait = FALSE)
fit <- kinf_collect(job)   # reports the state while pending/running, the log tail if it failed
```

With `wait = TRUE` (the default) the call blocks until the job finishes and returns the fit.

### Setting it up with an AI agent

The Claude Code skill [`.claude/skills/fastkar-kinf-install/SKILL.md`](.claude/skills/fastkar-kinf-install/SKILL.md)
walks an agent through surveying a cluster (GPUs, partitions, shared filesystems, quota), choosing the
`install_kinf()` arguments and checking the result. Claude Code finds it automatically when it runs inside a
clone of this repository; otherwise point your agent at the file on GitHub.

Convergence benchmark: `bench/kinf/` (small CCLE events with simulated truths; `run_kinf_bench.R`).

---

## Where to go next

* `R/graphstats.R` — sampling (`sample.gwalks`), karyotype hashing, walk-to-walk distances
* `R/forward.R` — forward simulation (`forward_simulate`)
* `R/comparisons.R` — likelihood and noise (`compmaps`, `make_noisymap`, `callloops`)
* `R/backwards.R` — fitting against observed Hi-C (`bestfit_search`, `prep_data_for_py`)
* `R/kinf.R`, `inst/python/kinf/` — karyotype inference by search (`install_kinf`, `infer_karyotype`)
* `R/tests.R` — synthetic graph constructors (`make_threeway_repdup`, `gimme_bfb`, `makerepdup`)
* `tests/testthat/` — regression and correctness tests
* `bench/` — benchmark harness and history (drop real `gGraph` `.rds` files in `bench/graphs/`)
