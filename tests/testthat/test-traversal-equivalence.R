# Verifies that the three traversal implementations:
#   * traverse_graph_cpp        (baseline)
#   * traverse_graph_v2_cpp     (O(V) refactor)
#   * traverse_graph_v2_batch_cpp (batched version)
# produce the same karyotype-partition for a fixed set of permutations.

shuffle_in_groups <- function(internal_edges, n_groups) {
  out <- internal_edges$right
  for (g in n_groups) {
    if (length(g) > 1) out[g] <- sample(out[g])
  }
  out
}

test_that("v1, v2, batch produce the same partition on a small BFB", {
  skip_if_no_loosefix()

  set.seed(42)
  gg <- gimme_bfb(1e6, 8, cycles = 5)$graph
  wiring <- gg.to.wiring(gg)
  ie <- wiring$internal.edges
  le <- wiring$loose.ends

  n_groups <- split(seq_len(nrow(ie)), ie$n)
  set.seed(1)
  N <- 50
  perms <- lapply(seq_len(N), function(i) shuffle_in_groups(ie, n_groups))

  walks_v1    <- lapply(perms, function(p) traverse_graph_cpp(    ie[, right := p], le))
  walks_v2    <- lapply(perms, function(p) traverse_graph_v2_cpp( ie[, right := p], le))
  walks_batch <- traverse_graph_v2_batch_cpp(ie, perms, le)

  hash <- function(w) hash_karyotype_cpp(w$snode.id, w$circular)
  h_v1    <- vapply(walks_v1,    hash, character(1))
  h_v2    <- vapply(walks_v2,    hash, character(1))
  h_batch <- vapply(walks_batch, hash, character(1))

  # Partition fingerprint: position of each element's first occurrence.
  # Two hash vectors that induce the same equivalence relation produce
  # identical match(h, h), regardless of what their hash strings spell.
  expect_identical(match(h_v1, h_v1), match(h_v2, h_v2))
  expect_identical(match(h_v1, h_v1), match(h_batch, h_batch))
})

test_that("v1, v2, batch agree on a small reciprocal-dup graph", {
  skip_if_no_loosefix()

  gg <- makerepdup(1e6, geometry = "cis1")$graph
  wiring <- gg.to.wiring(gg)
  ie <- wiring$internal.edges
  le <- wiring$loose.ends

  n_groups <- split(seq_len(nrow(ie)), ie$n)
  set.seed(7)
  N <- 50
  perms <- lapply(seq_len(N), function(i) shuffle_in_groups(ie, n_groups))

  walks_v1    <- lapply(perms, function(p) traverse_graph_cpp(    ie[, right := p], le))
  walks_v2    <- lapply(perms, function(p) traverse_graph_v2_cpp( ie[, right := p], le))
  walks_batch <- traverse_graph_v2_batch_cpp(ie, perms, le)

  hash <- function(w) hash_karyotype_cpp(w$snode.id, w$circular)
  h_v1    <- vapply(walks_v1,    hash, character(1))
  h_v2    <- vapply(walks_v2,    hash, character(1))
  h_batch <- vapply(walks_batch, hash, character(1))

  expect_identical(match(h_v1, h_v1), match(h_v2, h_v2))
  expect_identical(match(h_v1, h_v1), match(h_batch, h_batch))
})
