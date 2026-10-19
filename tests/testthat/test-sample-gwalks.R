# End-to-end checks of sample.gwalks(). The legacy implementation it was once
# compared against has been removed, so these check that the karyotype set does
# not depend on how permutations are batched (chunksize) and that every sampled
# walk set is structurally valid. Until 2026-10-19, N >= 2 * chunksize that chunksize
# did not divide drew only chunksize * (N %/% chunksize) permutations.

sample_hashes <- function(gg, N, chunksize) {
  set.seed(1)
  walks <- sample.gwalks(gg, N = N, return.gw = FALSE, remove.dups = TRUE,
                         verbose = FALSE, mc.cores = 1, chunksize = chunksize)
  sort(unname(vapply(walks, function(w) hash_karyotype_cpp(w$snode.id, w$circular), character(1))))
}

test_that("sample.gwalks karyotype set does not depend on chunksize on small BFB", {
  skip_if_no_loosefix()

  set.seed(42); gg <- gimme_bfb(1e6, 8, cycles = 5)$graph
  expect_identical(sample_hashes(gg, 50, chunksize = 1e3), sample_hashes(gg, 50, chunksize = 10))
})

test_that("sample.gwalks karyotype set does not depend on chunksize on repdup", {
  skip_if_no_loosefix()

  gg <- makerepdup(1e6, geometry = "trans")$graph
  expect_identical(sample_hashes(gg, 100, chunksize = 1e3), sample_hashes(gg, 100, chunksize = 20))
})

test_that("sample.gwalks draws exactly N karyotypes when chunksize does not divide N", {
  skip_if_no_loosefix()

  set.seed(42); gg <- gimme_bfb(1e6, 8, cycles = 5)$graph
  set.seed(1)
  walks <- sample.gwalks(gg, N = 25, return.gw = FALSE, remove.dups = FALSE,
                         verbose = FALSE, mc.cores = 1, chunksize = 10)
  expect_length(walks, 25)
  expect_identical(sample_hashes(gg, 55, chunksize = 1e3), sample_hashes(gg, 55, chunksize = 10))
})

test_that("walks reported by sample.gwalks have plausible structure", {
  skip_if_no_loosefix()

  set.seed(42); gg <- gimme_bfb(1e6, 8, cycles = 5)$graph
  V <- nrow(gg.to.wiring(gg)$internal.edges)

  set.seed(1); walks <- sample.gwalks(gg, N = 10, return.gw = FALSE,
                                      remove.dups = FALSE, verbose = FALSE,
                                      mc.cores = 1)
  expect_equal(length(walks), 10L)
  for (w in walks) {
    expect_true(is.list(w$snode.id))
    expect_true(is.logical(w$circular))
    expect_equal(length(w$snode.id), length(w$circular))
    # The walk-set should cover all V node-copies (junction balance)
    expect_equal(sum(vapply(w$snode.id, length, integer(1))), V)
  }
})
