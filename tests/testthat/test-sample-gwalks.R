# End-to-end: sample.gwalks(legacy = TRUE) and sample.gwalks(legacy = FALSE)
# must produce the same set of unique karyotypes for a fixed seed.

test_that("sample.gwalks legacy vs new produce identical karyotype sets on small BFB", {
  skip_if_no_loosefix()

  set.seed(42); gg <- gimme_bfb(1e6, 8, cycles = 5)$graph

  set.seed(1); walks_l <- sample.gwalks(gg, N = 50, return.gw = FALSE,
                                         remove.dups = TRUE, verbose = FALSE,
                                         mc.cores = 1, legacy = TRUE)
  set.seed(1); walks_n <- sample.gwalks(gg, N = 50, return.gw = FALSE,
                                         remove.dups = TRUE, verbose = FALSE,
                                         mc.cores = 1, legacy = FALSE)

  hash <- function(w) hash_karyotype_cpp(w$snode.id, w$circular)
  expect_identical(
    sort(vapply(walks_l, hash, character(1))),
    sort(vapply(walks_n, hash, character(1)))
  )
})

test_that("sample.gwalks legacy vs new produce identical karyotype sets on repdup", {
  skip_if_no_loosefix()

  gg <- makerepdup(1e6, geometry = "trans")$graph

  set.seed(1); walks_l <- sample.gwalks(gg, N = 100, return.gw = FALSE,
                                         remove.dups = TRUE, verbose = FALSE,
                                         mc.cores = 1, legacy = TRUE)
  set.seed(1); walks_n <- sample.gwalks(gg, N = 100, return.gw = FALSE,
                                         remove.dups = TRUE, verbose = FALSE,
                                         mc.cores = 1, legacy = FALSE)

  hash <- function(w) hash_karyotype_cpp(w$snode.id, w$circular)
  expect_identical(
    sort(vapply(walks_l, hash, character(1))),
    sort(vapply(walks_n, hash, character(1)))
  )
})

test_that("walks reported by sample.gwalks have plausible structure", {
  skip_if_no_loosefix()

  set.seed(42); gg <- gimme_bfb(1e6, 8, cycles = 5)$graph
  V <- nrow(gg.to.wiring(gg)$internal.edges)

  set.seed(1); walks <- sample.gwalks(gg, N = 10, return.gw = FALSE,
                                      remove.dups = FALSE, verbose = FALSE,
                                      mc.cores = 1, legacy = FALSE)
  expect_equal(length(walks), 10L)
  for (w in walks) {
    expect_true(is.list(w$snode.id))
    expect_true(is.logical(w$circular))
    expect_equal(length(w$snode.id), length(w$circular))
    # The walk-set should cover all V node-copies (junction balance)
    expect_equal(sum(vapply(w$snode.id, length, integer(1))), V)
  }
})
