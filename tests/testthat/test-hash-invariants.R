# Properties the canonical karyotype hash must satisfy.
# These don't depend on gGnome / loosefix / any graph code — just exercise
# hash_karyotype_cpp on hand-constructed walks.

test_that("linear walk hashes equal to its reverse complement", {
  W <- c(1L, 2L, 3L, 4L)
  W_rc <- as.integer(-rev(W))
  expect_identical(
    hash_karyotype_cpp(list(W),    c(FALSE)),
    hash_karyotype_cpp(list(W_rc), c(FALSE))
  )
})

test_that("multiset reorder leaves hash unchanged", {
  W1 <- c(1L, 2L, 3L)
  W2 <- c(4L, 5L)
  W3 <- c(6L, 7L, 8L)
  h_a <- hash_karyotype_cpp(list(W1, W2, W3), c(FALSE, FALSE, FALSE))
  h_b <- hash_karyotype_cpp(list(W3, W1, W2), c(FALSE, FALSE, FALSE))
  h_c <- hash_karyotype_cpp(list(W2, W3, W1), c(FALSE, FALSE, FALSE))
  expect_identical(h_a, h_b)
  expect_identical(h_a, h_c)
})

test_that("circular walk hashes equal to all its rotations", {
  W    <- c(1L, 2L, 3L, 4L)
  rot1 <- c(2L, 3L, 4L, 1L)
  rot2 <- c(3L, 4L, 1L, 2L)
  rot3 <- c(4L, 1L, 2L, 3L)
  h <- hash_karyotype_cpp(list(W), c(TRUE))
  expect_identical(h, hash_karyotype_cpp(list(rot1), c(TRUE)))
  expect_identical(h, hash_karyotype_cpp(list(rot2), c(TRUE)))
  expect_identical(h, hash_karyotype_cpp(list(rot3), c(TRUE)))
})

test_that("circular walk hashes equal to all rotations of its RC", {
  W <- c(1L, 2L, 3L, 4L)
  W_rc <- as.integer(-rev(W))
  rc_rot1 <- W_rc[c(2:length(W_rc), 1L)]
  rc_rot2 <- W_rc[c(3:length(W_rc), 1L, 2L)]
  h <- hash_karyotype_cpp(list(W), c(TRUE))
  expect_identical(h, hash_karyotype_cpp(list(W_rc),    c(TRUE)))
  expect_identical(h, hash_karyotype_cpp(list(rc_rot1), c(TRUE)))
  expect_identical(h, hash_karyotype_cpp(list(rc_rot2), c(TRUE)))
})

test_that("linear and circular versions of the same vector are NOT equivalent", {
  W <- c(1L, 2L, 3L)
  expect_false(identical(
    hash_karyotype_cpp(list(W), c(FALSE)),
    hash_karyotype_cpp(list(W), c(TRUE))
  ))
})

test_that("two structurally different karyotypes get different hashes", {
  W1 <- c(1L, 2L, 3L)
  W2 <- c(1L, 3L, 2L)  # different sequence, same elements
  expect_false(identical(
    hash_karyotype_cpp(list(W1), c(FALSE)),
    hash_karyotype_cpp(list(W2), c(FALSE))
  ))
})

test_that("empty karyotype hashes consistently", {
  expect_identical(
    hash_karyotype_cpp(list(),                logical(0)),
    hash_karyotype_cpp(list(),                logical(0))
  )
})
