# Regression test for the circular-walk canonicalization bug fixed in
# commit 6f85062: sort_snodes / hash_karyotype_cpp used to compare
# Booth(W) against rc(Booth(W)) instead of Booth(rc(W)), leaving
# RC-equivalent circular karyotypes presented in different rotations
# of their orbit with different canonical forms.

test_that("circular karyotype: W and a rotation-of-rc(W) hash identically (C++)", {
  W1 <- c(3L, 1L, 4L, 1L, 5L)
  # rotation of -rev(W1):
  #   -rev(W1) = c(-5, -1, -4, -1, -3)
  #   rotate left by 1 -> c(-1, -4, -1, -3, -5)
  W2 <- c(-1L, -4L, -1L, -3L, -5L)

  expect_identical(
    hash_karyotype_cpp(list(W1), c(TRUE)),
    hash_karyotype_cpp(list(W2), c(TRUE))
  )
})

test_that("circular karyotype: W and a rotation-of-rc(W) hash identically (R sort_snodes)", {
  W1 <- c(3L, 1L, 4L, 1L, 5L)
  W2 <- c(-1L, -4L, -1L, -3L, -5L)
  expect_identical(
    hash_snodelist(list(W1), c(TRUE)),
    hash_snodelist(list(W2), c(TRUE))
  )
})

test_that("circular karyotypes that are NOT in the same orbit hash differently", {
  W1 <- c(1L, 2L, 3L, 4L, 5L)
  W2 <- c(1L, 2L, 3L, 5L, 4L)  # different cycle, same nodes
  expect_false(identical(
    hash_karyotype_cpp(list(W1), c(TRUE)),
    hash_karyotype_cpp(list(W2), c(TRUE))
  ))
})

test_that("R sort_snodes and hash_karyotype_cpp agree on the partition of a small orbit", {
  # Generate every rotation + RC-rotation of a circular walk; all should map
  # to the same canonical form under each hasher.
  W <- c(7L, 2L, 9L, 4L, 1L)
  rotations <- function(v) lapply(seq_along(v), function(i) v[c(i:length(v), seq_len(i - 1))])
  W_rc <- as.integer(-rev(W))
  all_forms <- c(rotations(W), rotations(W_rc))

  h_cpp <- vapply(all_forms, function(v) hash_karyotype_cpp(list(v), c(TRUE)),
                  character(1))
  h_r   <- vapply(all_forms, function(v) hash_snodelist(list(v),     c(TRUE)),
                  character(1))
  expect_equal(length(unique(h_cpp)), 1L)
  expect_equal(length(unique(h_r)),   1L)
})
