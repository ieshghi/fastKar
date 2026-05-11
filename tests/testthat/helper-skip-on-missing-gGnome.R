skip_if_no_loosefix <- function() {
  if (!requireNamespace("gGnome", quietly = TRUE)) {
    testthat::skip("gGnome not installed")
  }
  if (!exists("loosefix", envir = asNamespace("gGnome"), inherits = FALSE)) {
    testthat::skip("gGnome::loosefix missing -- need the fork that exports it")
  }
}
