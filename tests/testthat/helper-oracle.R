## Shared helpers for the fastKarOptim equivalence tests.
## The oracle is the reference fastKar pinned at commit 0a5497e (plus the then-untracked
## R/identifiability.R), checked out as a detached worktree at ORACLE_PATH. It is loaded as
## namespace 'fastKar' next to the package under test ('fastKarOptim').

ORACLE_PATH = Sys.getenv('FASTKAR_ORACLE', path.expand('~/git/fastKar_oracle'))
FIXTURE_DIR = testthat::test_path('fixtures')
EQUIV_TOL = 1e-12

skip_if_no_oracle = function() {
    if (!dir.exists(ORACLE_PATH)) testthat::skip(paste('oracle not found at', ORACLE_PATH))
    if (!'fastKar' %in% loadedNamespaces()) {
        suppressMessages(pkgload::load_all(ORACLE_PATH, export_all = FALSE, attach = FALSE, quiet = TRUE))
    }
    invisible(TRUE)
}

## reference implementation of function `name`
ref = function(name) get(name, envir = asNamespace('fastKar'))
## implementation under test
opt = function(name) get(name, envir = asNamespace('fastKarOptim'))

fixture_names = function(toy_only = FALSE) {
    nm = sub('\\.rds$', '', list.files(FIXTURE_DIR, pattern = '\\.rds$'))
    if (toy_only) nm = grep('^toy_', nm, value = TRUE)
    nm
}

## fixture list(name, gg, ft, kars, note); kars rebuilt as gWalks on the fixture graph
load_fixture = function(name) {
    fx = readRDS(file.path(FIXTURE_DIR, paste0(name, '.rds')))
    fx$gw = lapply(fx$kars, function(k) {
        gw = gGnome::gW(graph = fx$gg, snode.id = k$snode.id, circular = k$circular)
        if (!is.null(k$cn)) gw$set(cn = k$cn)
        gw
    })
    fx
}

## pull a long-format pixel table out of a gMatrix, a list(gr, dat) or a data.table
as_pixdt = function(x) {
    if (inherits(x, 'gMatrix')) return(data.table::as.data.table(x$dat))
    if (data.table::is.data.table(x)) return(x)
    if (is.list(x) && !is.null(x$dat)) return(data.table::as.data.table(x$dat))
    stop('cannot extract a pixel table')
}

## structural columns identical, numeric columns equal to a relative tolerance
expect_pix_equiv = function(new, old, tol = EQUIV_TOL, cols = c('i', 'j', 'id'), info = NULL) {
    a = data.table::copy(as_pixdt(new)); b = data.table::copy(as_pixdt(old))
    cols = intersect(cols, intersect(names(a), names(b)))
    data.table::setorderv(a, cols); data.table::setorderv(b, cols)
    testthat::expect_identical(nrow(a), nrow(b), info = info)
    for (cc in cols) testthat::expect_identical(as.integer(a[[cc]]), as.integer(b[[cc]]), info = paste(info, cc))
    testthat::expect_equal(a$value, b$value, tolerance = tol, info = paste(info, 'value'))
}

## same, for a pair of gMatrix objects: also compare the bins
expect_gm_equiv = function(new, old, tol = EQUIV_TOL, info = NULL) {
    testthat::expect_identical(as.character(GenomicRanges::seqnames(new$gr)), as.character(GenomicRanges::seqnames(old$gr)), info = info)
    testthat::expect_identical(GenomicRanges::start(new$gr), GenomicRanges::start(old$gr), info = info)
    testthat::expect_identical(GenomicRanges::end(new$gr), GenomicRanges::end(old$gr), info = info)
    expect_pix_equiv(new, old, tol = tol, info = info)
}

## run f() twice under the same seed, once per implementation
with_seed = function(seed, f) { set.seed(seed); f() }
