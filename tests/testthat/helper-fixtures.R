## Shared helpers for the fixture-based property tests (fixtures built by tests/make_fixtures.R).

FIXTURE_DIR = testthat::test_path('fixtures')

## internal function `name` of the package under test
opt = function(name) get(name, envir = asNamespace('fastKar'))

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
