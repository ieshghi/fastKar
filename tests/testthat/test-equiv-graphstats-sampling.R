## Differential tests (vs the oracle) for the walk-space samplers in R/graphstats.R:
## gg.to.wiring, traverse_graph_cpp, local.sampling, sample.neighborhood, markov.gwalk,
## multi.markov, and sample.gwalks (unchanged; sanity check of the shared wiring).
## Slow cases (oracle traversal on the high-copy full CCLE graphs) need FASTKAR_SLOW_TESTS=1.

SLOW = Sys.getenv('FASTKAR_SLOW_TESTS') == '1'
## fixtures cheap enough for the oracle's quadratic traversal in the default run
FAST_FIX = c('toy_bfb', 'toy_dtrans', 'toy_repdup', 'toy_threeway', 'toy_tiny',
             'ccle_pyrgo', 'ccle_pyrgo_local', 'ccle_chromothripsis', 'ccle_chromothripsis_local')
sampling_fixtures = function() if (SLOW) setdiff(fixture_names(), 'toy_circ') else FAST_FIX

## walk lists that may hold gWalks: compare what defines them
walk_key = function(w) {
    if (inherits(w, 'gWalk')) return(list(snode.id = unname(w$snode.id), circular = unname(w$circular), hash = w$hash))
    w
}
walks_key = function(x) {
    if (inherits(x, 'gWalk') || (is.list(x) && !is.null(x$snode.id))) return(walk_key(x))
    if (is.list(x)) return(lapply(x, walks_key))
    x
}

test_that('gg.to.wiring is identical to the oracle on every fixture', {
    skip_if_no_oracle()
    for (nm in fixture_names()) {
        fx = load_fixture(nm)
        old = tryCatch(ref('gg.to.wiring')(fx$gg), error = function(e) e)
        new = tryCatch(opt('gg.to.wiring')(fx$gg), error = function(e) e)
        if (inherits(old, 'error')) {
            expect_s3_class(new, 'error')
            next
        }
        expect_identical(new$internal.edges, old$internal.edges, info = nm)
        expect_identical(attributes(new$internal.edges), attributes(old$internal.edges), info = nm)
        expect_identical(new$loose.ends, old$loose.ends, info = nm)
        expect_identical(new$gg$nodes$dt, old$gg$nodes$dt, info = nm)
        expect_identical(new$gg$edges$dt, old$gg$edges$dt, info = nm)
    }
})

test_that('gg.to.wiring matches the oracle on a graph without loose CN columns (loosefix path)', {
    skip_if_no_oracle()
    set.seed(3)
    gg = ref('gimme_bfb')(1e5, 6, cycles = 3)$graph
    expect_identical(opt('gg.to.wiring')(gg)$internal.edges, ref('gg.to.wiring')(gg)$internal.edges)
    gg = ref('makerepdup')(1e6, 'trans')$graph
    expect_identical(opt('gg.to.wiring')(gg)[1:2], ref('gg.to.wiring')(gg)[1:2])
})

test_that('traverse_graph_cpp is identical to the oracle on random wirings', {
    skip_if_no_oracle()
    for (nm in sampling_fixtures()) {
        fx = load_fixture(nm)
        w = ref('gg.to.wiring')(fx$gg)
        ie = data.table::copy(w$internal.edges)
        nperm = if (nrow(ie) > 1000) 2 else 15
        set.seed(11)
        for (k in seq_len(nperm)) {
            ie[, right := if (.N > 1) sample(right) else right, by = n]
            old = ref('traverse_graph_cpp')(ie, w$loose.ends)
            expect_identical(opt('traverse_graph_cpp')(ie, w$loose.ends), old, info = nm)
            expect_identical(opt('traverse_graph_fast_cpp')(ie$n, ie$left, ie$right, ie$id, w$loose.ends), old, info = nm)
        }
        ## rows out of order and ids that are not 1..n: visited is keyed on id, lookups in row order
        sh = ie[sample(.N)]
        expect_identical(opt('traverse_graph_cpp')(sh, w$loose.ends), ref('traverse_graph_cpp')(sh, w$loose.ends), info = nm)
        sh[, id := id * 7L + 3L]
        expect_identical(opt('traverse_graph_cpp')(sh, w$loose.ends), ref('traverse_graph_cpp')(sh, w$loose.ends), info = nm)
        ## duplicated ids (a visit marks every row sharing the id)
        dup = data.table::copy(ie)[, id := (id + 1L) %/% 2L]
        expect_identical(opt('traverse_graph_cpp')(dup, w$loose.ends), ref('traverse_graph_cpp')(dup, w$loose.ends), info = nm)
        ## a subset of the loose ends, with a value that matches no row
        le = c(sample(w$loose.ends, ceiling(length(w$loose.ends) / 2)), 10L^7)
        expect_identical(opt('traverse_graph_cpp')(ie, le), ref('traverse_graph_cpp')(ie, le), info = nm)
    }
})

test_that('local.sampling is identical to the oracle (seeds, sizes, flags, frozen nodes, starters)', {
    skip_if_no_oracle()
    for (nm in sampling_fixtures()) {
        fx = load_fixture(nm)
        big = nrow(ref('gg.to.wiring')(fx$gg)$internal.edges) > 1000
        nsteps = if (big) 2 else 4
        nwalk = if (big) 2 else 3
        run = function(pkg, ...) { f = if (pkg == 'ref') ref('local.sampling') else opt('local.sampling'); set.seed(7); f(fx$gg, ...) }
        old = run('ref', nsteps = nsteps, nwalk = nwalk)
        expect_identical(run('opt', nsteps = nsteps, nwalk = nwalk), old, info = nm)
        old = run('ref', nsteps = nsteps, nwalk = nwalk, onlyhash = TRUE, return.edges = TRUE)
        new = run('opt', nsteps = nsteps, nwalk = nwalk, onlyhash = TRUE, return.edges = TRUE)
        expect_identical(new, old, info = nm)
        expect_identical(attributes(new$starter_edges), attributes(old$starter_edges), info = nm)
        ## RNG stream left in the same state
        set.seed(7); ref('local.sampling')(fx$gg, nsteps = nsteps, nwalk = nwalk); r1 = runif(1)
        set.seed(7); opt('local.sampling')(fx$gg, nsteps = nsteps, nwalk = nwalk); r2 = runif(1)
        expect_identical(r2, r1, info = nm)
        if (big) next
        expect_identical(run('opt', nsteps = 1, nwalk = 2), run('ref', nsteps = 1, nwalk = 2), info = nm)
        expect_identical(run('opt', nsteps = 6, nwalk = 1, return.edges = TRUE),
                         run('ref', nsteps = 6, nwalk = 1, return.edges = TRUE), info = nm)
        ## frozen nodes: half of the swappable nodes, then all of them (no pivots left)
        ie = ref('gg.to.wiring')(fx$gg)$internal.edges
        piv = unique(ie[cn > 1]$n)
        if (length(piv)) {
            fr = piv[seq_len(ceiling(length(piv) / 2))]
            expect_identical(run('opt', nsteps = 4, nwalk = 2, frozen.nodes = fr),
                             run('ref', nsteps = 4, nwalk = 2, frozen.nodes = fr), info = nm)
            expect_identical(run('opt', nsteps = 3, nwalk = 2, frozen.nodes = piv, onlyhash = TRUE),
                             run('ref', nsteps = 3, nwalk = 2, frozen.nodes = piv, onlyhash = TRUE), info = nm)
        }
        ## starter_edges from a previous run (keyed), and an unkeyed copy
        st = run('ref', nsteps = 2, nwalk = 1, return.edges = TRUE)$starter_edges
        for (s0 in list(st, data.table::setkey(data.table::copy(st), NULL))) {
            a = data.table::copy(s0); b = data.table::copy(s0)
            old = run('ref', nsteps = 3, nwalk = 2, starter_edges = a, return.edges = TRUE)
            new = run('opt', nsteps = 3, nwalk = 2, starter_edges = b, return.edges = TRUE)
            expect_identical(new, old, info = nm)
            expect_identical(attributes(b), attributes(a), info = nm)
        }
        ## precomputed wiring: same result, wiring left untouched
        wiring = opt('gg.to.wiring')(fx$gg)
        w0 = data.table::copy(wiring$internal.edges)
        expect_identical(run('opt', nsteps = 3, nwalk = 2, wiring = wiring), run('ref', nsteps = 3, nwalk = 2), info = nm)
        expect_identical(wiring$internal.edges, w0, info = nm)
    }
})

test_that('local.sampling with return.gw matches the oracle on toys', {
    skip_if_no_oracle()
    for (nm in c('toy_bfb', 'toy_repdup', 'toy_threeway')) {
        fx = load_fixture(nm)
        set.seed(5); old = ref('local.sampling')(fx$gg, nsteps = 3, nwalk = 2, return.gw = TRUE, return.edges = TRUE)
        set.seed(5); new = opt('local.sampling')(fx$gg, nsteps = 3, nwalk = 2, return.gw = TRUE, return.edges = TRUE)
        expect_identical(walks_key(new$walks), walks_key(old$walks), info = nm)
        expect_identical(new$starter_edges, old$starter_edges, info = nm)
        set.seed(5); old = ref('local.sampling')(fx$gg, nsteps = 3, nwalk = 2, return.gw = TRUE, onlyhash = TRUE)
        set.seed(5); new = opt('local.sampling')(fx$gg, nsteps = 3, nwalk = 2, return.gw = TRUE, onlyhash = TRUE)
        expect_identical(new, old, info = nm)
    }
})

test_that('sample.neighborhood is identical to the oracle', {
    skip_if_no_oracle()
    for (nm in sampling_fixtures()) {
        fx = load_fixture(nm)
        big = nrow(ref('gg.to.wiring')(fx$gg)$internal.edges) > 1000
        nn = if (big) 2 else 6
        set.seed(9); old = ref('sample.neighborhood')(fx$gg, nn, return.gw = FALSE, return.edges = TRUE)
        set.seed(9); new = opt('sample.neighborhood')(fx$gg, nn, return.gw = FALSE, return.edges = TRUE)
        expect_identical(new, old, info = nm)
        if (big) next
        st = data.table::copy(old$starter_edges)
        set.seed(2); old = ref('sample.neighborhood')(fx$gg, nn, starter_edges = data.table::copy(st), return.gw = FALSE)
        set.seed(2); new = opt('sample.neighborhood')(fx$gg, nn, starter_edges = data.table::copy(st), return.gw = FALSE,
                                                     wiring = opt('gg.to.wiring')(fx$gg))
        expect_identical(new, old, info = nm)
        if (grepl('^toy_', nm)) {
            set.seed(4); old = ref('sample.neighborhood')(fx$gg, 3)
            set.seed(4); new = opt('sample.neighborhood')(fx$gg, 3)
            expect_identical(walks_key(new), walks_key(old), info = nm)
        }
    }
})

test_that('markov.gwalk and multi.markov are identical to the oracle', {
    skip_if_no_oracle()
    for (nm in sampling_fixtures()) {
        fx = load_fixture(nm)
        ie = ref('gg.to.wiring')(fx$gg)$internal.edges
        if (!nrow(ie[cn > 1])) next  # reference errors (sample() of an empty vector); see below
        big = nrow(ie) > 1000
        len = if (big) 2 else 4
        for (sa in c(FALSE, TRUE)) {
            old = ref('markov.gwalk')(fx$gg, len, self.avoid = sa, attempts = 3, seed = 13)
            new = opt('markov.gwalk')(fx$gg, len, self.avoid = sa, attempts = 3, seed = 13)
            expect_identical(new, old, info = paste(nm, sa))
        }
        if (big) next
        expect_identical(opt('markov.gwalk')(fx$gg, 3, attempts = 1, seed = 2), ref('markov.gwalk')(fx$gg, 3, attempts = 1, seed = 2), info = nm)
        expect_identical(opt('multi.markov')(fx$gg, 3, 3, seed = 21), ref('multi.markov')(fx$gg, 3, 3, seed = 21), info = nm)
        if (grepl('^toy_', nm)) {
            expect_identical(walks_key(opt('markov.gwalk')(fx$gg, 3, return.gw = TRUE, seed = 3)),
                             walks_key(ref('markov.gwalk')(fx$gg, 3, return.gw = TRUE, seed = 3)), info = nm)
        }
    }
    ## no swappable node: both error in sample()
    fx = load_fixture('toy_dtrans')
    if (!nrow(ref('gg.to.wiring')(fx$gg)$internal.edges[cn > 1])) {
        expect_error(ref('markov.gwalk')(fx$gg, 3, seed = 1))
        expect_error(opt('markov.gwalk')(fx$gg, 3, seed = 1))
    }
})

test_that('sample.gwalks (unchanged) still matches the oracle for a fixed seed', {
    skip_if_no_oracle()
    for (nm in c('toy_bfb', 'toy_repdup', 'toy_threeway', 'ccle_chromothripsis_local', 'ccle_pyrgo')) {
        fx = load_fixture(nm)
        for (oh in c(FALSE, TRUE)) {
            set.seed(17); old = ref('sample.gwalks')(fx$gg, 20, verbose = FALSE, return.gw = FALSE, onlyhash = oh)
            set.seed(17); new = opt('sample.gwalks')(fx$gg, 20, verbose = FALSE, return.gw = FALSE, onlyhash = oh)
            expect_identical(new, old, info = paste(nm, oh))
        }
    }
})
