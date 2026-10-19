## ecDNA caller helpers (R/ecfinder.R).

## boil_greedy() peels circles off a working copy of the graph; the graph passed in must be left as is.
## A graph read from rds shares its edge table with its $copy, so peeling on gg$copy used to strip
## edge copy number from the caller's graph (and from everything fastKar_ecfinder sampled from it next).

test_that('boil_greedy leaves an rds-loaded graph unchanged', {
    fx = load_fixture('ccle_dm') # the full graph: its $copy shares edges, and N = 200 finds circles to peel
    gg = fx$gg
    state = function(g) list(node.cn = g$nodes$dt$cn, edge.cn = g$edges$dt$cn)
    before = state(gg)
    set.seed(1)
    solns = boil_greedy(gg = gg, ft = fx$ft, N = 200, n_solns_try = 3)
    expect_true(length(solns) > 0)
    expect_identical(state(gg), before)
})

## squeeze() scores raw karyotypes (sample.gwalks(return.gw = FALSE)) and builds gWalks only for the winners,
## so the raw entropy must equal amp.entropy() on the gWalk, bit for bit (ties in the ranking depend on it)
test_that('amp.entropy.raw matches amp.entropy on the fixture karyotypes', {
    for (nm in c('ccle_dm', 'ccle_pyrgo', 'toy_bfb')) {
        fx = load_fixture(nm)
        amp = unique((fx$gg$nodes$gr %&% fx$ft)$node.id)
        for (i in seq_along(fx$kars)) {
            expect_identical(amp.entropy.raw(fx$kars[[i]], amp), amp.entropy(fx$gw[[i]], amp))
        }
    }
})
