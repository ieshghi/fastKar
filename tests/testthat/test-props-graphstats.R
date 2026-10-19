## Property tests for R/graphstats.R, independent of the oracle:
##  * local.sampling / sample.neighborhood / markov.gwalk walk sets are decompositions of the
##    graph: every node is visited exactly CN times;
##  * edit_dist_cpp is symmetric and zero on identical inputs;
##  * reads_fromwalk's total multiplicity equals the number of read starts
##    (2 x sum over walks of total + readL - 1 for a linear walk, total for a circular one).

node_visits_match_cn = function(walks, gg, info = NULL) {
    cn = gg$nodes$dt[, .(node.id, cn)]
    for (w in walks) {
        v = table(factor(abs(unlist(w$snode.id)), levels = cn$node.id))
        expect_identical(as.numeric(v), as.numeric(cn$cn), info = info)
    }
}

PROP_FIX = c('toy_bfb', 'toy_repdup', 'toy_threeway', 'toy_tiny', 'ccle_pyrgo_local', 'ccle_dm_local',
             'ccle_chromothripsis_local', 'ccle_tyfonas_local', 'ccle_pyrgo', 'ccle_tyfonas')

test_that('sampled walk sets visit every node CN times', {
    for (nm in PROP_FIX) {
        fx = load_fixture(nm)
        gg = opt('gg.to.wiring')(fx$gg)$gg
        set.seed(1)
        ls = opt('local.sampling')(fx$gg, nsteps = 5, nwalk = 3)
        for (chain in ls) node_visits_match_cn(chain, gg, info = nm)
        set.seed(2)
        sn = opt('sample.neighborhood')(fx$gg, 5, return.gw = FALSE)
        node_visits_match_cn(c(list(sn$orig_walk), sn$neighbor_walks), gg, info = nm)
        ie = opt('gg.to.wiring')(fx$gg)$internal.edges
        if (nrow(ie[cn > 1])) node_visits_match_cn(opt('markov.gwalk')(fx$gg, 4, seed = 3), gg, info = nm)
        ## every walk set hashes consistently with its walks
        for (w in ls[[1]]) expect_identical(w$hash, hash_karyotype_cpp(w$snode.id, w$circular), info = nm)
    }
})

test_that('local.sampling chains move by single swaps (consecutive wirings differ in two rows)', {
    fx = load_fixture('ccle_dm_local')
    w = opt('gg.to.wiring')(fx$gg)
    st = opt('.wiring_state')(w$internal.edges, w$loose.ends)
    set.seed(4)
    r0 = st$right
    piv = setdiff(st$n[which(st$cn > 1)], NULL)
    rows = lapply(piv, function(p) opt('.wiring_rows')(st, p))
    for (k in 1:50) {
        r1 = opt('.wiring_swap')(r0, rows[[sample.int(length(rows), 1)]])
        expect_identical(sum(r1 != r0), 2L)
        expect_identical(sort(r1), sort(r0))
        r0 = r1
    }
})

test_that('edit_dist_cpp is symmetric and zero on identical inputs', {
    for (nm in c('toy_bfb', 'toy_repdup', 'toy_threeway', 'ccle_pyrgo_local', 'ccle_dm_local', 'ccle_chromothripsis')) {
        fx = load_fixture(nm)
        gws = fx$gw
        nid = fx$gg$nodes$dt$node.id
        use = nid[seq(1, length(nid), by = 2)]
        for (k in seq_along(gws)[1:min(3, length(gws))]) {
            expect_equal(opt('edit_dist_cpp')(gws[[k]], gws[[k]]), 0, info = nm)
            expect_equal(opt('edit_dist_cpp')(gws[[k]], gws[[k]], use.nodes = use, constpenalty = TRUE), 0, info = nm)
            if (k > 1) {
                expect_equal(opt('edit_dist_cpp')(gws[[1]], gws[[k]]), opt('edit_dist_cpp')(gws[[k]], gws[[1]]),
                             tolerance = 1e-12, info = nm)
                expect_equal(opt('edit_dist_cpp')(gws[[1]], gws[[k]], use.nodes = use),
                             opt('edit_dist_cpp')(gws[[k]], gws[[1]], use.nodes = use), tolerance = 1e-12, info = nm)
                expect_gte(opt('edit_dist_cpp')(gws[[1]], gws[[k]]), 0)
            }
        }
    }
})

test_that('reads_fromwalk total multiplicity equals the number of read starts', {
    for (nm in c('toy_bfb', 'toy_circ', 'toy_tiny', 'toy_repdup', 'ccle_pyrgo_local', 'ccle_dm_local', 'ccle_chromothripsis')) {
        fx = load_fixture(nm)
        gw = fx$gw[[1]]
        wd = fx$gg$nodes$dt$width
        tot = vapply(gw$snode.id, function(v) sum(as.numeric(wd[abs(v)])), 1)
        for (L in c(1, 1e3, 1e5, 1e7)) {
            expected = 2 * sum(ifelse(gw$circular, tot, tot + L - 1))
            r = opt('reads_fromwalk')(gw, L)
            expect_identical(sum(r$nums), expected, info = paste(nm, L))
            expect_false(anyDuplicated(r$words) > 0)
            ## minsize / use.nodes only merge words, never change the mass
            r2 = opt('reads_fromwalk')(gw, L, minsize = stats::median(wd), use.nodes = fx$gg$nodes$dt$node.id[-1])
            expect_identical(sum(r2$nums), expected, info = paste(nm, L, 'filtered'))
        }
    }
})
