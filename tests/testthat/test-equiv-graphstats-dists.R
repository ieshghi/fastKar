## Differential tests (vs the oracle) for the walk distances in R/graphstats.R:
## compute_cost_matrix_cpp, split_walk_fragments, edit_dist_cpp (all modes), hic_kl,
## get_dists (pix.size > 0, readL > 0) and get_all_pair_distances.

SLOW = Sys.getenv('FASTKAR_SLOW_TESTS') == '1'

rand_walks = function(nw, max_id, max_len) {
    lapply(seq_len(nw), function(i) {
        len = sample.int(max_len, 1)
        as.integer(sample(c(-1, 1), len, TRUE) * sample.int(max_id, len, TRUE))
    })
}

test_that('compute_cost_matrix_cpp is identical to the oracle on random walk lists', {
    skip_if_no_oracle()
    set.seed(101)
    for (k in 1:40) {
        max_id = sample(c(3, 10, 50), 1)
        pen = -stats::runif(max_id) * sample(c(1, 1e5), 1)
        pen[sample.int(max_id, max_id %/% 4)] = 0
        x = rand_walks(sample(0:6, 1), max_id, 12)
        y = rand_walks(sample(0:6, 1), max_id, 12)
        if (length(x) && stats::runif(1) < 0.3) x[[1]] = integer(0)
        expect_identical(opt('compute_cost_matrix_cpp')(x, y, pen), ref('compute_cost_matrix_cpp')(x, y, pen), info = k)
    }
})

test_that('split_walk_fragments is identical to the oracle on random inputs', {
    skip_if_no_oracle()
    set.seed(202)
    for (k in 1:200) {
        max_id = sample(c(2, 5, 20), 1)
        nl = rand_walks(sample(0:8, 1), max_id, 10)
        circ = stats::runif(length(nl)) < 0.5
        pen = -as.numeric(stats::runif(max_id) < 0.7)
        li = stats::runif(max_id) < 0.6
        ri = stats::runif(max_id) < 0.6
        expect_identical(opt('split_walk_fragments')(nl, circ, pen, li, ri),
                         ref('split_walk_fragments')(nl, circ, pen, li, ri), info = k)
    }
    expect_identical(opt('split_walk_fragments')(list(), logical(0), -1, TRUE, TRUE), list())
    expect_identical(opt('split_walk_fragments')(list(1L), FALSE, 0, TRUE, TRUE), list())
})

## plain ranges (no seqinfo, unstranded) so fixtures with different seqlengths can be overlapped
bare = function(gr) GenomicRanges::GRanges(as.character(GenomicRanges::seqnames(gr)),
                                           IRanges::IRanges(GenomicRanges::start(gr), GenomicRanges::end(gr)))
in_ft = function(gr, ft) GenomicRanges::countOverlaps(bare(gr), bare(ft)) > 0
## left/right ends of each node inside the fixture footprint, named by node.id
split_ends = function(fx) {
    gr = bare(fx$gg$nodes$gr)
    l = GenomicRanges::resize(gr, 1, fix = 'start'); r = GenomicRanges::resize(gr, 1, fix = 'end')
    list(left = setNames(in_ft(l, fx$ft), fx$gg$nodes$dt$node.id),
         right = setNames(in_ft(r, fx$ft), fx$gg$nodes$dt$node.id))
}

test_that('edit_dist_cpp is identical to the oracle in every mode', {
    skip_if_no_oracle()
    nms = if (SLOW) setdiff(fixture_names(), 'toy_circ') else
        c('toy_bfb', 'toy_dtrans', 'toy_repdup', 'toy_threeway', 'toy_tiny', 'ccle_pyrgo', 'ccle_dm_local', 'ccle_chromothripsis')
    for (nm in nms) {
        fx = load_fixture(nm)
        gws = fx$gw
        if (length(gws) < 2) gws = c(gws, gws)
        nid = fx$gg$nodes$dt$node.id
        use = nid[seq(1, length(nid), by = 2)]
        nw = setNames(as.numeric(seq_along(use)) * 1000, use)
        se = split_ends(fx)
        modes = list(list(), list(constpenalty = TRUE), list(use.nodes = use), list(node.widths = nw),
                     list(split.ends = se), list(use.nodes = use, split.ends = se, node.widths = nw),
                     list(return_all = TRUE), list(return_all = TRUE, split.ends = se, constpenalty = TRUE))
        for (k in 2:min(3, length(gws))) for (md in modes) {
            info = paste(nm, k, paste(names(md), collapse = ','))
            old = do.call(ref('edit_dist_cpp'), c(list(gws[[1]], gws[[k]]), md))
            new = do.call(opt('edit_dist_cpp'), c(list(gws[[1]], gws[[k]]), md))
            expect_identical(new, old, info = info)
        }
        ## graph passed separately, raw traversal lists (as identifiability.R does)
        raw = list(snode.id = gws[[1]]$snode.id, circular = gws[[1]]$circular)
        expect_identical(opt('edit_dist_cpp')(raw, gws[[2]], graph = fx$gg, use.nodes = use),
                         ref('edit_dist_cpp')(raw, gws[[2]], graph = fx$gg, use.nodes = use), info = nm)
    }
})

## Hi-C KL is a sum of x log(x/y) terms that largely cancel, so the ~1e-14 pixel-level
## differences left by the reordered forward_simulate accumulator (F1) grow to ~1e-11 here.
KL_TOL = 1e-9
expect_dists_equiv = function(new, old, info = NULL) {
    expect_identical(names(new), names(old), info = info)
    other = setdiff(names(old), 'hic')
    expect_identical(as.data.frame(new)[other], as.data.frame(old)[other], info = info)
    expect_equal(new$hic, old$hic, tolerance = KL_TOL, info = info)
}

test_that('hic_kl matches the oracle', {
    skip_if_no_oracle()
    for (nm in c('toy_bfb', 'toy_threeway', 'ccle_pyrgo_local')) {
        fx = load_fixture(nm)
        a = fx$gw[[1]]; b = fx$gw[[2]]
        expect_equal(opt('hic_kl')(a, b, pix.size = 1e5, depth = 10), ref('hic_kl')(a, b, pix.size = 1e5, depth = 10), tolerance = KL_TOL, info = nm)
        expect_equal(opt('hic_kl')(a, b, target_region = fx$ft, pix.size = 2e5, depth = 3, theta = 5),
                     ref('hic_kl')(a, b, target_region = fx$ft, pix.size = 2e5, depth = 3, theta = 5), tolerance = KL_TOL, info = nm)
    }
})

test_that('get_dists matches the oracle (hic to KL_TOL; edit, longread exact)', {
    skip_if_no_oracle()
    cases = list(list(nm = 'toy_bfb', n = 4, pix = 1e5, readL = 1e5, use = NULL),
                 list(nm = 'toy_threeway', n = 3, pix = 2e5, readL = 3e4, use = NULL),
                 list(nm = 'ccle_pyrgo_local', n = 3, pix = 1e5, readL = 1e5, use = 'ft'),
                 list(nm = 'toy_repdup', n = 3, pix = 0, readL = 0, use = NULL))
    if (SLOW) cases = c(cases, list(list(nm = 'ccle_dm_local', n = 4, pix = 1e5, readL = 1e5, use = 'ft'),
                                    ## default-target call skipped: the first walk's footprint spans whole chromosomes and
                                    ## the oracle's rbindlist overflows (~6e10 rows) at 100 kb
                                    list(nm = 'ccle_chromothripsis', n = 3, pix = 1e5, readL = 1e5, use = 'ft', no_default = TRUE)))
    for (cs in cases) {
        fx = load_fixture(cs$nm)
        gw = fx$gw[seq_len(min(cs$n, length(fx$gw)))]
        use = if (identical(cs$use, 'ft')) fx$gg$nodes$dt$node.id[in_ft(fx$gg$nodes$gr, fx$ft)] else NULL
        old = ref('get_dists')(gw, target_region = fx$ft, pix.size = cs$pix, readL = cs$readL, use.nodes = use, depth = 10)
        new = opt('get_dists')(gw, target_region = fx$ft, pix.size = cs$pix, readL = cs$readL, use.nodes = use, depth = 10)
        expect_dists_equiv(new, old, info = cs$nm)
        ## default target (first walk's footprint) and an explicit graph
        if (cs$pix > 0 && !isTRUE(cs$no_default)) {
            old = ref('get_dists')(gw, graph = fx$gg, pix.size = cs$pix, readL = cs$readL, depth = 2)
            new = opt('get_dists')(gw, graph = fx$gg, pix.size = cs$pix, readL = cs$readL, depth = 2)
            expect_dists_equiv(new, old, info = cs$nm)
        }
    }
    ## the reads cache is scoped to the call
    expect_length(get('.reads_cache', envir = asNamespace('fastKarOptim'))$entries, 0)
})

test_that('get_all_pair_distances matches the oracle', {
    skip_if_no_oracle()
    fx = load_fixture('toy_bfb')
    gw = fx$gw[1:3]
    old = suppressMessages(ref('get_all_pair_distances')(gw, depth = 10, pix.size = 1e5))
    new = suppressMessages(opt('get_all_pair_distances')(gw, depth = 10, pix.size = 1e5))
    expect_equal(new, old, tolerance = KL_TOL)
})
