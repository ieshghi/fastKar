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

## sep_stats(): KS D and quantile margin of the simulated ecDNA vs HSR scores
test_that('sep_stats gives D = 1 and a positive margin for disjoint scores, D < 1 when they overlap', {
    s = sep_stats(ec = 10 + 1:50, hsr = -(1:50))
    expect_equal(s$sep_D, 1)
    expect_gt(s$sep_margin, 0)
    s = sep_stats(ec = 1:50, hsr = 26:75)
    expect_lt(s$sep_D, 1)
    expect_lt(s$sep_margin, 0)
})

## Without Hi-C (true_hic_path = NULL) a depth must be given
test_that('fastKar_ecfinder without Hi-C requires depth', {
    fx = load_fixture('toy_bfb')
    expect_error(fastKar_ecfinder(fx$gg, fx$ft, hic.res = -100, true_hic_path = NULL), 'depth is required')
})

## squeeze(keep.acen = FALSE): only karyotypes whose amplicon-heaviest walk carries a centromere (user, 2026-11-02)
test_that('amp_walk_centric looks at the walk with the most amplicon bp, copy-weighted', {
    w = c(10, 10, 100, 10) # node widths; amplicon = nodes 2, 3; centromere = node 4
    # walk 1: amplicon node 3 twice (200 bp) and no centromere; walk 2: node 2 (10 bp) and the centromere
    kar = list(snode.id = list(c(1, 3, -3), c(2, 4)), circular = c(FALSE, FALSE))
    expect_false(amp_walk_centric(kar, amp_nodes = 2:3, cen_nodes = 4, node_width = w))
    kar = list(snode.id = list(c(1, 3, -3, 4), c(2)), circular = c(FALSE, FALSE))
    expect_true(amp_walk_centric(kar, amp_nodes = 2:3, cen_nodes = 4, node_width = w))
    expect_false(amp_walk_centric(list(snode.id = list(c(1, 4)), circular = FALSE), 2:3, 4, w)) # no amplicon
})

test_that('fastKar_centromeres infers the build from chr1 and follows the graph chr naming', {
    mk = function(chr1len, prefix){
        nodes = GRanges(paste0(prefix, c('1', '8')), IRanges(1, c(chr1len, 1e8)))
        GenomeInfoDb::seqlengths(nodes) = c(chr1len, 145138636)
        gG(nodes = nodes)
    }
    c38 = fastKar_centromeres(mk(248956422, ''))
    expect_true('8' %in% as.character(seqnames(c38)))
    expect_equal(start(c38[seqnames(c38) == '8']), 43200001) # UCSC hg38 acen 8p11.1-q11.1
    c19 = fastKar_centromeres(mk(249250000, 'chr')) # JaBbA kb-rounded hg19 length
    expect_true('chr8' %in% as.character(seqnames(c19)))
    expect_equal(start(c19[seqnames(c19) == 'chr8']), 43100001) # UCSC hg19
    expect_error(fastKar_centromeres(mk(2e8, '')), 'cannot tell')
    expect_identical(c38, fastKar_centromeres(mk(2e8, ''), build = 'hg38'))
})

test_that('squeeze with keep.acen = FALSE returns only karyotypes with a centric amplicon walk', {
    fx = load_fixture('ccle_dm')
    amp = unique((fx$gg$nodes$gr %&% fx$ft)$node.id)
    # a fake centromere on one non-amplicon node that some but not all karyotypes route the amplicon through is hard to
    # guarantee, so check the invariant: every returned hypothesis passes amp_walk_centric for the centromeres used
    cen = fx$gg$nodes$gr[setdiff(fx$gg$nodes$dt[cn > 0]$node.id, amp)][1:20]
    set.seed(1)
    hs = suppressWarnings(squeeze(fx$gg, fx$ft, N = 200, k_return = 5, keep.acen = FALSE, centromeres = cen))
    cen_nodes = unique((fx$gg$nodes$gr %&% cen)$node.id)
    w = as.numeric(width(fx$gg$nodes$gr))
    for (h in hs) expect_true(amp_walk_centric(list(snode.id = h$snode.id, circular = h$circular), amp, cen_nodes, w))
    set.seed(1)
    expect_warning(none <- squeeze(fx$gg, fx$ft, N = 50, k_return = 5, keep.acen = FALSE, centromeres = GRanges('22', IRanges(1, 2))), 'centric')
    expect_length(none, 0)
    set.seed(1)
    expect_length(squeeze(fx$gg, fx$ft, N = 50, k_return = 5, keep.acen = TRUE), 5)
})
