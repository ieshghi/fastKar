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

## ---- paired ecDNA / HSR hypotheses (boil_paired, outside_signature, group leave-out) ----

## chr1 = nodes 1-10 (100 kb tiles), chr2 = nodes 11-12; amplicon = nodes 4, 5 (CN 3, ALT edge 5 -> 4 CN 2) on the one
## linear walk W, the chr2 walk is untouched by the amplicon
toy_paired = function() {
    tiles = gr.tile(parse.gr('1:1-1e6,2:1-2e5'), 1e5)
    W = c(1, 2, 3, 4, 5, 4, 5, 4, 5, 6, 7, 8, 9, 10)
    gg = makewalk(tiles, list(W, 11:12))$graph
    ft = GRanges('1', IRanges(3e5 + 1, 5e5))
    hsr = gW(graph = gg, snode.id = list(W, 11:12), circular = c(FALSE, FALSE))
    list(gg = gg, ft = ft, hsr = hsr, amp = unique((gg$nodes$gr %&% ft)$node.id), W = W)
}

## chr1 = nodes 1-5, chr2 = nodes 6-10; the amplicon (nodes 3 and 8) sits on two linear walks, one per chromosome, joined by
## ALT edges 3 -> 8 and 8 -> 3: the circle (3, 8) uses both walks
toy_two_walks = function() {
    tiles = gr.tile(parse.gr('1:1-5e5,2:1-5e5'), 1e5)
    X = c(1, 2, 3, 8, 9, 10)
    Y = c(6, 7, 8, 3, 4, 5)
    gg = makewalk(tiles, list(X, Y))$graph
    ft = c(GRanges('1', IRanges(2e5 + 1, 3e5)), GRanges('2', IRanges(2e5 + 1, 3e5)))
    list(gg = gg, ft = ft, hsr = gW(graph = gg, snode.id = list(X, Y), circular = c(FALSE, FALSE)),
         amp = unique((gg$nodes$gr %&% ft)$node.id))
}

test_that('touch_subgraph CN = traversal counts of the amplicon-touching walks, nodes and edges', {
    tp = toy_paired()
    before = list(n = tp$gg$nodes$dt$cn, e = tp$gg$edges$dt$cn)
    sub = touch_subgraph(tp$hsr, tp$amp)
    expect_equal(sub$nodes$dt$cn, tabulate(abs(tp$W), 12)) # chr2 walk is not touched: CN 0
    e = table(abs(unlist(tp$hsr$sedge.id[1])))
    expect_equal(sub$edges$dt$cn, as.integer(ifelse(is.na(e[as.character(sub$edges$dt$edge.id)]), 0, e[as.character(sub$edges$dt$edge.id)])))
    expect_equal(sub$nodes$dt$cn[4:5], c(3, 3))
    expect_equal(sub$edges$dt[n1 == 5 & n2 == 4]$cn, 2)
    expect_equal(sub$edges$dt[n1 == 11]$cn, 0)
    expect_identical(list(n = tp$gg$nodes$dt$cn, e = tp$gg$edges$dt$cn), before) # the input graph is untouched
    # a repeated visit counts once per traversal, also with two touching walks
    tw = toy_two_walks()
    sub = touch_subgraph(tw$hsr, tw$amp)
    expect_equal(sub$nodes$dt$cn, tabulate(abs(unlist(tw$hsr$snode.id)), 10))
    expect_equal(sub$edges$dt[(n1 == 3 & n2 == 8) | (n1 == 8 & n2 == 3)]$cn, c(1, 1))
})

test_that('boil_paired: every ecDNA hypothesis keeps the other walks of its HSR partner unchanged', {
    tp = toy_paired()
    set.seed(1)
    circles = find_circles(tp$gg, tp$ft, N = 50)
    hyps = boil_paired(tp$hsr, circles, tp$ft, n_boil = 3)
    expect_gt(length(hyps), 0)
    for (h in hyps) {
        expect_true(any(vapply(h$snode.id, function(s) isTRUE(all.equal(as.numeric(s), c(11, 12))), logical(1))))
        expect_true(any(h$circular))
        expect_identical(h$graph, tp$hsr$graph) # walks are on the parent graph, not the sub-graph
    }
    expect_false(anyDuplicated(vapply(hyps, function(h) h$hash, character(1))) > 0)
    # and on a real graph: walks not touching the amplicon come first and unchanged
    fx = load_fixture('ccle_dm')
    amp = unique((fx$gg$nodes$gr %&% fx$ft)$node.id)
    set.seed(1)
    hsr = suppressWarnings(squeeze(fx$gg, fx$ft, N = 100, k_return = 1, keep.acen = TRUE))[[1]]
    circ = suppressWarnings(find_circles(fx$gg, fx$ft, N = 200))
    others = hsr$snode.id[!vapply(hsr$snode.id, function(s) any(abs(s) %in% amp), logical(1))]
    hyps = suppressWarnings(boil_paired(hsr, circ, fx$ft, n_boil = 2))
    expect_gt(length(hyps), 0)
    for (h in hyps) expect_identical(unname(h$snode.id[seq_along(others)]), unname(others))
})

test_that('the first boil is the deterministic peel, with parent-graph node ids', {
    tp = toy_paired()
    set.seed(1)
    circles = find_circles(tp$gg, tp$ft, N = 50)
    sub = touch_subgraph(tp$hsr, tp$amp)
    p = greedypeel(sub, circles, tp$ft, deterministic = TRUE)
    expect_identical(p$rest$nodes$gr[, c()], tp$gg$nodes$gr[, c()]) # loosefix does not renumber nodes
    expect_identical(p$rest$nodes$dt$node.id, tp$gg$nodes$dt$node.id)
    hyps = boil_paired(tp$hsr, circles, tp$ft, n_boil = 3)
    k = length(p$circles)
    h1 = hyps[[1]]
    circ_walks = h1$snode.id[h1$circular]
    expect_equal(unname(tail(h1$snode.id, k)), unname(p$circles))
    expect_equal(length(circ_walks) >= k, TRUE)
})

test_that('a toy amplicon split over two linear walks is boiled from both', {
    tw = toy_two_walks()
    set.seed(1)
    circles = find_circles(tw$gg, tw$ft, N = 50)
    sub = touch_subgraph(tw$hsr, tw$amp)
    p = greedypeel(sub, circles, tw$ft, deterministic = TRUE)
    expect_length(p$circles, 1)
    expect_setequal(abs(p$circles[[1]]), c(3, 8)) # the circle uses nodes of both walks
    hyps = boil_paired(tw$hsr, circles, tw$ft, n_boil = 2)
    expect_gt(length(hyps), 0)
    h = hyps[[1]]
    expect_true(any(vapply(h$snode.id[h$circular], function(s) setequal(abs(s), c(3, 8)), logical(1))))
    expect_equal(amp_circular_fraction(h, tw$ft), 0.5) # one of the two copies of the amplicon
})

test_that('boil_paired leaves an rds-loaded graph unchanged', {
    fx = load_fixture('ccle_dm')
    gg = fx$gg
    state = function(g) list(node.cn = g$nodes$dt$cn, edge.cn = g$edges$dt$cn)
    before = state(gg)
    set.seed(1)
    hsr = suppressWarnings(squeeze(gg, fx$ft, N = 100, k_return = 1, keep.acen = TRUE))[[1]]
    circ = suppressWarnings(find_circles(gg, fx$ft, N = 200))
    hyps = suppressWarnings(boil_paired(hsr, circ, fx$ft, n_boil = 2))
    expect_gt(length(hyps), 0)
    expect_identical(state(gg), before)
    expect_identical(state(hsr$graph), before)
})

test_that('outside_signature: equal when only the amplicon differs, orientation-normalised, different when rewired outside', {
    tp = toy_paired()
    set.seed(1)
    circles = find_circles(tp$gg, tp$ft, N = 50)
    hyp = boil_paired(tp$hsr, circles, tp$ft, n_boil = 1)[[1]]
    s_hsr = outside_signature(tp$hsr, tp$amp)
    s_hyp = outside_signature(hyp, tp$amp)
    expect_identical(as.character(s_hsr), as.character(s_hyp))
    expect_identical(as.numeric(attr(s_hsr, 'bp')), as.numeric(attr(s_hyp, 'bp')))
    expect_true(outside_mismatch(s_hsr, s_hyp)$match)
    expect_equal(outside_mismatch(s_hsr, s_hyp)$frac_bp, 0)
    # walk ends and amplicon boundaries are part of the key
    expect_true(all(grepl('^[EA]\\|.*\\|[EA]$', as.character(s_hsr))))
    # reversed walks give the same signature
    rev_hyp = gW(graph = tp$gg, snode.id = lapply(tp$hsr$snode.id, function(s) rev(-s)), circular = c(FALSE, FALSE))
    expect_identical(as.character(outside_signature(rev_hyp, tp$amp)), as.character(s_hsr))
    # rewired outside: node 2 (CN 2) is paired 1 -> 2 -> 3 and 5 -> 2 -> 6 in A, 1 -> 2 -> 6 and 5 -> 2 -> 3 in B; amplicon = node 7
    tiles = gr.tile(parse.gr('1:1-7e5'), 1e5)
    gg = makewalk(tiles, list(c(1, 2, 3, 4), c(5, 2, 6), 7))$graph
    A = gW(graph = gg, snode.id = list(c(1, 2, 3, 4), c(5, 2, 6), 7), circular = c(FALSE, FALSE, FALSE))
    B = gW(graph = gg, snode.id = list(c(1, 2, 6), c(5, 2, 3, 4), 7), circular = c(FALSE, FALSE, FALSE))
    sA = outside_signature(A, 7); sB = outside_signature(B, 7)
    expect_false(identical(as.character(sA), as.character(sB)))
    m = outside_mismatch(sA, sB)
    expect_false(m$match)
    expect_gt(m$frac_bp, 0)
    expect_lte(m$frac_bp, 1)
    expect_true(outside_mismatch(sA, outside_signature(A, 7))$match)
})

test_that('outside_signature handles circular walks: rotation, orientation, and amplicon-free circles', {
    tiles = gr.tile(parse.gr('1:1-4e5'), 1e5)
    gw = gW(grl = GRangesList(tiles[1:4]), circular = TRUE)
    gg = gw$graph
    mk = function(s) gW(graph = gg, snode.id = list(s), circular = TRUE)
    s0 = outside_signature(mk(c(1, 2, 3, 4)), 3)
    expect_identical(as.character(outside_signature(mk(c(2, 3, 4, 1)), 3)), as.character(s0)) # rotation
    expect_identical(as.character(outside_signature(mk(c(3, 4, 1, 2)), 3)), as.character(s0))
    expect_identical(as.character(outside_signature(mk(c(-4, -3, -2, -1)), 3)), as.character(s0)) # reverse complement
    expect_identical(as.character(s0), 'A|-2,-1,-4|A') # one run (4, 1, 2) although it wraps around the walk start; reverse complement is the smaller key
    c0 = outside_signature(mk(c(1, 2, 3, 4)), 99) # no amplicon node: one canonical key
    expect_length(c0, 1)
    expect_true(startsWith(as.character(c0), 'C|'))
    expect_identical(as.character(outside_signature(mk(c(3, 4, 1, 2)), 99)), as.character(c0))
    expect_identical(as.character(outside_signature(mk(c(-1, -4, -3, -2)), 99)), as.character(c0))
    expect_equal(as.numeric(attr(c0, 'bp')), 4e5)
})

test_that('amp_telo_frac and amp_centric_frac use linear walks, telomere ends and centromere nodes', {
    tp = toy_paired()
    expect_equal(amp_telo_frac(tp$hsr, tp$ft), 1) # walk 1..10 runs telomere to telomere
    expect_equal(amp_centric_frac(tp$hsr, tp$ft, cen_nodes = 1), 1)
    expect_equal(amp_centric_frac(tp$hsr, tp$ft, cen_nodes = 11), 0) # only the chr2 walk has node 11, no amplicon on it
    cut = gW(graph = tp$gg, snode.id = list(c(2, 3, 4, 5, 4, 5, 4, 5, 6, 7, 8, 9, 10), 11:12), circular = c(FALSE, FALSE))
    expect_equal(amp_telo_frac(cut, tp$ft), 0) # starts at node 2, a loose end
    rev = gW(graph = tp$gg, snode.id = list(rev(-tp$W), 11:12), circular = c(FALSE, FALSE))
    expect_equal(amp_telo_frac(rev, tp$ft), 1)
})

test_that('leave-out helpers: group, solution and none modes', {
    ec_group = c(1, 1, 2, 3, 3)
    expect_identical(pick_leave_out_mode(ec_group, 5, 3), 'group')
    expect_identical(pick_leave_out_mode(c(1, 1, 1), 3, 2), 'solution')
    expect_identical(pick_leave_out_mode(1, 1, 3), 'none')
    expect_identical(pick_leave_out_mode(c(1, 1), 2, 1), 'none')
    expect_identical(leave_out_ix('hsr', 1, ec_group, 'group'), list(hsr = 1, ec = 1:2))
    expect_identical(leave_out_ix('ec', 4, ec_group, 'group'), list(hsr = 3, ec = 4:5))
    expect_identical(leave_out_ix('ec', 3, ec_group, 'group'), list(hsr = 2, ec = 3L))
    expect_identical(leave_out_ix('hsr', 2, c(1, 1, 3), 'solution'), list(hsr = 2, ec = integer(0)))
    expect_identical(leave_out_ix('ec', 2, c(1, 1, 3), 'solution'), list(hsr = integer(0), ec = 2))
    expect_identical(leave_out_ix('ec', 2, 1, 'none'), list(hsr = integer(0), ec = integer(0)))
})

test_that('outside_mismatch counts the bp of runs present in different numbers', {
    a = structure(c('E|1|A', 'E|1|A', 'E|2|E'), bp = c(10, 10, 5))
    b = structure(c('E|1|A', 'E|3|E'), bp = c(10, 7))
    m = outside_mismatch(a, b)
    expect_false(m$match)
    expect_equal(m$frac_bp, (10 + 5 + 7) / (25 + 17))
    expect_true(outside_mismatch(a, a)$match)
    expect_equal(outside_mismatch(structure(character(0), bp = numeric(0)), structure(character(0), bp = numeric(0)))$frac_bp, 0)
})

## fastKar_ecfinder without Hi-C on the toy: scalar outputs, group leave-out, HSR hypotheses without a partner kept,
## CI early return when no ecDNA hypothesis passes
test_that('fastKar_ecfinder (paired): scalar outputs, partnerless HSR hypotheses kept, CI when no ecDNA hypothesis passes', {
    tp = toy_paired()
    run = function(...) suppressWarnings(suppressMessages(fastKar_ecfinder(tp$gg, tp$ft, hic.res = -10, true_hic_path = NULL, depth = 1,
        n_ec = 3, n_sample = 6, n_circle_samples = 50, n_hsr_samples = 100, keep.acen = TRUE, return_all = TRUE, seed = 1, hypotheses = "paired", ...)))
    scalars = c('n_ec', 'n_hsr', 'n_groups', 'leave_out', 'leave_out_mode', 'max_ecdna_circ', 'pair_match_frac', 'hsr_top_cov',
                'hsr_centric_frac', 'ec_centric_frac', 'hsr_telo_frac')
    o = run()
    for (nm in scalars) expect_length(o[[nm]], 1)
    expect_identical(o$leave_out_mode, 'group')
    expect_true(o$leave_out)
    expect_equal(o$n_groups, 3)
    expect_equal(o$pair_match_frac, 1)
    expect_equal(o$hyp$ec_group, sort(o$hyp$ec_group))
    expect_length(o$hyp$pair_outside_match, o$n_ec)
    expect_length(o$hyp$hsr_cov, o$n_hsr)
    expect_equal(o$hsr_telo_frac, 1)
    expect_s3_class(data.table::as.data.table(o[setdiff(names(o), c('hyp'))]), 'data.table')
    # squeeze() can return the same karyotype several times; use three distinct HSR hypotheses (the walk W cut at different
    # places) and let HSR hypothesis 2 get no ecDNA partner: it is kept, only groups 1 and 3 remain
    cuts = list(list(tp$W, 11:12), list(tp$W[1:9], tp$W[10:14], 11:12), list(tp$W[1:11], tp$W[12:14], 11:12))
    three = lapply(cuts, function(k) gW(graph = tp$gg, snode.id = k, circular = rep(FALSE, length(k))))
    skip_hash = three[[2]]$hash
    orig = boil_paired
    testthat::local_mocked_bindings(
        squeeze = function(...) three,
        boil_paired = function(gw, circles, ft, n_boil = 3) {
            if (gw$hash == skip_hash) list() else orig(gw, circles, ft, n_boil)
        })
    o2 = run()
    expect_equal(o2$n_hsr, 3)
    expect_equal(o2$n_groups, 2)
    expect_false(2 %in% o2$hyp$ec_group)
    # nothing passes the gate: CI with the scalar diagnostics and the hypotheses
    o3 = run(ec_fraction_threshold = 1.01)
    expect_identical(o3$call, 'CI')
    expect_identical(o3$reason, 'No good ecDNA solutions')
    expect_equal(o3$n_groups, 0)
    expect_equal(o3$n_hsr, 3)
    expect_true(is.na(o3$ec_centric_frac))
    expect_true(all(c('hsr_solns', 'ec_group', 'hsr_cov', 'ecdna_candidates', 'ecdna_circ') %in% names(o3$hyp)))
    for (nm in c('max_ecdna_circ', 'n_groups', 'pair_match_frac', 'hsr_top_cov', 'hsr_telo_frac')) expect_length(o3[[nm]], 1)
})

## call rule (user, 2026-11-02): quantiles of the simulated scores, so one outlying sample does not decide the call
test_that('ecfinder_call uses the call_q quantiles, not the extremes', {
    hsr = c(-200, 1:99 - 150) # one far-left HSR sample
    ec = c(-1000, 101:199)    # one far-left ecDNA sample, the rest above every HSR sample
    expect_identical(ecfinder_call(-120, ec, hsr, TRUE), 'CI')      # inside HSR; min/max rule: ecDNA (> -1000)
    expect_identical(ecfinder_call(150, ec, hsr, TRUE), 'ecDNA')
    expect_identical(ecfinder_call(0, ec = 1:100 - 20, hsr = 1:100 - 80, separated = TRUE), 'None') # between ecDNA q05 and HSR q95
    expect_identical(ecfinder_call(150, ec, hsr, FALSE), 'None')    # not separated
})

test_that('fastKar_ecfinder hypotheses = independent runs and reports no pairing', {
    tp = toy_paired()
    o = suppressWarnings(suppressMessages(fastKar_ecfinder(tp$gg, tp$ft, hic.res = -10, true_hic_path = NULL, depth = 1,
        n_ec = 3, n_sample = 6, n_circle_samples = 50, n_hsr_samples = 100, ec_fraction_threshold = 0, seed = 1)))
    expect_identical(o$hypotheses, 'independent')
    expect_true(is.na(o$n_groups))
    expect_true(is.na(o$pair_match_frac))
    expect_identical(o$leave_out_mode, if (o$n_ec > 1 & o$n_hsr > 1) 'solution' else 'none')
})
