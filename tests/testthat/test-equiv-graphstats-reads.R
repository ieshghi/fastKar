## Differential tests (vs the oracle) for reads_fromwalk: words, multiplicities and row order
## must be identical, for linear and circular walks, minsize, use.nodes, plain-list walks and
## read lengths shorter and longer than the walks. Slow cases need FASTKAR_SLOW_TESTS=1.

SLOW = Sys.getenv('FASTKAR_SLOW_TESTS') == '1'

expect_reads_identical = function(new, old, info = NULL) {
    expect_identical(new, old, info = info)
    ## and after canonical sorting (what callers that key by word rely on)
    a = data.table::copy(new); b = data.table::copy(old)
    data.table::setorder(a, words); data.table::setorder(b, words)
    expect_identical(a, b, info = info)
}

reads_cases = function(fx) {
    nid = fx$gg$nodes$dt$node.id
    wd = fx$gg$nodes$dt$width
    list(list(minsize = 0, use.nodes = NULL),
         list(minsize = stats::median(wd), use.nodes = NULL),
         list(minsize = 0, use.nodes = nid[seq(1, length(nid), by = 2)]),
         list(minsize = stats::quantile(wd, 0.25, names = FALSE), use.nodes = nid[seq_len(ceiling(length(nid) / 2))]))
}

test_that('reads_fromwalk is identical to the oracle on toys (1e3, 1e5, 1e7; minsize; use.nodes)', {
    skip_if_no_oracle()
    for (nm in fixture_names(toy_only = TRUE)) {
        fx = load_fixture(nm)
        for (k in seq_along(fx$gw)) {
            gw = fx$gw[[k]]
            for (L in c(1, 1e3, 1e5, 1e7)) for (cs in reads_cases(fx)) {
                info = paste(nm, k, L, cs$minsize, length(cs$use.nodes))
                old = ref('reads_fromwalk')(gw, L, minsize = cs$minsize, use.nodes = cs$use.nodes)
                new = opt('reads_fromwalk')(gw, L, minsize = cs$minsize, use.nodes = cs$use.nodes)
                expect_reads_identical(new, old, info = info)
            }
        }
    }
})

test_that('reads_fromwalk is identical on circular walks around the walk length', {
    skip_if_no_oracle()
    fx = load_fixture('toy_circ')
    gw = fx$gw[[1]]
    tot = sum(fx$gg$nodes$dt$width[abs(gw$snode.id[[1]])])  # circular walk length
    for (L in c(tot - 1, tot, tot + 1, 2 * tot, 199999, 200000, 200001)) {
        expect_reads_identical(opt('reads_fromwalk')(gw, L), ref('reads_fromwalk')(gw, L), info = L)
    }
    ## a BFB with a circular walk rebuilt from its nodes, and a single-node cycle
    fx = load_fixture('toy_bfb')
    w = list(graph = fx$gg, snode.id = list(c(1L, 2L, -2L, 3L), 4L, c(-5L, -4L)), circular = c(TRUE, TRUE, FALSE))
    for (L in c(1e3, 1e5, 3e5, 1e7)) {
        expect_reads_identical(opt('reads_fromwalk')(w, L), ref('reads_fromwalk')(w, L), info = L)
        expect_reads_identical(opt('reads_fromwalk')(w, L, use.nodes = 2:4), ref('reads_fromwalk')(w, L, use.nodes = 2:4), info = L)
    }
})

test_that('reads_fromwalk is identical on CCLE fixtures (plain-list walks as in identifiability.R)', {
    skip_if_no_oracle()
    nms = grep('^ccle_', fixture_names(), value = TRUE)
    if (!SLOW) nms = c('ccle_pyrgo_local', 'ccle_dm_local', 'ccle_chromothripsis_local', 'ccle_tyfonas_local', 'ccle_pyrgo')
    for (nm in nms) {
        fx = load_fixture(nm)
        gw = fx$gw[[1]]
        use = fx$gg$nodes$dt[GenomicRanges::start(fx$gg$nodes$gr) <= max(GenomicRanges::end(fx$ft))]$node.id
        Ls = if (SLOW) c(1e3, 1e5, 1e7) else c(1e3, 1e5)
        for (L in Ls) {
            expect_reads_identical(opt('reads_fromwalk')(gw, L), ref('reads_fromwalk')(gw, L), info = paste(nm, L))
            w = list(graph = gw$graph, snode.id = gw$snode.id, circular = gw$circular)
            expect_reads_identical(opt('reads_fromwalk')(w, L, use.nodes = use, minsize = 1e4),
                                   ref('reads_fromwalk')(w, L, use.nodes = use, minsize = 1e4), info = paste(nm, L, 'use'))
        }
    }
})

test_that('reads_fromwalk at 1e7 on a full CCLE graph', {
    skip_if_no_oracle()
    fx = load_fixture('ccle_chromothripsis')
    gw = fx$gw[[2]]
    expect_reads_identical(opt('reads_fromwalk')(gw, 1e7), ref('reads_fromwalk')(gw, 1e7))
})

test_that('reads_fromwalk falls back to the reference for non-integer or tiny read lengths', {
    skip_if_no_oracle()
    fx = load_fixture('toy_repdup')
    gw = fx$gw[[1]]
    for (L in c(12345.5, 0.5)) {
        old = tryCatch(ref('reads_fromwalk')(gw, L), error = function(e) e)
        new = tryCatch(opt('reads_fromwalk')(gw, L), error = function(e) e)
        if (inherits(old, 'error')) expect_s3_class(new, 'error') else expect_reads_identical(new, old, info = L)
    }
    expect_reads_identical(opt('reads_fromwalk')(gw, 12345L), ref('reads_fromwalk')(gw, 12345L))
})

test_that('longread_probdist and longread_kl agree with the oracle', {
    skip_if_no_oracle()
    for (nm in c('toy_bfb', 'toy_threeway', 'ccle_dm_local')) {
        fx = load_fixture(nm)
        a = fx$gw[[1]]; b = fx$gw[[2]]
        for (L in c(1e4, 1e5)) {
            expect_identical(opt('longread_probdist')(list(a, b), L), ref('longread_probdist')(list(a, b), L), info = nm)
            expect_identical(opt('longread_kl')(a, b, readL = L, depth = 10, use.nodes = NULL),
                             ref('longread_kl')(a, b, readL = L, depth = 10, use.nodes = NULL), info = nm)
        }
    }
})
