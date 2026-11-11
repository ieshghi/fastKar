## Amplicon callers (R/integration.R): fastKar_ecdna() and fastKar_integration_tiles(), on a toy 1 Mb mcool written to a tempdir

fastkar_toy = function(){
    # 10 chromosomes x 20 Mb, CN 2 (chr5: CN 3 after 6 Mb, a breakpoint); amplicon chr1:10.0-10.3 Mb at CN 20 (< 1 Mb: its one tile is the amplicon tile)
    chs = as.character(1:10); tl = data.table(chrom = rep(chs, each = 20), start = rep(0:19, 10)*1e6)[, end := start + 1e6]
    nodes = GRanges(c(chs[-c(1, 5)], '1', '1', '1', '5', '5'), IRanges(c(rep(1, 8), 1, 10e6+1, 10.3e6+1, 1, 6e6+1), c(rep(20e6, 8), 10e6, 10.3e6, 20e6, 6e6, 20e6)))
    gg = list(nodes = list(gr = nodes, dt = data.table(seqnames = as.character(seqnames(nodes)), start = start(nodes), end = end(nodes), cn = c(rep(2, 8), 2, 20, 2, 2, 3))))
    cnbp = rep(2e6, nrow(tl)); a = which(tl$chrom == '1' & tl$start == 10e6); cnbp[a] = 0.3e6*20 + 0.7e6*2; cnbp[tl$chrom == '5' & tl$start >= 6e6] = 3e6
    mk = function(rate_fun){  # reads per tile pair (all pairs i < j) = 1e-11 x CN.bp x CN.bp x rate
        f = tempfile(fileext = '.mcool'); G = 'resolutions/1000000'
        rhdf5::h5createFile(f); for (g in c('resolutions', G, paste0(G, '/bins'), paste0(G, '/pixels'))) rhdf5::h5createGroup(f, g)
        rhdf5::h5write(paste0('chr', tl$chrom), f, paste0(G, '/bins/chrom')); rhdf5::h5write(tl$start, f, paste0(G, '/bins/start')); rhdf5::h5write(tl$end, f, paste0(G, '/bins/end'))
        p = CJ(i = seq_len(nrow(tl)), j = seq_len(nrow(tl)))[i < j]
        p[, v := round(1e-11*cnbp[i]*cnbp[j]*rate_fun(i, j))]
        rhdf5::h5write(p$i - 1L, f, paste0(G, '/pixels/bin1_id')); rhdf5::h5write(p$j - 1L, f, paste0(G, '/pixels/bin2_id')); rhdf5::h5write(p$v, f, paste0(G, '/pixels/count'))
        rhdf5::h5closeAll(); f}
    list(tl = tl, gg = gg, a = a, mk = mk, amp = GRanges('1', IRanges(10e6+1, 10.3e6)), site = which(tl$chrom == '5' & tl$start %in% c(5e6, 6e6)))
}

test_that('fastKar_ecdna calls a genome-wide amplicon excess ecDNA and a single-site excess not', {
    skip_if_not_installed('rhdf5'); y = fastkar_toy(); a = y$a; site = y$site
    # ecDNA: every genome tile contacts the amplicon 2x more per copy-bp than it contacts the genome -> r = 1 on all 9 trans chromosomes
    z = fastKar_ecdna(y$gg, y$mk(function(i, j) ifelse(i == a | j == a, 2, 1)), y$amp)
    expect_equal(z$n_amp_tiles, 1); expect_equal(z$n_chr, 9); expect_equal(z$n_chr_pos, 9)
    expect_equal(unique(round(z$tiles$r, 6)), 1); expect_equal(z$call, 'ecDNA')
    expect_false(any(z$tiles$seqnames == '1'))  # chr1 genome tiles are cis to the only amplicon tile
    # integration: amplicon contacts 0.9x everywhere except a 2 Mb site on chr5 at 20x -> pooled score > 0, but every chromosome median < 0
    z = fastKar_ecdna(y$gg, y$mk(function(i, j) ifelse(i == a | j == a, ifelse(i %in% site | j %in% site, 20, 0.9), 1)), y$amp)
    expect_gt(z$score, 0); expect_equal(z$n_chr_pos, 0); expect_equal(z$call, 'not ecDNA')
    # mask: a tile with >= 500 kb masked is dropped
    z = fastKar_ecdna(y$gg, y$mk(function(i, j) 1), y$amp, mask = GRanges('3', IRanges(4e6+1, 4.6e6)))
    expect_false(any(z$tiles$seqnames == '3' & z$tiles$start == 4e6))
})

test_that('fastKar_integration_tiles calls a breakpoint site and absorbs a uniform excess', {
    skip_if_not_installed('rhdf5'); y = fastkar_toy(); a = y$a; site = y$site
    # uniform 2x excess with the amplicon: absorbed by the NB scale r, no site (depth = 0: no cis term in the toy)
    z = fastKar_integration_tiles(y$gg, y$mk(function(i, j) ifelse(i == a | j == a, 2, 1)), y$amp, depth = 0)
    expect_equal(z$call, 'no site'); expect_equal(z$amp_cn, 20)
    # a 2 Mb site on chr5 at the 6 Mb breakpoint (20x vs 0.9x elsewhere): both tiles called, one locus
    z = fastKar_integration_tiles(y$gg, y$mk(function(i, j) ifelse(i == a | j == a, ifelse(i %in% site | j %in% site, 20, 0.9), 1)), y$amp, depth = 0)
    expect_equal(z$call, 'Integrated'); expect_equal(sort(z$sites$start), c(5e6, 6e6)); expect_equal(nrow(z$loci), 1)
    # the same excess away from any breakpoint (chr7, no CN change) is not tested
    s7 = which(y$tl$chrom == '7' & y$tl$start %in% c(5e6, 6e6))
    z = fastKar_integration_tiles(y$gg, y$mk(function(i, j) ifelse(i == a | j == a, ifelse(i %in% s7 | j %in% s7, 20, 0.9), 1)), y$amp, depth = 0)
    expect_equal(z$call, 'no site'); expect_false(any(z$tiles[seqnames == '7' & start %in% c(5e6, 6e6)]$tested))
})

## ---- boil_greedy / squeeze and helpers (restored 2026-11-11 from the archived test-ecfinder.R) ----

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
