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
