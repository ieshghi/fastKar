## rebin_matrix() treats pixel values as counts: rebinning conserves the total, nested bins are
## summed exactly, and split pixels are shared in proportion to the overlap.

pix_value = function(gm, a, b) {
    d = gm$dat
    v = d[i == a & j == b]$value
    if (length(v)) v else 0
}

## contacts at positions (x <= y) counted into upper-triangle pixels of bins gr
bin_contacts = function(x, y, gr) {
    bx = gr.match(GRanges(1, IRanges(x, width = 1)), gr)
    by = gr.match(GRanges(1, IRanges(y, width = 1)), gr)
    d = data.table(i = pmin(bx, by), j = pmax(bx, by))[, .(value = as.double(.N)), by = .(i, j)]
    rebin_matrix(gM(gr = gr, dat = d), gr)
}

test_that('a half-bin shift splits pixels by overlap and conserves the total', {
    src = GRanges(1, IRanges(c(5001, 15001), width = 1e4))
    tgt = GRanges(1, IRanges(c(1, 10001, 20001), width = 1e4))
    d = data.table(i = c(1L, 1L, 2L), j = c(1L, 2L, 2L), value = c(4, 8, 4))
    m = gM(gr = src, dat = d)
    r = rebin_matrix(m, tgt)
    expect_equal(sum(r$dat$value), 16)
    ## each source bin is half in two target bins (f = 1/2): a diagonal pixel V gives V/4 to each
    ## within-bin target pixel and V/2 to the pixel between them; an off-diagonal pixel V gives V/4 to each of 4
    expect_equal(pix_value(r, 1, 1), 1)
    expect_equal(pix_value(r, 1, 2), 2 + 2)
    expect_equal(pix_value(r, 2, 2), 1 + 2 + 1)
    expect_equal(pix_value(r, 1, 3), 2)
    expect_equal(pix_value(r, 2, 3), 2 + 2)
    expect_equal(pix_value(r, 3, 3), 1)
})

test_that('coarsening onto nested bins equals binning the contacts directly', {
    set.seed(1)
    n = 5000
    x = sample.int(1e5, n, replace = TRUE); y = sample.int(1e5, n, replace = TRUE)
    lo = pmin(x, y); hi = pmax(x, y)
    fine = GRanges(1, IRanges(c(1, 3001, 10001, 12001, 25001, 40001, 41001, 70001), c(3000, 10000, 12000, 25000, 40000, 41000, 70000, 1e5)))
    coarse = GRanges(1, IRanges(c(1, 25001, 70001), c(25000, 70000, 1e5)))
    from_fine = rebin_matrix(bin_contacts(lo, hi, fine), coarse)
    direct = bin_contacts(lo, hi, coarse)
    expect_equal(from_fine$dat[order(i, j)]$value, direct$dat[order(i, j)]$value)
    expect_equal(sum(from_fine$dat$value), n)
})

test_that('unaligned and finer target bins conserve the total', {
    set.seed(2)
    x = sample.int(1e5, 3000, replace = TRUE); y = sample.int(1e5, 3000, replace = TRUE)
    src = bin_contacts(pmin(x, y), pmax(x, y), gr.tile(GRanges(1, IRanges(1, 1e5)), 7000))
    shifted = GRanges(1, IRanges(seq(-2999, 1e5, by = 1e4), width = 1e4))
    shifted = shifted[end(shifted) > 0]; start(shifted)[1] = 1
    expect_equal(sum(rebin_matrix(src, shifted)$dat$value), 3000)
    expect_equal(sum(rebin_matrix(src, gr.tile(GRanges(1, IRanges(1, 1e5)), 2500))$dat$value), 3000)
})

test_that('parts outside the target are dropped and same bins are returned unchanged', {
    src = GRanges(1, IRanges(c(1, 10001), width = 1e4))
    d = data.table(i = c(1L, 1L, 2L), j = c(1L, 2L, 2L), value = c(3, 5, 7))
    m = gM(gr = src, dat = d)
    r = rebin_matrix(m, src[1])
    expect_equal(r$dat$value, 3)
    same = rebin_matrix(m, src)
    expect_equal(same$dat[order(i, j)]$value, c(3, 5, 7))
})

test_that('a simulated map rebinned onto cooler-style bins keeps its total', {
    fx = load_fixture('ccle_dm_local')
    res = 1e4
    ctx = gUtils::streduce(fx$ft)
    bins = gr.tile(GRanges(seqnames(ctx), IRanges(floor((start(ctx) - 1) / res) * res + 1, ceiling(end(ctx) / res) * res)), res)
    s = suppressWarnings(forward_simulate(fx$gw[[1]], target_region = ctx, pix.size = res, depth = 10))
    inside = sum(suppressWarnings(rebin_matrix(s, bins))$dat$value)
    expect_equal(inside, sum(s$dat$value), tolerance = 1e-10)
    ## simulating on the bins themselves nests every tile in one bin: exact sum, same total
    s2 = suppressWarnings(forward_simulate(fx$gw[[1]], target_region = bins, pix.size = res, depth = 10))
    expect_equal(sum(suppressWarnings(rebin_matrix(s2, bins))$dat$value), sum(s2$dat$value), tolerance = 1e-10)
})
