## prep_data_for_py(): the problem folder for the Python search (pure R).

test_that('prep_data_for_py writes every table and the truth files', {
    skip_if_no_loosefix()
    data.table::setDTthreads(1)
    x = kinf_toy()
    ex = kinf_export(x)
    fs = c('tiled_target', 'dedup_tiles', 'nodes', 'internal_edges', 'external_edges', 'bins', 'gm_dat', 'depth', 'purity')
    expect_true(all(file.exists(file.path(ex$dir, paste0(fs, '.csv')))))
    expect_true(file.exists(file.path(ex$dir, 'truth_gwalks.json')))
    expect_true(file.exists(file.path(ex$dir, 'truth_expected.csv')))
    ## every half-node copy appears once; internal edges pair left and right copies of each node
    nd = ex$out$nodes
    expect_equal(nrow(nd), 2 * sum(x$gg$nodes$dt$cn))
    expect_equal(sort(nd$unique.id), seq_len(nrow(nd)))
    ie = ex$out$internal_edges
    expect_equal(nrow(ie), sum(x$gg$nodes$dt$cn))
    expect_true(all(ie$copies == x$gg$nodes$dt$cn[match(ie$node_id, x$gg$nodes$dt$node.id)]))
})

test_that('the observed counts are exported on the data bins as they are, and every tile lies in one bin', {
    skip_if_no_loosefix()
    data.table::setDTthreads(1)
    ## data on 300 kb bins, simulation tiles of 100 kb: the bins are scored, the tiles nest in them
    x = kinf_toy(pix = 1e5, sim_pix = 3e5)
    ex = kinf_export(x)
    expect_equal(nrow(ex$out$bins), length(x$gm$gr))
    expect_gt(nrow(ex$out$tiled_target), nrow(ex$out$bins))
    expect_equal(sum(ex$out$gm_dat$value), sum(x$gm$dat$value), tolerance = 1e-10)
    obs = merge(x$gm$dat[, .(i, j, value)], ex$out$gm_dat[, .(i, j, v2 = value)], by = c('i', 'j'), all = TRUE)
    expect_equal(obs[!is.na(value)]$value, obs[!is.na(value)]$v2) # no split, no rescale
    tt = ex$out$tiled_target; b = ex$out$bins
    expect_true(all(tt$start >= b$start[tt$bin.id] & tt$end <= b$end[tt$bin.id] & tt$seqnames == b$seqnames[tt$bin.id]))
    expect_equal(sum(tt$width), sum(b$width))
})

test_that('truth_gwalks.json round-trips through read_kinf_gwalks', {
    skip_if_no_loosefix()
    data.table::setDTthreads(1)
    for (nm in c('toy_repdup', 'ccle_pyrgo_local')) {
        x = if (nm == 'toy_repdup') kinf_toy(nm, pix = 5e4) else kinf_toy(nm, k = 'circ', pix = NULL)
        ex = kinf_export(x)
        raw = read_kinf_gwalks(file.path(ex$dir, 'truth_gwalks.json'))
        expect_equal(raw$snode.id, lapply(x$truth$snode.id, as.integer))
        expect_equal(raw$circular, as.logical(x$truth$circular))
        gw = read_kinf_gwalks(file.path(ex$dir, 'truth_gwalks.json'), x$gg)
        expect_s3_class(gw, 'gWalk')
        expect_equal(gw$circular, x$truth$circular)
        expect_equal(lapply(gw$snode.id, as.integer), lapply(x$truth$snode.id, as.integer))
    }
})

test_that('infer_karyotype rejects purity < 1 and unbalanced graphs', {
    skip_if_no_loosefix()
    x = kinf_toy()
    expect_error(infer_karyotype(x$gg, x$gm, depth = 10, purity = 0.8, python = '/bin/true'), 'purity')
    bad = x$gg$copy
    bad$edges$mark(cn = bad$edges$dt$cn + 5L)
    expect_error(opt('.kinf_check_graph')(bad), 'not balanced')
    expect_no_error(opt('.kinf_check_graph')(x$gg))
})
