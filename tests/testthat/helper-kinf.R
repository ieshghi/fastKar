## Shared helpers for the kinf tests: a toy fixture with a simulated noisy map, exported with a truth.

## k = 'circ': the first karyotype with a circular walk; pix = NULL: ~60 tiles over the footprint
kinf_toy = function(name = 'toy_repdup', k = 1, depth = 10, pix = 1e5, sim_pix = pix, seed = 1) {
    fx = load_fixture(name)
    if (identical(k, 'circ')) k = which(vapply(fx$gw, function(w) any(w$circular), logical(1)))[1]
    if (is.null(pix)) pix = sim_pix = max(5e3, round(sum(width(fx$ft)) / 60 / 1e3) * 1e3)
    truth = fx$gw[[k]]
    set.seed(seed)
    gm = make_noisymap(forward_simulate(truth, target_region = fx$ft, depth = depth, pix.size = sim_pix), 2)[[1]]
    list(gg = fx$gg, ft = fx$ft, truth = truth, gm = gm, depth = depth, pix = pix)
}

kinf_export = function(x, dir = tempfile('kinf_test_')) {
    dir.create(dir, recursive = TRUE)
    out = prep_data_for_py(gg = x$gg, gm = x$gm, depth = x$depth, pix.size = x$pix, target_region = x$ft,
                           write.to.folder = paste0(dir, '/'), truth = x$truth)
    list(dir = dir, out = out)
}

## R reference NLL of an exported folder: area-weighted NB NLL (compdats, theta 2, area0 = median tile width^2)
## of the observed map vs truth_expected.csv on the export tiling (as KI scripts/ref_nll.R)
kinf_ref_nll = function(dir) {
    tt = data.table::fread(file.path(dir, 'tiled_target.csv'))
    gr = GenomicRanges::GRanges(tt$seqnames, IRanges::IRanges(tt$start, tt$end))
    tmpl = data.table::as.data.table(make_template_dat_cpp(data.table::data.table(width = as.numeric(width(gr)))))
    obs = data.table::fread(file.path(dir, 'gm_dat.csv'))
    ex = data.table::fread(file.path(dir, 'truth_expected.csv'))
    d1 = data.table::copy(tmpl); d1[, value := 0]; d1[obs, value := i.value, on = c('i', 'j')]
    d2 = data.table::copy(tmpl); d2[, value := 0]; d2[ex, value := i.value, on = c('i', 'j')]
    compdats(d1, d2, theta = 2, ifsum = TRUE, area0 = stats::median(width(gr))^2)
}
