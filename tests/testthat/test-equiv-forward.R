## Differential tests of the forward simulator against the pinned reference (oracle).
## Default set: every fixture at a cheap resolution over the argument grid; the heavy CCLE
## cases at 1e4 run only with FASTKAR_SLOW_TESTS=1 (see test-equiv-forward-slow.R).

fwd_grid = function(full = TRUE) {
    g = data.table::CJ(if.sum = c(TRUE, FALSE), gm.out = c(TRUE, FALSE), model = c(0, 1e5),
                       if.interchr = c(TRUE, FALSE), purity = c(1, 0.6, 0))
    g = g[!(if.sum == FALSE & purity < 1)] # reference: purity < 1 needs if.sum = TRUE (both error; checked separately)
    if (!full) g = g[gm.out == FALSE & (if.interchr == TRUE | model == 0)]
    g
}

run_fwd_grid = function(fx, gw, target, pix, grid, info) {
    for (r in seq_len(nrow(grid))) {
        a = grid[r]
        args = list(gw, target_region = target, pix.size = pix, depth = 10, if.sum = a$if.sum,
                    gm.out = a$gm.out, model = a$model, if.interchr = a$if.interchr, purity = a$purity)
        expect_both_equiv('forward_simulate', args, info = paste(info, pix, paste(names(a), unlist(a), sep = '=', collapse = ' ')))
    }
}

small_grid = data.table::data.table(if.sum = c(TRUE, TRUE, TRUE, FALSE), gm.out = c(TRUE, FALSE, FALSE, FALSE),
                                    model = c(0, 0, 1e5, 0), if.interchr = c(TRUE, FALSE, TRUE, TRUE), purity = c(1, 0.6, 1, 1))

test_that('forward_simulate matches the oracle on the toy fixtures (1e4/1e5/1e6, full grid, all targets)', {
    skip_if_no_oracle()
    for (nm in fixture_names(toy_only = TRUE)) {
        fx = load_fixture(nm)
        tg = fixture_targets(fx)
        ## full grid at 1e5 on the fixtures with walk CN / circular walks / revisited tiny nodes
        full = nm %in% c('toy_bfb', 'toy_circ', 'toy_tiny') || slow_tests()
        run_fwd_grid(fx, fx$gw[[1]], fx$ft, 1e5, fwd_grid(full = full), paste(nm, 'kar 1 footprint'))
        for (pix in c(1e4, 1e6)) run_fwd_grid(fx, fx$gw[[1]], fx$ft, pix, if (full) fwd_grid(full = FALSE) else small_grid, paste(nm, 'kar 1 footprint'))
        if (length(fx$gw) > 1) run_fwd_grid(fx, fx$gw[[2]], fx$ft, 1e5, small_grid, paste(nm, 'kar 2 footprint'))
        for (tn in setdiff(names(tg), 'footprint')) {
            run_fwd_grid(fx, fx$gw[[1]], tg[[tn]], 1e5, small_grid, paste(nm, 'kar 1', tn))
            run_fwd_grid(fx, fx$gw[[1]], tg[[tn]], 1e4, small_grid[2], paste(nm, 'kar 1', tn))
        }
    }
})

test_that('forward_simulate matches the oracle on the CCLE fixtures at 1e5 (and 1e4 for local graphs / pyrgo)', {
    skip_if_no_oracle()
    for (nm in grep('^ccle_', fixture_names(), value = TRUE)) {
        fx = load_fixture(nm)
        tg = fixture_targets(fx)
        heavy = grepl('^ccle_dm', nm) # ~4 s per reference call at 1e5
        run_fwd_grid(fx, fx$gw[[1]], fx$ft, 1e5, if (heavy) small_grid else fwd_grid(full = FALSE), paste(nm, 'kar 1 footprint'))
        run_fwd_grid(fx, fx$gw[[2]], fx$ft, 1e5, small_grid[if (heavy) 2 else 1:4], paste(nm, 'kar 2 footprint'))
        for (tn in setdiff(names(tg), 'footprint')) run_fwd_grid(fx, fx$gw[[1]], tg[[tn]], 1e5, small_grid[2], paste(nm, tn))
        if (nm %in% c('ccle_pyrgo', 'ccle_pyrgo_local', 'ccle_tyfonas_local')) {
            run_fwd_grid(fx, fx$gw[[1]], fx$ft, 1e4, small_grid, paste(nm, 'kar 1 footprint'))
        } else if (nm == 'ccle_dm_local') {
            run_fwd_grid(fx, fx$gw[[1]], fx$ft, 1e4, small_grid[2], paste(nm, 'kar 1 footprint'))
        }
    }
})

test_that('forward_simulate matches on walk-set edge cases: CN, minus strand, walks missing the target', {
    skip_if_no_oracle()
    fx = load_fixture('toy_circ') # circular walk on nodes 3:6 (cn 1) + linear walk on 1:10 (cn 2)
    gw = fx$gw[[1]]
    ## a window that only the linear walk reaches: the circular walk misses the target
    far = GenomicRanges::GRanges('1', IRanges::IRanges(1.4e6 + 1, 2e6))
    run_fwd_grid(fx, gw, far, 1e5, fwd_grid(full = TRUE)[purity == 1 | if.sum == TRUE], 'toy_circ far window')
    ## walks given as a plain list (no CN column), one reversed, plus a walk that misses the target
    lst = list(graph = fx$gg, snode.id = list(-rev(gw$snode.id[[2]]), gw$snode.id[[1]]), circular = c(FALSE, TRUE))
    run_fwd_grid(fx, lst, far, 1e5, fwd_grid(full = FALSE), 'toy_circ list walks')
    run_fwd_grid(fx, lst, fx$ft, 1e5, fwd_grid(full = FALSE), 'toy_circ list walks footprint')
    ## no walk hits the target: both error
    none = GenomicRanges::GRanges('1', IRanges::IRanges(1, 1e5))
    lst2 = list(graph = fx$gg, snode.id = list(gw$snode.id[[1]]), circular = TRUE)
    r = run_both('forward_simulate', list(lst2, target_region = none, pix.size = 1e5))
    expect_true(inherits(r$old, 'try-error'))
    expect_true(inherits(r$new, 'try-error'))
    ## purity < 1 with if.sum = FALSE errors in both
    r = run_both('forward_simulate', list(gw, target_region = fx$ft, pix.size = 1e5, if.sum = FALSE, purity = 0.5))
    expect_true(inherits(r$old, 'try-error') && inherits(r$new, 'try-error'))
    ## tiny 2 bp node revisited back to back on both strands, minus-strand steps
    fx = load_fixture('toy_tiny')
    run_fwd_grid(fx, fx$gw[[1]], fx$ft, 1e4, fwd_grid(full = TRUE), 'toy_tiny')
    lst = list(graph = fx$gg, snode.id = list(-rev(fx$gw[[1]]$snode.id[[1]])), circular = FALSE)
    run_fwd_grid(fx, lst, fx$ft, 1e4, fwd_grid(full = FALSE), 'toy_tiny reversed')
})

test_that('prep_for_sim matches the oracle (identical tilings and width tables)', {
    skip_if_no_oracle()
    for (nm in fixture_names()) {
        fx = load_fixture(nm)
        tg = fixture_targets(fx)
        pixes = if (grepl('^toy_', nm)) c(1e4, 1e5, 1e6) else c(1e4, 1e5)
        tns = if (grepl('^toy_', nm)) names(tg) else c('footprint', 'past')
        for (pix in pixes) for (tn in tns) {
            r = run_both('prep_for_sim', list(fx$gw[[1]], tg[[tn]], pix))
            info = paste(nm, pix, tn)
            if (inherits(r$old, 'try-error')) { expect_true(inherits(r$new, 'try-error'), info = info); next }
            expect_identical(r$new, r$old, info = info)
            expect_identical(data.table::key(r$new$tiled.target), data.table::key(r$old$tiled.target), info = info)
        }
        ## target_region = NULL falls back to the walks' footprint
        expect_identical(opt('prep_for_sim')(fx$gw[[1]], NULL, 1e5), ref('prep_for_sim')(fx$gw[[1]], NULL, 1e5), info = nm)
    }
})

test_that('simulate_diploid_map and make_template_dat match the oracle', {
    skip_if_no_oracle()
    for (nm in c(fixture_names(toy_only = TRUE), 'ccle_pyrgo', 'ccle_dm_local', 'ccle_chromothripsis')) {
        fx = load_fixture(nm)
        tg = fixture_targets(fx)
        for (tn in c('footprint', 'multi', 'past')) {
            pt = ref('prep_for_sim')(fx$gw[[1]], tg[[tn]], 1e5)$tiled.target
            for (ic in c(TRUE, FALSE)) for (gmo in c(TRUE, FALSE)) {
                expect_both_equiv('simulate_diploid_map', list(data.table::copy(pt), depth = 3.5, model = 0, if.interchr = ic, gm.out = gmo),
                                  info = paste(nm, tn, ic, gmo))
            }
            ## long-read model: the reference errors (data.table looks up same_chrom as a variable)
            expect_both_equiv('simulate_diploid_map', list(data.table::copy(pt), model = 1e5, gm.out = FALSE), info = paste(nm, tn, 'model 1e5'))
            ## the input table is not modified
            pt0 = data.table::copy(pt)
            invisible(opt('simulate_diploid_map')(pt, gm.out = FALSE))
            expect_identical(pt, pt0)
            ## overlapping target ranges give duplicate tiles at distance 0 off the diagonal
            if (tn == 'footprint') {
                dup = rbind(pt[1:min(4, .N)], pt[1:min(4, .N)])[, tile.id := .I]
                expect_both_equiv('simulate_diploid_map', list(dup, gm.out = FALSE), info = paste(nm, 'duplicate tiles'))
            }
            mt_new = opt('make_template_dat')(pt); mt_old = ref('make_template_dat')(pt)
            expect_identical(mt_new, mt_old, info = paste(nm, tn, 'template dt'))
            expect_identical(opt('make_template_dat')(dt2gr(pt)), ref('make_template_dat')(dt2gr(pt)), info = paste(nm, tn, 'template gr'))
            expect_identical(opt('make_template_dat_cpp')(pt), ref('make_template_dat_cpp')(pt), info = paste(nm, tn, 'template cpp'))
        }
    }
})

test_that('simulate_walks matches the oracle when called directly (as identifiability.R does)', {
    skip_if_no_oracle()
    for (nm in c('toy_bfb', 'toy_repdup', 'ccle_pyrgo', 'ccle_dm_local', 'ccle_tyfonas')) {
        fx = load_fixture(nm)
        p = ref('prep_for_sim')(fx$gw[[1]], fx$ft, 1e5)
        for (k in seq_len(min(if (nm == 'ccle_dm_local') 1 else 3, length(fx$gw)))) {
            for (s in c(TRUE, FALSE)) {
                expect_both_equiv('simulate_walks', list(fx$gw[[k]], p$tiled.target, p$widthdt, if.sum = s, depth = 5, gm.out = FALSE),
                                  info = paste(nm, k, s))
            }
        }
    }
})

test_that('sum_matrices and calculate_interchrom match the oracle', {
    skip_if_no_oracle()
    fx = load_fixture('toy_threeway')
    p = ref('prep_for_sim')(fx$gw[[1]], fx$ft, 1e5)
    maps = ref('simulate_walks')(fx$gw[[1]], p$tiled.target, p$widthdt, if.sum = FALSE, gm.out = FALSE)
    expect_fwd_equiv(opt('sum_matrices')(lapply(maps, data.table::copy)), ref('sum_matrices')(lapply(maps, data.table::copy)), info = 'sum_matrices')
    cn = rep(1, length(maps))
    expect_fwd_equiv(opt('calculate_interchrom')(lapply(maps, data.table::copy), NULL, p$tiled.target, cn),
                     ref('calculate_interchrom')(lapply(maps, data.table::copy), NULL, p$tiled.target, cn), info = 'calculate_interchrom')
    ## tables on different pixel sets take the general (merge) path
    part = lapply(maps, function(m) data.table::copy(m)[id %% 3 != 0])
    part[[1]] = part[[1]][id > 5]
    expect_fwd_equiv(opt('sum_matrices')(lapply(part, data.table::copy)), ref('sum_matrices')(lapply(part, data.table::copy)), info = 'sum_matrices sparse')
})
