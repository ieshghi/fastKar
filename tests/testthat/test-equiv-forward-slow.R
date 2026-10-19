## Slow differential tests (FASTKAR_SLOW_TESTS=1): the heavy CCLE cases at 1e4, where the
## reference takes 10-40 s per call, and a wider argument grid on the CCLE graphs at 1e5.

test_that('forward_simulate matches the oracle on every CCLE fixture at 1e4', {
    skip_if_no_oracle()
    if (!slow_tests()) skip('set FASTKAR_SLOW_TESTS=1 for the heavy cases')
    grid = data.table::data.table(if.sum = c(TRUE, TRUE, TRUE), gm.out = c(TRUE, FALSE, FALSE),
                                  model = c(0, 0, 1e5), if.interchr = c(TRUE, FALSE, TRUE), purity = c(1, 0.6, 1))
    for (nm in grep('^ccle_', fixture_names(), value = TRUE)) {
        fx = load_fixture(nm)
        for (k in 1:2) {
            for (r in seq_len(nrow(grid))) {
                a = grid[r]
                args = list(fx$gw[[k]], target_region = fx$ft, pix.size = 1e4, depth = 10, if.sum = a$if.sum,
                            gm.out = a$gm.out, model = a$model, if.interchr = a$if.interchr, purity = a$purity)
                expect_both_equiv('forward_simulate', args, info = paste(nm, k, 1e4, paste(unlist(a), collapse = ',')))
                if (k == 2) break
            }
        }
        ## per-walk maps on the local graphs (one table per walk)
        if (grepl('_local$', nm) && nm != 'ccle_chromothripsis_local')
            expect_both_equiv('forward_simulate', list(fx$gw[[1]], target_region = fx$ft, pix.size = 1e4, if.sum = FALSE, gm.out = FALSE),
                              info = paste(nm, 'per walk 1e4'))
    }
})

test_that('forward_simulate matches the oracle on the CCLE fixtures at 1e5 over the full grid', {
    skip_if_no_oracle()
    if (!slow_tests()) skip('set FASTKAR_SLOW_TESTS=1 for the heavy cases')
    g = data.table::CJ(if.sum = c(TRUE, FALSE), gm.out = c(TRUE, FALSE), model = c(0, 1e5),
                       if.interchr = c(TRUE, FALSE), purity = c(1, 0.6, 0))[!(if.sum == FALSE & purity < 1)]
    for (nm in grep('^ccle_', fixture_names(), value = TRUE)) {
        fx = load_fixture(nm)
        tg = fixture_targets(fx)
        for (tn in names(tg)) for (r in seq_len(nrow(g))) {
            a = g[r]
            if (tn != 'footprint' && a$gm.out) next
            args = list(fx$gw[[3]], target_region = tg[[tn]], pix.size = 1e5, depth = 10, if.sum = a$if.sum,
                        gm.out = a$gm.out, model = a$model, if.interchr = a$if.interchr, purity = a$purity)
            expect_both_equiv('forward_simulate', args, info = paste(nm, tn, paste(unlist(a), collapse = ',')))
        }
    }
})
