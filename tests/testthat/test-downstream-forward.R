## T5: downstream end-to-end regression against the oracle on two CCLE local fixtures, fixed seeds.
## Candidate order must be identical and NLLs agree to 1e-10.

downstream_cases = c('ccle_pyrgo_local', 'ccle_dm_local')
NLL_TOL = 1e-10

## a minimal identifiability context: the fixture's local graph is already cut to its footprint, so
## use it directly with the footprint as target and identifiability_context()'s resolution rule
## (identifiability_context() itself needs gGnome's gGraph %&% method, not available in every session)
downstream_ctx = function(fx) {
    ft = gUtils::streduce(fx$ft)
    target_width = sum(as.numeric(GenomicRanges::width(ft)))
    list(graph = fx$gg, target_region = ft, pix.size = 1e3 * ceiling(target_width / (1e3 * 100)),
         target_width = target_width)
}

expect_kars_identical = function(new, old, info) {
    expect_identical(length(new), length(old), info = info)
    expect_identical(lapply(new, function(k) unname(k$snode.id)), lapply(old, function(k) unname(k$snode.id)), info = info)
    expect_identical(lapply(new, function(k) unname(k$circular)), lapply(old, function(k) unname(k$circular)), info = info)
}

test_that('bestfit_search(topk = -1) matches the oracle, with and without rescale', {
    skip_if_no_oracle()
    for (nm in downstream_cases) {
        fx = load_fixture(nm)
        ctx = downstream_ctx(fx)
        set.seed(11)
        hic = ref('make_noisymap')(ref('forward_simulate')(fx$gw[[1]], target_region = ctx$target_region, pix.size = ctx$pix.size, depth = 10), 1)[[1]]
        nsample = if (slow_tests()) 12 else 6
        combos = list(c(FALSE, 1), c(TRUE, 1), c(TRUE, 0.7))
        if (slow_tests()) combos = c(combos, list(c(FALSE, 0.7)))
        for (cb in combos) {
            rescale = as.logical(cb[1]); purity = cb[2]
            args = list(fx$gg, hic, nsample = nsample, depth = 10, pix.size = ctx$pix.size, target_region = ctx$target_region,
                        topk = -1, purity = purity, rescale = rescale)
            set.seed(5); old = do.call(ref('bestfit_search'), args); seed_old = .Random.seed
            set.seed(5); new = do.call(opt('bestfit_search'), args); seed_new = .Random.seed
            info = paste(nm, 'rescale', rescale, 'purity', purity)
            expect_kars_identical(new$kars, old$kars, info)
            expect_equal(new$nlls, old$nlls, tolerance = NLL_TOL, info = info)
            expect_identical(order(new$nlls), order(old$nlls), info = info)
            expect_equal(new$ratio, old$ratio, tolerance = NLL_TOL, info = info)
            expect_identical(new$purity, old$purity)
            expect_identical(length(new$sims), length(old$sims))
            for (k in seq_along(old$sims)) expect_gm_equiv(new$sims[[k]], old$sims[[k]], info = paste(info, 'sim', k))
            ## the random stream after the call is the same too
            expect_identical(seed_new, seed_old, info = info)
        }
        ## topk > 0 and default pix.size / target_region
        set.seed(6); old = ref('bestfit_search')(fx$gg, hic, nsample, 10, topk = 3)
        set.seed(6); new = opt('bestfit_search')(fx$gg, hic, nsample, 10, topk = 3)
        expect_kars_identical(new$kars, old$kars, paste(nm, 'topk 3'))
        expect_equal(new$nlls, old$nlls, tolerance = NLL_TOL)
    }
})

test_that('identifiability_hic_nlls and local_hic_nll_matrix match the oracle', {
    skip_if_no_oracle()
    for (nm in downstream_cases) {
        fx = load_fixture(nm)
        ctx = downstream_ctx(fx)
        set.seed(21)
        truth = ref('make_identifiability_truth')(ctx, depth = 10)
        set.seed(22)
        kars = c(list(truth$true_gw), ref('sample.gwalks')(ctx$graph, if (slow_tests()) 10 else 5, verbose = FALSE))
        old = ref('identifiability_hic_nlls')(kars, truth, ctx, depth = 10)
        new = opt('identifiability_hic_nlls')(kars, truth, ctx, depth = 10)
        expect_equal(new, old, tolerance = NLL_TOL, info = nm)
        expect_identical(order(new), order(old), info = nm)
        set.seed(23); old = ref('local_hic_nll_matrix')(kars, ctx, n_obs = 5, depth = 10)
        set.seed(23); new = opt('local_hic_nll_matrix')(kars, ctx, n_obs = 5, depth = 10)
        expect_equal(new, old, tolerance = NLL_TOL, info = nm)
        expect_identical(apply(new, 1, order), apply(old, 1, order), info = nm)
    }
})
