## Differential tests of R/comparisons.R on one large map: ccle_chromothripsis_local at 10 kb
## (4131 tiles, 8.5 M pixels). One reference simulation (~50 s) is shared by all cases; the
## second map is a deterministic perturbation of the first. Slower cases need
## FASTKAR_SLOW_TESTS=1.

large_maps = local({
    cache = NULL
    function() {
        if (is.null(cache)) {
            fx = load_fixture('ccle_chromothripsis_local')
            a = suppressWarnings(ref('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = 1e4, depth = 10))
            d = a$dat
            b = GxG::gM(gr = a$gr, dat = d[, .(i, j, value = value * (1 + 0.3 * sin(i + 2 * j)))])
            set.seed(99); obs = ref('make_noisymap')(a, 1)[[1]]
            cache <<- list(fx = fx, a = a, b = b, obs = obs)
        }
        cache
    }
})
slow = Sys.getenv('FASTKAR_SLOW_TESTS') == '1'

test_that('compmaps and compdats match the oracle on a large map', {
    skip_if_no_oracle()
    m = large_maps()
    expect_gt(length(m$a$gr), 4000)
    for (ar in list(list(), list(return_kl = TRUE), if (slow) list(theta = 0, if.diag = FALSE)))
        if (!is.null(ar)) expect_equal(do.call(opt('compmaps'), c(list(m$obs, m$b, ifsum = TRUE), ar)),
                                       do.call(ref('compmaps'), c(list(m$obs, m$b, ifsum = TRUE), ar)), tolerance = EQUIV_TOL)
    expect_gm_equiv(opt('compmaps')(m$obs, m$b), ref('compmaps')(m$obs, m$b))
    t1 = m$obs$dat; d1 = m$b$dat; t2 = data.table::copy(t1); d2 = data.table::copy(d1)
    expect_equal(opt('compdats')(t1, d1), ref('compdats')(t2, d2), tolerance = EQUIV_TOL)
    expect_identical(as.list(t1), as.list(t2)); expect_identical(data.table::key(t1), data.table::key(t2))
    expect_identical(data.table::key(d1), data.table::key(d2))
})

test_that('make_noisydat, hic_cov and rebin_matrix match the oracle on a large map', {
    skip_if_no_oracle()
    m = large_maps()
    new = with_seed(8, function() opt('make_noisydat')(m$a, 2)); old = with_seed(8, function() ref('make_noisydat')(m$a, 2))
    expect_identical(names(new), names(old))
    for (k in seq_along(old)) { expect_identical(as.list(new[[k]]), as.list(old[[k]])); expect_identical(names(new[[k]]), names(old[[k]])) }
    rm(new, old); invisible(gc())
    hn = opt('hic_cov')(m$a); ho = ref('hic_cov')(m$a)
    expect_identical(gUtils::gr2dt(hn)[, -'cov'], gUtils::gr2dt(ho)[, -'cov']); expect_equal(hn$cov, ho$cov, tolerance = EQUIV_TOL)
    expect_identical(opt('rebin_matrix')(m$a, m$a$gr), m$a)
    if (slow) {
        coarse = suppressWarnings(ref('forward_simulate')(m$fx$gw[[1]], target_region = m$fx$ft, pix.size = 1e5, depth = 10))
        expect_gm_equiv(suppressWarnings(opt('rebin_matrix')(m$a, coarse$gr)), suppressWarnings(ref('rebin_matrix')(m$a, coarse$gr)))
        expect_equal(suppressWarnings(opt('compmaps')(coarse, m$b, ifsum = TRUE)), suppressWarnings(ref('compmaps')(coarse, m$b, ifsum = TRUE)), tolerance = EQUIV_TOL)
    }
})

test_that('liktest_separable matches the oracle on a large map (slow)', {
    skip_if_no_oracle()
    if (!slow) skip('set FASTKAR_SLOW_TESTS=1')
    m = large_maps()
    new = with_seed(4, function() opt('liktest_separable')(m$a, m$b, nsamp = 2))
    old = with_seed(4, function() ref('liktest_separable')(m$a, m$b, nsamp = 2))
    expect_identical(new, old)
})
