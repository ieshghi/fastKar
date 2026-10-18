## Differential tests of R/comparisons.R (compdats, compmaps, rebin_matrix, make_noisydat,
## make_noisymap, hic_cov, kl_nb) against the reference fastKar. Inputs are simulated with the
## reference forward_simulate so they do not depend on the package under test.

.cmp_cache = new.env()
## reference maps for fixture `nm` at pix: list(a, b) = first two karyotypes (b = a if only one)
cmp_maps = function(nm, pix) {
    key = paste(nm, pix)
    if (is.null(.cmp_cache[[key]])) {
        fx = load_fixture(nm)
        sim = function(gw) suppressWarnings(ref('forward_simulate')(gw, target_region = fx$ft, pix.size = pix, depth = 10))
        a = sim(fx$gw[[1]])
        b = if (length(fx$gw) > 1) sim(fx$gw[[2]]) else a
        .cmp_cache[[key]] = list(a = a, b = b, fx = fx)
    }
    .cmp_cache[[key]]
}
noisy = function(m, seed = 11, nsamp = 1, theta = 2) { set.seed(seed); ref('make_noisymap')(m, nsamp, theta) }
## data.tables equal column by column (types, values, order), same key
expect_dt_identical = function(new, old, info = NULL) {
    expect_identical(names(new), names(old), info = info)
    expect_identical(data.table::key(new), data.table::key(old), info = info)
    expect_identical(as.list(new), as.list(old), info = info)
}
sparse_map = function(m, keep_frac = 0.5, seed = 5) {
    d = m$dat
    set.seed(seed)
    d = d[sort(sample.int(nrow(d), floor(nrow(d) * keep_frac)))]
    GxG::gM(gr = m$gr, dat = d[, .(i, j, value)])
}

small_cases = list(c('toy_dtrans', 1e5), c('toy_repdup', 1e4), c('toy_bfb', 1e4), c('ccle_pyrgo_local', 1e4), c('ccle_dm_local', 1e5))

test_that('compmaps matches the oracle (NB, Poisson, KL, if.diag, area0, mask, maskzeros, ifsum)', {
    skip_if_no_oracle()
    for (cs in small_cases) {
        mm = cmp_maps(cs[1], as.numeric(cs[2])); a = mm$a; b = mm$b
        obs = noisy(a)[[1]]
        mask = a$gr[ceiling(length(a$gr) / 3)] + 1
        args = list(list(), list(theta = 0), list(return_kl = TRUE), list(if.diag = FALSE),
                    list(area0 = 1e9), list(maskzeros = TRUE), list(mask = mask),
                    list(return_kl = TRUE, if.diag = FALSE, mask = mask), list(theta = 0, if.diag = FALSE))
        for (ar in args) for (ifsum in c(TRUE, FALSE)) {
            info = paste(cs[1], cs[2], paste(names(ar), collapse = ','), ifsum)
            new = do.call(opt('compmaps'), c(list(obs, b, ifsum = ifsum), ar))
            old = do.call(ref('compmaps'), c(list(obs, b, ifsum = ifsum), ar))
            if (ifsum) expect_equal(new, old, tolerance = EQUIV_TOL, info = info)
            else expect_gm_equiv(new, old, info = info)
        }
        ## expected map against itself, and noisy against noisy
        expect_equal(opt('compmaps')(a, a, ifsum = TRUE), ref('compmaps')(a, a, ifsum = TRUE), tolerance = EQUIV_TOL)
        expect_equal(opt('compmaps')(obs, obs, return_kl = TRUE, ifsum = TRUE), ref('compmaps')(obs, obs, return_kl = TRUE, ifsum = TRUE), tolerance = EQUIV_TOL)
    }
})

test_that('compmaps matches the oracle with missing pixels and with different bins', {
    skip_if_no_oracle()
    for (cs in small_cases[1:4]) {
        mm = cmp_maps(cs[1], as.numeric(cs[2])); a = mm$a; b = mm$b
        obs = noisy(a)[[1]]
        sobs = sparse_map(obs)            # observed with missing pixels
        sb = sparse_map(b, 0.7, seed = 9) # and a sparse expected map
        for (ifsum in c(TRUE, FALSE)) for (kl in c(FALSE, TRUE)) {
            info = paste(cs[1], ifsum, kl)
            new = opt('compmaps')(sobs, b, ifsum = ifsum, return_kl = kl); old = ref('compmaps')(sobs, b, ifsum = ifsum, return_kl = kl)
            if (ifsum) expect_equal(new, old, tolerance = EQUIV_TOL, info = info) else expect_gm_equiv(new, old, info = info)
            new = opt('compmaps')(obs, sb, ifsum = ifsum, return_kl = kl); old = ref('compmaps')(obs, sb, ifsum = ifsum, return_kl = kl)
            if (ifsum) expect_equal(new, old, tolerance = EQUIV_TOL, info = info) else expect_gm_equiv(new, old, info = info)
        }
    }
    ## observed at 1e5 vs simulation at 1e4 (rebinned through disjoin/agg)
    for (nm in c('toy_dtrans', 'toy_repdup')) {
        coarse = noisy(cmp_maps(nm, 1e5)$a)[[1]]
        fine = cmp_maps(nm, 1e4)$b
        for (ifsum in c(TRUE, FALSE)) {
            new = suppressWarnings(opt('compmaps')(coarse, fine, ifsum = ifsum))
            old = suppressWarnings(ref('compmaps')(coarse, fine, ifsum = ifsum))
            if (ifsum) expect_equal(new, old, tolerance = EQUIV_TOL, info = nm) else expect_gm_equiv(new, old, info = nm)
        }
    }
})

test_that('compdats matches the oracle, including its setkey side effect and template error', {
    skip_if_no_oracle()
    for (cs in small_cases) {
        mm = cmp_maps(cs[1], as.numeric(cs[2])); a = mm$a; b = mm$b
        set.seed(3); nd = ref('make_noisydat')(a, 2)[[2]]           # has widthprod
        gmd = noisy(a)[[1]]$dat                                     # gM table: no widthprod
        set.seed(4); perm = sample.int(nrow(nd))
        for (tst in list(nd, gmd)) for (ar in list(list(), list(theta = 0), list(return_kl = TRUE), list(if.diag = FALSE),
                                                   list(area0 = 1e9), list(return_kl = TRUE, if.diag = FALSE))) for (ifsum in c(TRUE, FALSE)) {
            info = paste(cs[1], 'widthprod' %in% names(tst), paste(names(ar), collapse = ','), ifsum)
            ## shuffled inputs: both implementations must sort the caller's tables the same way
            t1 = data.table::copy(tst)[perm]; d1 = data.table::copy(b$dat)[rev(seq_len(nrow(b$dat)))]
            t2 = data.table::copy(t1); d2 = data.table::copy(d1)
            new = do.call(opt('compdats'), c(list(t1, d1, ifsum = ifsum), ar))
            old = do.call(ref('compdats'), c(list(t2, d2, ifsum = ifsum), ar))
            if (ifsum) expect_equal(new, old, tolerance = EQUIV_TOL, info = info)
            else { expect_pix_equiv(new, old, info = info); expect_identical(names(new), names(old)); expect_identical(data.table::key(new), data.table::key(old), info = info) }
            expect_dt_identical(t1, t2, info = paste(info, 'dat_test side effect'))
            expect_dt_identical(d1, d2, info = paste(info, 'dat_true side effect'))
        }
        ## NA and zero expectations (value.y == 0 -> logprob 0, NA dropped by na.rm)
        d = data.table::copy(b$dat); d[c(1, 3), value := 0]; d[2, value := NA]
        for (ifsum in c(TRUE, FALSE)) {
            new = opt('compdats')(data.table::copy(nd), data.table::copy(d), ifsum = ifsum); old = ref('compdats')(data.table::copy(nd), data.table::copy(d), ifsum = ifsum)
            if (ifsum) expect_equal(new, old, tolerance = EQUIV_TOL) else expect_identical(as.list(new), as.list(old))
        }
        ## different templates: both refuse
        expect_error(opt('compdats')(data.table::copy(nd)[-1], data.table::copy(b$dat)), 'Equalize')
        expect_error(ref('compdats')(data.table::copy(nd)[-1], data.table::copy(b$dat)), 'Equalize')
    }
})

test_that('make_noisydat and make_noisymap give identical draws under the same seed', {
    skip_if_no_oracle()
    for (cs in small_cases) {
        a = cmp_maps(cs[1], as.numeric(cs[2]))$a
        for (src in list(a, sparse_map(a))) for (theta in c(2, 0)) for (nsamp in c(1, 7)) {
            info = paste(cs[1], theta, nsamp, nrow(src$dat))
            new = with_seed(42, function() opt('make_noisydat')(src, nsamp, theta)); r_new = runif(1)
            old = with_seed(42, function() ref('make_noisydat')(src, nsamp, theta)); r_old = runif(1)
            expect_identical(names(new), names(old), info = info)
            for (k in seq_along(old)) expect_dt_identical(new[[k]], old[[k]], info = paste(info, k))
            expect_identical(r_new, r_old, info = paste(info, 'RNG state after'))
        }
        new = with_seed(7, function() opt('make_noisymap')(a, 3, 2)); old = with_seed(7, function() ref('make_noisymap')(a, 3, 2))
        expect_identical(names(new), names(old))
        for (k in seq_along(old)) expect_gm_equiv(new[[k]], old[[k]], info = paste(cs[1], 'noisymap', k))
    }
    ## a gm.out = FALSE list map (bins as a data.table)
    fx = load_fixture('toy_dtrans')
    lm = ref('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = 1e5, depth = 10, gm.out = FALSE)
    new = with_seed(1, function() opt('make_noisydat')(lm, 2)); old = with_seed(1, function() ref('make_noisydat')(lm, 2))
    for (k in seq_along(old)) expect_dt_identical(new[[k]], old[[k]])
})

test_that('the fast paths are taken on simulated and noisy maps', {
    skip_if_no_oracle()
    for (cs in small_cases) {
        mm = cmp_maps(cs[1], as.numeric(cs[2]))
        obs = noisy(mm$a)[[1]]
        n = length(mm$a$gr)
        for (m in list(mm$a, mm$b, obs)) expect_true(opt('.is_template')(opt('.gm_rawdat')(m), n))
        expect_true(opt('.aligned_pair')(opt('.gm_rawdat')(obs), opt('.gm_rawdat')(mm$b)))
        expect_false(is.null(opt('.noisy_prep')(mm$a, 3)))
        expect_false(is.null(opt('.noisy_prep')(sparse_map(mm$a), 3)))
        expect_true(is.null(opt('.noisy_prep')(mm$a, 2.5)))   # non-integer nsamp -> reference code
        set.seed(2); nd = opt('make_noisydat')(mm$a, 1)[[1]]
        expect_true(opt('.aligned_pair')(nd, opt('.gm_rawdat')(mm$b)))
    }
})

test_that('rebin_matrix matches the oracle on both paths', {
    skip_if_no_oracle()
    for (nm in c('toy_dtrans', 'toy_repdup', 'ccle_pyrgo_local')) {
        a = cmp_maps(nm, 1e4)$a
        expect_identical(opt('rebin_matrix')(a, a$gr), a)                   # already on the bins: same object
        s = sparse_map(a)
        expect_gm_equiv(opt('rebin_matrix')(s, s$gr), ref('rebin_matrix')(s, s$gr), info = nm)   # fill missing pixels
        coarse = cmp_maps(nm, 1e5)$a
        expect_gm_equiv(suppressWarnings(opt('rebin_matrix')(a, coarse$gr)), suppressWarnings(ref('rebin_matrix')(a, coarse$gr)), info = nm)
        expect_gm_equiv(suppressWarnings(opt('rebin_matrix')(s, coarse$gr)), suppressWarnings(ref('rebin_matrix')(s, coarse$gr)), info = nm)
    }
})

test_that('hic_cov and kl_nb match the oracle', {
    skip_if_no_oracle()
    for (cs in small_cases) {
        a = cmp_maps(cs[1], as.numeric(cs[2]))$a
        for (m in list(a, sparse_map(a), noisy(a)[[1]])) {
            new = opt('hic_cov')(m); old = ref('hic_cov')(m)
            expect_identical(gUtils::gr2dt(new)[, -'cov'], gUtils::gr2dt(old)[, -'cov'], info = cs[1])
            expect_equal(new$cov, old$cov, tolerance = EQUIV_TOL, info = cs[1])
            expect_identical(GenomeInfoDb::seqlengths(new), GenomeInfoDb::seqlengths(old))
        }
        b = cmp_maps(cs[1], as.numeric(cs[2]))$b
        for (r in c(2, 0.5)) for (ifsum in c(TRUE, FALSE))
            expect_equal(opt('kl_nb')(a$dat$value, b$dat$value, r = r, ifsum = ifsum), ref('kl_nb')(a$dat$value, b$dat$value, r = r, ifsum = ifsum), tolerance = EQUIV_TOL)
    }
})
