## Property tests of R/comparisons.R on the package under test alone (no oracle): they catch
## errors the optimized and the reference code could share.

prop_maps = local({
    cache = list()
    function(nm, pix = 1e5) {
        key = paste(nm, pix)
        if (is.null(cache[[key]])) {
            fx = load_fixture(nm)
            sim = function(gw) suppressWarnings(opt('forward_simulate')(gw, target_region = fx$ft, pix.size = pix, depth = 10))
            cache[[key]] <<- list(a = sim(fx$gw[[1]]), b = sim(fx$gw[[2]]))
        }
        cache[[key]]
    }
})
scaled = function(m, f) GxG::gM(gr = m$gr, dat = m$dat[, .(i, j, value = value * f)])

test_that('expected NLL is minimized at the true mean', {
    for (nm in c('toy_dtrans', 'ccle_dm_local')) {
        m = prop_maps(nm); mu = m$a
        for (theta in c(2, 0)) {
            set.seed(1); obs = opt('make_noisydat')(mu, 30, theta)
            mean_nll = function(cand) mean(vapply(obs, function(o) opt('compdats')(data.table::copy(o), cand$dat, theta = theta), 1))
            at_truth = mean_nll(mu)
            for (f in c(0.7, 0.9, 1.1, 1.4)) expect_lt(at_truth, mean_nll(scaled(mu, f)))
            expect_lt(at_truth, mean_nll(m$b))
            ## and through compmaps, on a gMatrix observation
            o = opt('make_noisymap')(mu, 1, theta)[[1]]
            expect_lt(opt('compmaps')(o, mu, theta = theta, ifsum = TRUE), opt('compmaps')(o, scaled(mu, 2), theta = theta, ifsum = TRUE))
        }
    }
})

test_that('KL is non-negative and zero for identical maps', {
    for (nm in c('toy_dtrans', 'toy_threeway', 'ccle_dm_local')) {
        m = prop_maps(nm)
        kl = opt('compmaps')(m$a, m$b, return_kl = TRUE)
        expect_true(all(kl$dat$value >= 0))
        expect_gte(opt('compmaps')(m$a, m$b, return_kl = TRUE, ifsum = TRUE), 0)
        expect_identical(opt('compmaps')(m$a, m$a, return_kl = TRUE, ifsum = TRUE), 0)
        expect_identical(opt('compmaps')(m$b, m$b, return_kl = TRUE, ifsum = TRUE, if.diag = FALSE), 0)
        expect_identical(opt('kl_nb')(m$a$dat$value, m$a$dat$value), 0)
        expect_true(all(opt('kl_nb')(m$a$dat$value, m$b$dat$value, ifsum = FALSE) >= 0))
    }
})

test_that('compmaps and compdats do not depend on the row order of the observed table', {
    for (nm in c('toy_dtrans', 'ccle_dm_local')) {
        m = prop_maps(nm)
        set.seed(2); obs = opt('make_noisymap')(m$a, 1)[[1]]
        set.seed(3); d = obs$dat[sample.int(nrow(obs$dat))]
        shuffled = list(gr = obs$gr, dat = d)
        for (ar in list(list(), list(return_kl = TRUE), list(theta = 0, if.diag = FALSE))) {
            expect_identical(do.call(opt('compmaps'), c(list(shuffled, m$b, ifsum = TRUE), ar)),
                             do.call(opt('compmaps'), c(list(obs, m$b, ifsum = TRUE), ar)))
            expect_identical(do.call(opt('compdats'), c(list(data.table::copy(d), m$b$dat), ar)),
                             do.call(opt('compdats'), c(list(obs$dat, m$b$dat), ar)))
        }
        ## ifsum = TRUE is the NA-dropping sum of the per-pixel values
        expect_equal(opt('compdats')(obs$dat, m$b$dat), sum(opt('compdats')(obs$dat, m$b$dat, ifsum = FALSE)$value, na.rm = TRUE), tolerance = 1e-14)
    }
})

test_that('make_noisydat draws are reproducible and centred on the map', {
    m = prop_maps('ccle_dm_local')
    set.seed(5); x = opt('make_noisydat')(m$a, 3)
    set.seed(5); y = opt('make_noisydat')(m$a, 3)
    expect_identical(x, y)
    expect_false(identical(x[[1]]$value, x[[2]]$value))
    mu = m$a$dat$value
    sd_tot = sqrt(sum(mu + mu^2 / 2))   # NB(size = 2) variance mu + mu^2 / size
    for (s in x) expect_lt(abs(sum(s$value) - sum(mu)), 5 * sd_tot)
})
