## Differential tests of the seeded comparison functions in R/comparisons.R: same seed -> same
## draws and the same results as the reference fastKar.

sim_ref = function(fx, gw, pix) suppressWarnings(ref('forward_simulate')(gw, target_region = fx$ft, pix.size = pix, depth = 10))
both_seeded = function(seed, name, ...) {
    new = with_seed(seed, function() opt(name)(...)); r_new = runif(1)
    old = with_seed(seed, function() ref(name)(...)); r_old = runif(1)
    list(new = new, old = old, rng_same = identical(r_new, r_old))
}

test_that('liktest_separable matches the oracle under the same seed', {
    skip_if_no_oracle()
    for (cs in list(c('toy_dtrans', 1e5), c('toy_threeway', 1e5), c('ccle_pyrgo_local', 1e4), c('ccle_dm_local', 1e5))) {
        fx = load_fixture(cs[1]); pix = as.numeric(cs[2])
        a = sim_ref(fx, fx$gw[[1]], pix); b = sim_ref(fx, fx$gw[[2]], pix)
        for (theta in c(2, 0)) for (nsamp in c(1, 6)) {
            r = both_seeded(17, 'liktest_separable', a, b, theta = theta, nsamp = nsamp)
            expect_identical(r$new, r$old, info = paste(cs[1], theta, nsamp))
            expect_true(r$rng_same)
        }
        ## depth-scaled copy: maps close enough that the test sometimes errs
        b2 = GxG::gM(gr = a$gr, dat = a$dat[, .(i, j, value = value * 1.02)])
        r = both_seeded(5, 'liktest_separable', a, b2, nsamp = 10)
        expect_identical(r$new, r$old, info = cs[1])
        ## sparse map: falls back to the reference path (compdats refuses unequal templates)
        sb = GxG::gM(gr = b$gr, dat = b$dat[seq(1, .N, 2)])
        expect_error(with_seed(1, function() opt('liktest_separable')(a, sb, nsamp = 2)), 'Equalize')
        expect_error(with_seed(1, function() ref('liktest_separable')(a, sb, nsamp = 2)), 'Equalize')
    }
})

test_that('longread_probdist matches the oracle (use.nodes, minsize)', {
    skip_if_no_oracle()
    for (nm in c('toy_dtrans', 'toy_threeway', 'toy_circ', 'toy_tiny', 'ccle_pyrgo_local')) {
        fx = load_fixture(nm)
        gws = fx$gw[seq_len(min(3, length(fx$gw)))]
        nodes = fx$gg$nodes$dt$node.id
        for (readL in c(5e4, 3e5)) for (ar in list(list(), list(minsize = 5e4), list(use.nodes = nodes[seq(1, length(nodes), 2)]),
                                                    list(use.nodes = nodes[1], minsize = 1e3), list(background = 1e-5))) {
            info = paste(nm, readL, paste(names(ar), collapse = ','))
            new = do.call(opt('longread_probdist'), c(list(gws, readL), ar))
            old = do.call(ref('longread_probdist'), c(list(gws, readL), ar))
            expect_identical(names(new), names(old), info = info)
            for (k in seq_along(old)) {
                expect_identical(names(new[[k]]), names(old[[k]]), info = info)
                expect_equal(unname(new[[k]]), unname(old[[k]]), tolerance = EQUIV_TOL, info = info)
            }
        }
        expect_equal(opt('longread_probdist')(gws[[1]], 1e5), ref('longread_probdist')(gws[[1]], 1e5), tolerance = EQUIV_TOL)
    }
})

test_that('liktest_separable_lr matches the oracle under the same seed', {
    skip_if_no_oracle()
    for (nm in c('toy_dtrans', 'toy_threeway', 'toy_repdup', 'ccle_pyrgo_local')) {
        fx = load_fixture(nm)
        nodes = fx$gg$nodes$dt$node.id
        for (readL in c(1e5, 1e6)) for (un in list(NULL, nodes[seq(1, length(nodes), 2)], nodes[1])) {
            info = paste(nm, readL, length(un))
            r = both_seeded(23, 'liktest_separable_lr', fx$gw[[1]], fx$gw[[2]], readL = readL, depth = 10, nsamp = 8, use.nodes = un, return.all = TRUE)
            expect_identical(r$new, r$old, info = info)
            expect_true(r$rng_same, info = info)
            r = both_seeded(23, 'liktest_separable_lr', fx$gw[[1]], fx$gw[[2]], readL = readL, depth = 0.5, nsamp = 3, use.nodes = un)
            expect_identical(r$new, r$old, info = info)
            new = opt('liktest_separable_lr')(fx$gw[[1]], fx$gw[[2]], readL = readL, depth = 10, use.nodes = un, return.kl = TRUE)
            old = ref('liktest_separable_lr')(fx$gw[[1]], fx$gw[[2]], readL = readL, depth = 10, use.nodes = un, return.kl = TRUE)
            expect_equal(new, old, tolerance = EQUIV_TOL, info = info)
        }
    }
})

test_that('group_klcomp matches the oracle (symmetric and two-set)', {
    skip_if_no_oracle()
    fx = load_fixture('toy_threeway')
    gws = c(fx$gw[1:3], fx$gw[1])   # repeated walk -> equal hashes -> 0 without simulating
    new = opt('group_klcomp')(gws, hic_resolution = 2e5); old = ref('group_klcomp')(gws, hic_resolution = 2e5)
    expect_identical(new[, .(i, j)], old[, .(i, j)]); expect_equal(new$kldiv, old$kldiv, tolerance = EQUIV_TOL)
    new = opt('group_klcomp')(fx$gw[1:2], fx$gw[2:4], target_region = fx$ft, hic_resolution = 1e5)
    old = ref('group_klcomp')(fx$gw[1:2], fx$gw[2:4], target_region = fx$ft, hic_resolution = 1e5)
    expect_identical(names(new), names(old))
    expect_identical(new[, .(i, j)], old[, .(i, j)]); expect_equal(new$kldiv, old$kldiv, tolerance = EQUIV_TOL)
    ## default resolution, sub-window target
    fx = load_fixture('toy_repdup')
    sub = fx$ft[1]
    new = opt('group_klcomp')(fx$gw[1:3], target_region = sub); old = ref('group_klcomp')(fx$gw[1:3], target_region = sub)
    expect_identical(new[, .(i, j)], old[, .(i, j)]); expect_equal(new$kldiv, old$kldiv, tolerance = EQUIV_TOL)
})

test_that('neighbor_separable matches the oracle under the same seed', {
    skip_if_no_oracle()
    ## neighbor_separable restricts each walk with an unqualified `gWalk %&% footprint`. In this
    ## environment that cannot run for either package: the search path finds a copy of the generic
    ## without the gWalk method, and gGnome's own gWalk method fails inside aggregate.formula()
    ## under current R. The footprint used here is the whole graph, so every walk is kept; shim
    ## `%&%` as the identity on gWalks (the global env precedes the search path in the lookup).
    assign('%&%', function(x, y) if (inherits(x, 'gWalk')) x else gUtils::`%&%`(x, y), envir = globalenv())
    withr::defer(rm('%&%', envir = globalenv()))
    for (nm in c('toy_dtrans', 'toy_threeway')) {
        fx = load_fixture(nm)
        nodes = fx$gg$nodes$dt$node.id
        r = suppressWarnings(both_seeded(31, 'neighbor_separable', fx$gg, n_init = 2, n_neighbor = 3, use.nodes = nodes,
                                         sep_samples = 4, resolution = c(2e5, -1e5, 5e5), depth = 10))
        expect_identical(r$new, r$old, info = nm)
        expect_true(r$rng_same, info = nm)
    }
    ## quirk kept: frozen.nodes is undefined when use.nodes is NULL
    fx = load_fixture('toy_dtrans')
    expect_error(opt('neighbor_separable')(fx$gg, 1, 1, resolution = 1e5), 'frozen.nodes')
    expect_error(ref('neighbor_separable')(fx$gg, 1, 1, resolution = 1e5), 'frozen.nodes')
})

test_that('f1score_comparemaps matches the oracle without a mask', {
    skip_if_no_oracle()
    for (cs in list(c('toy_dtrans', 1e5), c('ccle_pyrgo_local', 1e4))) {
        fx = load_fixture(cs[1]); pix = as.numeric(cs[2])
        a = sim_ref(fx, fx$gw[[1]], pix); b = sim_ref(fx, fx$gw[[2]], pix)
        b2 = GxG::gM(gr = a$gr, dat = a$dat[, .(i, j, value = value * 1.03)])
        for (hyp in list(b, b2)) for (theta in c(2, 0)) for (rs in c(FALSE, TRUE)) {
            r = both_seeded(3, 'f1score_comparemaps', a, hyp, theta = theta, significance = 0.25, return_samples = rs)
            expect_identical(r$new, r$old, info = paste(cs[1], theta, rs))
            expect_true(r$rng_same)
        }
    }
    ## quirk kept: the mask branch refers to an undefined `dat`
    fx = load_fixture('toy_dtrans'); a = sim_ref(fx, fx$gw[[1]], 1e5)
    ## (`dat` resolves to whatever the search path offers; either way both implementations fail alike)
    e_new = tryCatch(opt('f1score_comparemaps')(a, a, mask = a$gr[2]), error = conditionMessage)
    e_old = tryCatch(ref('f1score_comparemaps')(a, a, mask = a$gr[2]), error = conditionMessage)
    expect_type(e_old, 'character')
    expect_identical(e_new, e_old)
})
