## Property tests of the forward simulator (oracle-independent): they catch errors the
## reference and the optimized code could share.

fs = function(...) opt('forward_simulate')(...)
vals = function(x) { d = as_pixdt(x); d[order(id)]$value }

prop_cases = c('toy_bfb', 'toy_circ', 'toy_repdup', 'toy_tiny', 'ccle_pyrgo', 'ccle_dm_local', 'ccle_tyfonas_local')
prop_pix = function(nm) if (grepl('^toy_', nm)) 1e4 else 1e5

## the walks of a gWalk as a plain list
as_walk_list = function(gw, snode.id = gw$snode.id, cn = gw$dt$cn) {
    out = list(graph = gw$graph, snode.id = unname(snode.id), circular = unname(gw$circular))
    if (!is.null(cn)) out$dt = data.table::data.table(cn = cn)
    out
}

test_that('maps are non-negative and finite (Hi-C and long-read, summed and per walk)', {
    for (nm in prop_cases) {
        fx = load_fixture(nm)
        for (model in c(0, 1e5)) {
            v = vals(fs(fx$gw[[1]], target_region = fx$ft, pix.size = prop_pix(nm), depth = 10, model = model, gm.out = FALSE))
            expect_true(all(is.finite(v)) && all(v >= 0), info = paste(nm, model))
            expect_gt(sum(v), 0)
            per = fs(fx$gw[[1]], target_region = fx$ft, pix.size = prop_pix(nm), if.sum = FALSE, model = model, gm.out = FALSE)
            expect_true(all(vapply(per, function(d) all(is.finite(d$value)) && all(d$value >= 0), logical(1))), info = paste(nm, model))
        }
    }
})

test_that('maps are linear in depth and in walk CN', {
    for (nm in prop_cases) {
        fx = load_fixture(nm)
        gw = fx$gw[[1]]
        pix = prop_pix(nm)
        a = vals(fs(gw, target_region = fx$ft, pix.size = pix, depth = 10, gm.out = FALSE))
        b = vals(fs(gw, target_region = fx$ft, pix.size = pix, depth = 25, gm.out = FALSE))
        expect_equal(b, 2.5 * a, tolerance = 1e-13, info = nm)
        ## walk CN: without the interchromosomal term (which depends on node CN only) the map scales
        cn0 = if (is.null(gw$dt$cn)) rep(1, length(gw)) else gw$dt$cn
        m1 = vals(fs(as_walk_list(gw, cn = cn0), target_region = fx$ft, pix.size = pix, if.interchr = FALSE, gm.out = FALSE))
        m3 = vals(fs(as_walk_list(gw, cn = 3 * cn0), target_region = fx$ft, pix.size = pix, if.interchr = FALSE, gm.out = FALSE))
        expect_equal(m3, 3 * m1, tolerance = 1e-13, info = nm)
        ## doubling every walk equals CN 2
        lst = as_walk_list(gw, cn = cn0)
        dbl = list(graph = lst$graph, snode.id = c(lst$snode.id, lst$snode.id), circular = c(lst$circular, lst$circular),
                   dt = data.table::data.table(cn = c(cn0, cn0)))
        m2 = vals(fs(as_walk_list(gw, cn = 2 * cn0), target_region = fx$ft, pix.size = pix, if.interchr = FALSE, gm.out = FALSE))
        expect_equal(vals(fs(dbl, target_region = fx$ft, pix.size = pix, if.interchr = FALSE, gm.out = FALSE)), m2, tolerance = 1e-13, info = nm)
    }
})

test_that('per-walk maps plus the interchromosomal term sum to the if.sum = TRUE map', {
    for (nm in prop_cases) {
        fx = load_fixture(nm)
        pix = prop_pix(nm)
        depth = 7
        total = fs(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = depth, gm.out = FALSE)
        per = fs(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = depth, if.sum = FALSE, gm.out = FALSE)
        tt = total$gr
        tm = make_template_dat(tt)
        inter = data.table(i = tm$i, j = tm$j, value = depth / 2 * tm$widthprod * fastKar::interchrom_density * tt$cn[tm$i] * tt$cn[tm$j])
        inter = double_samelocus_pixels(inter, tt)$value # toy_circ has overlapping nodes
        s =Reduce(`+`, lapply(per, function(d) d[order(id)]$value))
        expect_equal(s + inter, vals(total), tolerance = 1e-12, info = nm)
        intra = vals(fs(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = depth, if.interchr = FALSE, gm.out = FALSE))
        expect_equal(s, intra, tolerance = 1e-12, info = nm)
    }
})

test_that('reverse-complemented walks and rotated circular walks give the same map', {
    for (nm in prop_cases) {
        fx = load_fixture(nm)
        gw = fx$gw[[1]]
        pix = prop_pix(nm)
        base = vals(fs(as_walk_list(gw), target_region = fx$ft, pix.size = pix, gm.out = FALSE))
        rc = lapply(gw$snode.id, function(w) -rev(w))
        expect_equal(vals(fs(as_walk_list(gw, snode.id = rc), target_region = fx$ft, pix.size = pix, gm.out = FALSE)), base,
                     tolerance = 1e-12, info = paste(nm, 'revcomp'))
        if (any(gw$circular)) {
            rot = lapply(seq_along(gw$snode.id), function(k) {
                w = gw$snode.id[[k]]
                if (gw$circular[k] && length(w) > 1) { s = ceiling(length(w) / 2); c(w[-(1:s)], w[1:s]) } else w
            })
            expect_equal(vals(fs(as_walk_list(gw, snode.id = rot), target_region = fx$ft, pix.size = pix, gm.out = FALSE)), base,
                         tolerance = 1e-12, info = paste(nm, 'rotation'))
        }
    }
    ## a hand-made circular walk: every rotation, both orientations
    fx = load_fixture('toy_circ')
    gw = fx$gw[[1]]
    w = gw$snode.id[[1]]
    base = vals(fs(gw, target_region = fx$ft, pix.size = 1e5, gm.out = FALSE))
    for (s in seq_along(w)[-1]) {
        rw = c(w[s:length(w)], w[seq_len(s - 1)])
        for (orient in c(1, -1)) {
            rr = if (orient == 1) rw else -rev(rw)
            v = vals(fs(as_walk_list(gw, snode.id = list(rr, gw$snode.id[[2]])), target_region = fx$ft, pix.size = 1e5, gm.out = FALSE))
            expect_equal(v, base, tolerance = 1e-12, info = paste('toy_circ rotation', s, orient))
        }
    }
})

test_that('purity mixing is linear', {
    for (nm in prop_cases) {
        fx = load_fixture(nm)
        pix = prop_pix(nm)
        m = function(p) vals(fs(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10, purity = p, gm.out = FALSE))
        m1 = m(1); m0 = m(0)
        for (p in c(0.25, 0.6)) expect_equal(m(p), p * m1 + (1 - p) * m0, tolerance = 1e-12, info = paste(nm, p))
        ## purity 0 is the diploid reference on the same tiling
        tt = opt('prep_for_sim')(fx$gw[[1]], fx$ft, pix)$tiled.target
        dip = opt('simulate_diploid_map')(tt, depth = 5, gm.out = FALSE)$dat
        expect_equal(m0, dip[order(id)]$value, tolerance = 1e-12, info = nm)
    }
})

test_that('results do not depend on mc.cores', {
    for (nm in c('toy_bfb', 'ccle_dm_local', 'ccle_tyfonas')) {
        fx = load_fixture(nm)
        for (s in c(TRUE, FALSE)) {
            a = fs(fx$gw[[1]], target_region = fx$ft, pix.size = 1e5, if.sum = s, gm.out = FALSE, mc.cores = 1)
            b = fs(fx$gw[[1]], target_region = fx$ft, pix.size = 1e5, if.sum = s, gm.out = FALSE, mc.cores = 4)
            expect_identical(a, b, info = paste(nm, s))
        }
    }
})

test_that('a reused prep_for_sim() result gives the same output', {
    for (nm in c('toy_bfb', 'ccle_pyrgo', 'ccle_dm_local')) {
        fx = load_fixture(nm)
        pp = opt('prep_for_sim')(fx$gw[[1]], fx$ft, 1e5)
        for (k in seq_len(min(3, length(fx$gw)))) {
            for (p in c(1, 0.6)) {
                a = fs(fx$gw[[k]], target_region = fx$ft, pix.size = 1e5, purity = p, gm.out = FALSE)
                b = fs(fx$gw[[k]], target_region = fx$ft, pix.size = 1e5, purity = p, gm.out = FALSE, prepped = pp)
                expect_identical(a, b, info = paste(nm, k, p))
            }
        }
    }
})

test_that('a target with gaps gives the same pixels as the full target (gaps inside nodes keep their length)', {
    ## the sub-target is every other block of 3 tiles of the full tiling, so many gaps fall inside nodes; walk_rows()
    ## used to drop the untiled part of such nodes, which made tiles on either side of a gap neighbours
    for (nm in c('toy_bfb', 'toy_repdup', 'toy_circ')) {
        fx = load_fixture(nm)
        full = streduce(fx$gg$nodes$gr)
        a = fs(fx$gw[[1]], target_region = full, pix.size = 1e4, gm.out = FALSE)
        tiles = dt2gr(a$gr)
        grid = gr.tile(full, 1e4) # blocks of whole grid tiles, so the sub-target tiles are full-target tiles
        sub = streduce(grid[((seq_along(grid) - 1) %/% 3) %% 2 == 0])
        b = fs(fx$gw[[1]], target_region = sub, pix.size = 1e4, gm.out = FALSE)
        bt = dt2gr(b$gr)
        m = match(paste(seqnames(bt), start(bt), end(bt), bt$node.id), paste(seqnames(tiles), start(tiles), end(tiles), tiles$node.id)) # toy_circ has overlapping nodes
        expect_false(anyNA(m), info = nm)
        expect_lt(length(bt), length(tiles))
        ad = data.table::as.data.table(a$dat)[, .(i, j, va = value)]
        bd = data.table::as.data.table(b$dat)[, .(i = pmin(m[i], m[j]), j = pmax(m[i], m[j]), vb = value)]
        cmp = merge(bd, ad, by = c('i', 'j'), all.x = TRUE)
        expect_equal(cmp$vb, cmp$va, tolerance = 1e-10, info = nm)
    }
})
