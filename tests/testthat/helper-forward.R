## Helpers for the forward-simulator equivalence, property and downstream tests.

slow_tests = function() Sys.getenv('FASTKAR_SLOW_TESTS') == '1'

## structure of a long-format table: column names, key, column classes
expect_same_shape = function(new, old, info = NULL) {
    testthat::expect_identical(names(new), names(old), info = paste(info, 'names'))
    testthat::expect_identical(data.table::key(new), data.table::key(old), info = paste(info, 'key'))
    testthat::expect_identical(lapply(new, class), lapply(old, class), info = paste(info, 'classes'))
}

## any forward_simulate()/simulate_walks() output: gMatrix, list(gr, dat), data.table, or a list of those
expect_fwd_equiv = function(new, old, info = NULL, tol = EQUIV_TOL) {
    if (inherits(old, 'gMatrix')) {
        testthat::expect_true(inherits(new, 'gMatrix'), info = info)
        expect_gm_equiv(new, old, tol = tol, info = info)
        return(invisible())
    }
    if (data.table::is.data.table(old)) {
        expect_same_shape(new, old, info)
        expect_pix_equiv(new, old, tol = tol, info = info)
        ## other columns (widthprod of per-walk tables) in pixel order
        for (cc in setdiff(names(old), c('i', 'j', 'id', 'value')))
            testthat::expect_identical(new[[cc]][order(new$id)], old[[cc]][order(old$id)], info = paste(info, cc))
        return(invisible())
    }
    if (is.list(old) && !is.null(old$dat) && !is.null(old$gr)) {
        testthat::expect_identical(names(new), names(old), info = info)
        testthat::expect_identical(new$gr, old$gr, info = paste(info, 'gr'))
        expect_fwd_equiv(new$dat, old$dat, info = info, tol = tol)
        return(invisible())
    }
    testthat::expect_identical(length(new), length(old), info = paste(info, 'length'))
    for (k in seq_along(old)) expect_fwd_equiv(new[[k]], old[[k]], info = paste(info, 'walk', k), tol = tol)
}

## run a function from both packages on the same arguments; both must error or both succeed
run_both = function(name, args) {
    old = try(do.call(ref(name), args), silent = TRUE)
    new = try(do.call(opt(name), args), silent = TRUE)
    list(new = new, old = old)
}

expect_both_equiv = function(name, args, info = NULL, tol = EQUIV_TOL) {
    r = run_both(name, args)
    if (inherits(r$old, 'try-error')) {
        testthat::expect_true(inherits(r$new, 'try-error'), info = paste(info, 'reference errors:', r$old))
        return(invisible(r))
    }
    if (inherits(r$new, 'try-error')) testthat::fail(paste(info, 'new errors:', r$new))
    else expect_fwd_equiv(r$new, r$old, info = info, tol = tol)
    invisible(r)
}

## target regions for a fixture: footprint, strict sub-window, window past the graph, disjoint ranges
fixture_targets = function(fx) {
    ft = GenomicRanges::reduce(fx$ft[, c()])
    gg_gr = fx$gg$gr
    big = ft[which.max(GenomicRanges::width(ft))]
    w = GenomicRanges::width(big)
    sub = GenomicRanges::GRanges(GenomicRanges::seqnames(big),
                                 IRanges::IRanges(GenomicRanges::start(big) + floor(w / 4), GenomicRanges::start(big) + floor(3 * w / 4)))
    past = GenomicRanges::GRanges(GenomicRanges::seqnames(big),
                                  IRanges::IRanges(max(1, GenomicRanges::start(big) - 5e5), GenomicRanges::end(big) + 1e6))
    ## two disjoint pieces on the two chromosomes with the most graph sequence
    chr = names(sort(tapply(as.numeric(GenomicRanges::width(gg_gr)), as.character(GenomicRanges::seqnames(gg_gr)), sum), decreasing = TRUE))
    pieces = lapply(chr[seq_len(min(2, length(chr)))], function(cc) {
        g = gg_gr[as.character(GenomicRanges::seqnames(gg_gr)) == cc]
        s = min(GenomicRanges::start(g)); e = max(GenomicRanges::end(g))
        GenomicRanges::GRanges(cc, IRanges::IRanges(s + floor((e - s) / 3), s + floor((e - s) / 3) + max(floor((e - s) / 4), 1e5)))
    })
    multi = suppressWarnings(do.call(c, pieces))
    if (length(chr) < 2) {
        multi = c(GenomicRanges::GRanges(GenomicRanges::seqnames(big), IRanges::IRanges(GenomicRanges::start(big), GenomicRanges::start(big) + floor(w / 5))),
                  GenomicRanges::GRanges(GenomicRanges::seqnames(big), IRanges::IRanges(GenomicRanges::end(big) - floor(w / 5), GenomicRanges::end(big))))
    }
    list(footprint = fx$ft, sub = sub, past = past, multi = multi)
}
