## obs_words in longread_probdist() and readdist_probdist() (read-length mixture)

## longread_probdist / .probdist_from_reads as on master before obs_words existed,
legacy_probdist = function(gwlist, readL, minsize = 0, background = 1e-10, use.nodes = NULL) {
    if (is(gwlist, 'gWalk')) gwlist = list(gwlist)
    reads_list = lapply(gwlist, function(x) opt('reads_fromwalk')(x, readL, minsize = minsize, use.nodes = use.nodes))
    all_words = unique(unname(unlist(lapply(reads_list, function(x) x$words))))
    lapply(reads_list, function(x) {
        cx = setNames(unname(x$nums)[match(all_words, x$words)], all_words)
        cx[is.na(cx)] = 0
        (cx + background) / (sum(cx) + background * length(all_words))
    })
}

rd_fixture = function() {
    fx = load_fixture('toy_dtrans')
    list(gw = fx$gw[1:2], fx = fx)
}

test_that('obs_words = NULL leaves longread_probdist unchanged', {
    g = rd_fixture()$gw
    for (L in c(50, 1e3, 1e5)) {
        new = opt('longread_probdist')(g, L)
        old = legacy_probdist(g, L)
        expect_equal(new, old, tolerance = 1e-12)
        expect_identical(names(new[[1]]), names(new[[2]]))
        expect_identical(new, opt('longread_probdist')(g, L, obs_words = NULL))
    }
})

test_that('obs_words adds unseen words with background probability and keeps distributions normalized', {
    g = rd_fixture()$gw
    L = 1e3
    base = opt('longread_probdist')(g, L)
    extra = c('zz|yy|xx', 'ww|vv')
    bg = 1e-10
    out = opt('longread_probdist')(g, L, background = bg, obs_words = c(extra, names(base[[1]])[1]))
    for (p in out) {
        expect_equal(sum(p), 1, tolerance = 1e-12)
        expect_true(all(extra %in% names(p)))
        expect_false(anyDuplicated(names(p)) > 0)
    }
    ## unseen words: background / (total count + background * vocabulary size)
    cnt = sum(opt('reads_fromwalk')(g[[1]], L)$nums)
    expect_equal(unname(out[[1]][extra[1]]), bg / (cnt + bg * length(out[[1]])), tolerance = 1e-12)
    ## words already in the vocabulary are not duplicated; vocabulary grows by exactly the new words
    expect_equal(length(out[[1]]), length(base[[1]]) + length(extra))
    ## adding "" / NA observed words must not give them counts
    out2 = opt('longread_probdist')(g, L, obs_words = c('', NA))
    for (p in out2) expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that('readdist_probdist with a single read length equals longread_probdist', {
    g = rd_fixture()$gw
    for (L in c(100, 5e3)) {
        a = opt('readdist_probdist')(g, rep(L, 25))
        b = opt('longread_probdist')(g, L)
        expect_length(a, length(g))
        for (k in seq_along(g)) {
            expect_equal(sort(names(a[[k]])), sort(names(b[[k]])))
            expect_equal(a[[k]][names(b[[k]])], b[[k]], tolerance = 1e-12)
        }
    }
})

test_that('readdist_probdist is a normalized mixture of per-bin distributions', {
    g = rd_fixture()$gw
    set.seed(1)
    rl = round(exp(runif(200, log(50), log(1e4))))
    extra = 'zz|yy'
    out = opt('readdist_probdist')(g, rl, bins = 4, obs_words = extra)
    expect_length(out, length(g))
    for (p in out) {
        expect_equal(sum(p), 1, tolerance = 1e-9)
        expect_true(extra %in% names(p))
        expect_false(anyNA(p))
        expect_true(all(p > 0))
    }
    ## use.nodes is passed through
    nd = rd_fixture()$fx$gg$nodes$dt$node.id[-1]
    out2 = opt('readdist_probdist')(g, rl, bins = 4, use.nodes = nd)
    for (p in out2) expect_equal(sum(p), 1, tolerance = 1e-9)
    ## a single gWalk is accepted
    expect_length(opt('readdist_probdist')(g[[1]], rl, bins = 3), 1)
})

test_that('readdist_probdist does not break when cut() creates empty bins', {
    g = rd_fixture()$gw
    ## two well separated clusters of read lengths: the middle bins are empty
    rl = c(rep(100, 30), rep(1e5, 10))
    out = opt('readdist_probdist')(g, rl, bins = 8)
    for (p in out) {
        expect_equal(sum(p), 1, tolerance = 1e-9)
        expect_false(anyNA(p))
    }
    ## equals the explicit two-component mixture
    pa = opt('longread_probdist')(g, 100)
    pb = opt('longread_probdist')(g, 1e5)
    for (k in seq_along(g)) {
        w = unique(c(names(pa[[k]]), names(pb[[k]])))
        fa = pa[[k]][w]; fa[is.na(fa)] = 0
        fb = pb[[k]][w]; fb[is.na(fb)] = 0
        expect_equal(unname(out[[k]][w]), unname(0.75 * fa + 0.25 * fb), tolerance = 1e-12)
    }
})
