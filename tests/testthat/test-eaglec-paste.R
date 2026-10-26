## Loose ends joined by Hi-C junction calls (read_eaglec / paste_loose_ends_eaglec, R/graphstats.R).

## a REF chain 1:1001-6000 of five 1 kb nodes at CN 2, 3, 2, 4, 2 (REF edges CN 2), so the CN 3 node has a loose end of
## CN 1 on each side and the CN 4 node one of CN 2 on each side: right of node 2 at 1:3000 ('-'), left of node 4 at 1:4001 ('+')
toy_loose_graph = function(){
    nodes = GenomicRanges::GRanges('1', IRanges::IRanges(seq(1001, 5001, 1000), width = 1000), cn = c(2, 3, 2, 4, 2))
    edges = data.table(n1 = 1:4, n2 = 2:5, n1.side = 'right', n2.side = 'left', cn = 2, type = 'REF')
    gg = loosefix(gGnome::gG(nodes = nodes, edges = edges))
    gg$set(y.field = 'cn')
    gg
}

## an EagleC combined table: chr1 pos1 chr2 pos2 then P(++), P(+-), P(-+), P(--)
write_eaglec = function(rows){
    f = tempfile(fileext = '.txt')
    data.table::fwrite(rows, f, sep = '\t', col.names = FALSE)
    f
}

test_that('read_eaglec flips EagleC orientation to gGnome strand and filters by probability', {
    f = write_eaglec(data.table(c('chr1', 'chr2'), c(2000, 10), c('chr1', 'chr2'), c(5001, 50), c(0, 0.05), c(0.9, 0), c(0.2, 0), c(0, 0)))
    calls = read_eaglec(f, min_p = 0.1)
    expect_equal(nrow(calls), 2)
    expect_equal(calls$p, c(0.9, 0.2)) #highest first, the chr2 call (0.05) dropped
    expect_equal(calls$chr1, c('1', '1'))
    expect_equal(calls[1, c(str1, str2)], c('-', '+')) #EagleC +- = right side of the first end, left side of the second
    expect_equal(calls[2, c(str1, str2)], c('+', '-'))
})

test_that('paste_loose_ends_eaglec joins only loose ends a call links, at the smaller residual CN', {
    gg = toy_loose_graph()
    ## node 2's right end (1:3000, '-') to node 4's left end (1:4001, '+'); a second call from node 2's left end to a non-loose position
    f = write_eaglec(data.table(c('chr1', 'chr1'), c(3000, 2001), c('chr1', 'chr1'), c(4001, 5500), c(0, 0), c(0.9, 0), c(0, 0), c(0, 0.8)))
    calls = read_eaglec(f, min_p = 0.1)
    pe = paste_loose_ends_eaglec(gg, calls, tol = 50, return_all = TRUE) #toy nodes are 1 kb, the default tol 10 kb
    alt = pe$gg$edges$dt[type == 'ALT']
    expect_equal(nrow(alt), 1)
    expect_equal(alt$cn, 1) #min(1, 2)
    expect_equal(sort(pe$calls$partner), c('loose end in region', 'no loose end'))
    lo = gr2dt(pe$gg$loose)[cn > 0 & !terminal] #not the chain ends
    ## left over: node 2 left 1, node 4 left 2 - 1 = 1, node 4 right 2; node 2 right used up
    expect_equal(sum(lo$cn), 4)
    expect_equal(nrow(gg$edges$dt), 4) #input graph untouched
})

test_that('paste_loose_ends_eaglec needs one end in region and leaves everything loose without calls', {
    gg = toy_loose_graph()
    f = write_eaglec(data.table('chr1', 3000, 'chr1', 4001, 0, 0.9, 0, 0))
    calls = read_eaglec(f)
    away = GenomicRanges::GRanges('1', IRanges::IRanges(8001, 9000))
    expect_equal(nrow(paste_loose_ends_eaglec(gg, calls, region = away, tol = 50)$edges$dt[type == 'ALT']), 0)
    nb = paste_loose_ends_eaglec(gg, calls, region = GenomicRanges::GRanges('1', IRanges::IRanges(4001, 4100)), tol = 50)
    expect_equal(nrow(nb$edges$dt[type == 'ALT']), 1) #node 4 is in region, node 2 outside: still joined
    none = paste_loose_ends_eaglec(gg, calls[0])
    expect_equal(sum(gr2dt(none$loose)[cn > 0 & !terminal]$cn), 6)
})
