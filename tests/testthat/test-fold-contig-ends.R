## fold_contig_ends() (R/graphstats.R): loose ends on non-canonical contigs inside the footprint fold back onto themselves.

## chr1 node (CN 2) joined right -> left to a contig node (CN 5) at CN 2: the contig has loose ends of CN 3 (left) and
## CN 5 (right), the chr1 node one of CN 2 (left)
toy_contig_graph = function(){
    nodes = GenomicRanges::GRanges(c('1', 'GL000220.1'), IRanges::IRanges(c(1001, 1), width = 1000), cn = c(2, 5))
    edges = data.table(n1 = 1, n2 = 2, n1.side = 'right', n2.side = 'left', cn = 2, type = 'ALT')
    gg = loosefix(gGnome::gG(nodes = nodes, edges = edges))
    gg$set(y.field = 'cn')
    gg
}
loose_cn = function(gg){
    lo = gr2dt(gg$loose)[cn > 0]
    setNames(lo$cn, paste(lo$node.id, lo$orientation))
}

test_that('fold_contig_ends closes contig loose ends at floor(cn / 2) and leaves canonical ones', {
    gg = toy_contig_graph()
    before = gg$edges$dt$cn
    g2 = fold_contig_ends(gg, GenomicRanges::GRanges(c('1', 'GL000220.1'), IRanges::IRanges(1, 1e4)))
    fb = g2$edges$dt[type == 'ALT' & n1 == 2 & n2 == 2]
    expect_equal(nrow(fb), 2)
    expect_equal(sort(fb$cn), c(1, 2)) # left 3 -> 1 (one copy stays loose), right 5 -> 2
    expect_true(all(fb$n1.side == fb$n2.side))
    lc = loose_cn(g2)
    expect_equal(unname(lc[grepl('^2 ', names(lc))]), c(1, 1))
    expect_equal(sum(lc[grepl('^1 ', names(lc))]), 2) # canonical loose end untouched
    expect_identical(gg$edges$dt$cn, before) # input graph unchanged
})

test_that('fold_contig_ends ignores contig loose ends outside the footprint', {
    gg = toy_contig_graph()
    g2 = fold_contig_ends(gg, GenomicRanges::GRanges('1', IRanges::IRanges(1, 1e4)))
    expect_equal(nrow(g2$edges$dt), nrow(gg$edges$dt))
    expect_equal(loose_cn(g2), loose_cn(gg))
})
