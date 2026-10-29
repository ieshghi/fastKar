## The buffered footprint must not depend on DEFAULT_GENOME / DEFAULT_BSGENOME (graph_seqlengths()).
## Without either, parse.gr() ends each chromosome at the footprint's own end, and the buffer used to be cut on
## the right (found on HMF, run on a machine without the env set).

no_default_genome = function(code) {
    old = Sys.getenv(c('DEFAULT_GENOME', 'DEFAULT_BSGENOME'), unset = NA)
    Sys.unsetenv(c('DEFAULT_GENOME', 'DEFAULT_BSGENOME'))
    on.exit(for (v in names(old)) if (!is.na(old[[v]])) do.call(Sys.setenv, as.list(old[v])))
    suppressWarnings(code)
}

test_that('identifiability_context buffers the footprint on both sides without a default genome', {
    fx = load_fixture('ccle_pyrgo')
    fp = gUtils::gr.string(gUtils::streduce(fx$ft))
    ctx = no_default_genome(opt('identifiability_context')(fx$gg, footprint = gUtils::parse.gr(fp)))
    ft0 = ctx$footprint0
    sl = GenomeInfoDb::seqlengths(fx$gg$nodes$gr)[as.character(GenomicRanges::seqnames(ft0))]
    expect_equal(min(GenomicRanges::start(ctx$target_region)), max(1, min(GenomicRanges::start(ft0)) - ctx$footprint_buffer))
    expect_equal(max(GenomicRanges::end(ctx$target_region)), min(sl, max(GenomicRanges::end(ft0)) + ctx$footprint_buffer))
    expect_gt(max(GenomicRanges::end(ctx$target_region)), max(GenomicRanges::end(ft0)))
})

test_that('graph_seqlengths gives a parsed range the graph lengths', {
    fx = load_fixture('ccle_pyrgo')
    gr = no_default_genome(gUtils::parse.gr(gUtils::gr.string(fx$ft[1])))
    chr = as.character(GenomicRanges::seqnames(gr))
    expect_equal(GenomeInfoDb::seqlengths(opt('graph_seqlengths')(gr, fx$gg))[[chr]],
                 GenomeInfoDb::seqlengths(fx$gg$nodes$gr)[[chr]])
})
