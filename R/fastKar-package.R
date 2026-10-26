#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table data.table
## usethis namespace: end
NULL

## Imports, so fastKar works from an installed library() and not only under
## devtools::load_all() with its dependencies attached. BiocGenerics is
## imported whole: its S4 generics (start, end, width, strand, unique, unlist,
## ...) mask base functions, as they do when GenomicRanges is attached. %&% is
## gGnome's (the gUtils generic with the gGraph/gWalk/gNode methods), so
## gGraph %&% GRanges dispatches as when gGnome is attached last.
#' @import BiocGenerics
#' @importFrom S4Vectors mcols
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges GRangesList seqnames reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom data.table CJ as.data.table copy dcast first foverlaps fread is.data.table
#' @importFrom data.table merge.data.table rbindlist set setDT setattr setkey setkeyv setnames setorder
#' @importFrom gUtils %O% %Q% dt2gr gr.chr gr.end gr.findoverlaps gr.merge gr.nochr gr.start
#' @importFrom gUtils gr.stripstrand gr.tile gr2dt hg_seqlengths parse.gr si2gr streduce
#' @importFrom gGnome %&% balance embedloops gG gW
#' @importFrom GxG cooler gM straw
#' @importFrom Matrix sparseMatrix
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom pracma ceil mod
#' @importFrom skitools dodo.call ppdf qq_pval
#' @importFrom ggplot2 aes coord_equal element_blank element_text facet_wrap geom_hline geom_point
#' @importFrom ggplot2 geom_tile geom_vline ggplot labs theme theme_classic
#' @importFrom ggforce geom_sina
#' @importFrom viridis scale_fill_viridis
#' @rawNamespace exportPattern("^[[:alpha:]]")
NULL

## gGnome's loosefix() is internal (not exported).
loosefix <- function(...) utils::getFromNamespace('loosefix', 'gGnome')(...)
