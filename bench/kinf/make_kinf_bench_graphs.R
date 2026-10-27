# make_kinf_bench_graphs.R: cut graphs for the kinf convergence benchmark (bench/kinf/run_kinf_bench.R).
# Events: the 5 v3 CCLE events of the karyotype-inference benchmark (Fig. 4 local identifiability set, lines with
# Hi-C). Graph prep as the Fig. 4 local analysis / KI scripts/export_v3.R load_event_graph(): buffered footprint
# (0.2 x width, >= 10 kb) = target region; graph = loosefix(gg$copy$disjoin(ft) %&% ft); pix.size = ~150 tiles
# (multiple of 1 kb, >= 5 kb). Saved as node/edge tables (a few KB; a serialized gGraph is ~1 MB); rebuild with
# gG(nodes = x$nodes, edges = x$edges) (checked: identical nodes, edges and sample.gwalks() draws).
# Usage (from the repo root): Rscript bench/kinf/make_kinf_bench_graphs.R
suppressMessages(devtools::load_all('.', quiet = TRUE))
data.table::setDTthreads(1)
assets = '~/projects/karyotyping/data/figure4_local_identifiability_flow/ccle_v3_local_swap_event_assets.rds'
events = c('G32831.A2058_event_928', 'G32831.A2058_event_929', 'G32831.HCC38_event_1423',
           'G32831.HCC1937_event_1365', 'G32831.HCC38_event_1427')
ev = as.data.table(readRDS(assets))[event_id %in% events]
stopifnot(nrow(ev) == length(events))
for (r in seq_len(nrow(ev))){
  e = ev[r]
  gg = readRDS(normalizePath(path.expand(e$gg), mustWork = TRUE))
  ft0 = streduce(parse.gr(strsplit(e$footprint, ';', fixed = TRUE)[[1]]))
  buffer = max(sum(width(ft0)) * 0.2, 1e4)
  ft = streduce(ft0 + buffer)
  graph = loosefix(gg$copy$disjoin(ft) %&% ft)
  pix = max(5e3, round(sum(width(ft)) / 150 / 1e3) * 1e3)
  nd = graph$nodes$dt
  cat(e$event_id, e$type, 'nodes', nrow(nd), 'CN>1', sum(nd$cn > 1), 'max CN', max(nd$cn), 'pix', pix, '\n')
  saveRDS(list(nodes = graph$nodes$gr[, 'cn'], edges = graph$edges$dt[, .(n1, n2, n1.side, n2.side, cn, type)], ft = ft, pix.size = pix, event_id = e$event_id, type = e$type, pair = e$pair),
          file.path('bench/kinf/graphs', paste0(e$event_id, '.rds')))
}
