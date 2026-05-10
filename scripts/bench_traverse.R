# scripts/bench_traverse.R
#
# Compare three configurations on a battery of synthetic graphs (and any
# .rds graphs in bench/graphs/), verify equivalence by partition, append to
# bench/history.rds.
#
#   v1     = traverse_graph_cpp     + hash_snodelist           (baseline)
#   v2t    = traverse_graph_v2_cpp  + hash_snodelist           (traversal-only)
#   v2t_h2 = traverse_graph_v2_cpp  + hash_karyotype_cpp       (full)
#
# Run from the fastKar package root:
#     Rscript scripts/bench_traverse.R
#
# Drop gGraph .rds files into bench/graphs/ to include real graphs.

suppressPackageStartupMessages({
  library(Rcpp)         # forward.R uses cppFunction at top level
  library(devtools)
  library(data.table)
  library(gGnome)       # must provide loosefix(); use a fork that exports it
  library(skitools)
})

devtools::load_all(".")

results <- bench_all(Ns = c(50, 500))

cat("\n=== This run (seconds) ===\n")
print(bench_summary(results))

cat("\n=== Speedups (vs v1 baseline) ===\n")
wide <- bench_summary(results)
wide[, traverse_speedup_v2t    := v1_traverse / v2t_traverse]
wide[, traverse_speedup_v2t_h2 := v1_traverse / v2t_h2_traverse]
wide[, hash_speedup_v2t        := v1_hash     / v2t_hash]
wide[, hash_speedup_v2t_h2     := v1_hash     / v2t_h2_hash]
wide[, total_speedup_v2t       := v1_total    / v2t_total]
wide[, total_speedup_v2t_h2    := v1_total    / v2t_h2_total]
print(wide[, .(graph, N,
               traverse_speedup_v2t, traverse_speedup_v2t_h2,
               hash_speedup_v2t,     hash_speedup_v2t_h2,
               total_speedup_v2t,    total_speedup_v2t_h2)])

cat("\nHistory appended to bench/history.rds\n")
