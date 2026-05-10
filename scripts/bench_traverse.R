# scripts/bench_traverse.R
#
# Compare traverse_graph_cpp (baseline) vs traverse_graph_v2_cpp on a battery
# of synthetic graphs (and any *.rds graphs in bench/graphs/), verify output
# equivalence, append to bench/history.rds.
#
# Run from the fastKar package root:
#     Rscript scripts/bench_traverse.R
#
# To include real graphs, drop gGraph .rds files into bench/graphs/ — the
# harness will pick them up automatically.

suppressPackageStartupMessages({
  library(Rcpp)         # forward.R uses cppFunction at top level
  library(devtools)
  library(data.table)
  library(gGnome)       # must provide loosefix(); use a fork that exports it
  library(skitools)
})

devtools::load_all(".")

results <- bench_all(
  Ns    = c(1000, 10000),
  impls = list(
    v1 = traverse_graph_cpp,
    v2 = traverse_graph_v2_cpp
  )
)

cat("\n=== This run (seconds) ===\n")
print(bench_summary(results))

cat("\n=== Speedup (v1 / v2) per phase ===\n")
spd <- dcast(results, graph + N ~ version + phase, value.var = "seconds")
spd[, traverse_speedup := v1_traverse / v2_traverse]
spd[, hash_speedup     := v1_hash     / v2_hash]
print(spd[, .(graph, N,
              v1_traverse, v2_traverse, traverse_speedup,
              v1_hash,     v2_hash,     hash_speedup)])

cat("\nHistory appended to bench/history.rds\n")
