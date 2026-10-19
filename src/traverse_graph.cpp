#include <Rcpp.h>
using namespace Rcpp;

List traverse_graph_core(const int* n, const int* left, const int* right, const int* id,
                         int nrow, const NumericVector& loose_ends);

// given a wiring, with internal edges and loose ends, traverse graph and return paths and cycles.
// The original O(nrow^2) scan (fastKar 0a5497e) was replaced by the exact linear-time
// traverse_graph_core in graphstats_traverse.cpp; output is identical (tested against the oracle
// in tests/testthat/test-equiv-graphstats-sampling.R).
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
List traverse_graph_cpp(DataFrame A, NumericVector loose_ends) {
  IntegerVector n = A["n"];
  IntegerVector left = A["left"];
  IntegerVector right = A["right"];
  IntegerVector id = A["id"];
  return traverse_graph_core(n.begin(), left.begin(), right.begin(), id.begin(), n.size(), loose_ends);
}
