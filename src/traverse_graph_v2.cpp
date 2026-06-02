#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// O(V) per-sample replacement for traverse_graph_cpp.
// Same return shape: list("snode.id" = list_of_int_vectors, "circular" = logical).
// Kept alongside the original traverse_graph_cpp so we can A/B benchmark.
//
// Key changes vs traverse_graph_cpp:
//   * Per-step row lookup is O(1) via precomputed subid -> (row, side) tables
//     instead of a linear scan over all rows.
//   * Visited-row tracking is a flat std::vector<char> instead of std::set.
//   * Each linear walk is opened exactly once: the terminating loose end is
//     marked consumed when reached, so we don't redundantly re-enter from it.
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
List traverse_graph_v2_cpp(DataFrame A, NumericVector loose_ends) {
  IntegerVector n_col     = A["n"];
  IntegerVector left_col  = A["left"];
  IntegerVector right_col = A["right"];
  const int nrow = n_col.size();

  // ---------------- 1. edge_subid -> (row, side) lookup ----------------
  // Each edge subid attaches to at most 2 rows of the wiring:
  //   non-loose edges: 2 (one as left of some node-copy, one as right of another)
  //   loose edges:     1
  // side encoding: 0 = left of node-copy, 1 = right.

  int max_subid = 0;
  for (int i = 0; i < nrow; i++) {
    if (left_col[i]  > max_subid) max_subid = left_col[i];
    if (right_col[i] > max_subid) max_subid = right_col[i];
  }
  // Also accommodate any subids that appear only in loose_ends (defensive).
  for (int k = 0; k < loose_ends.size(); k++) {
    int s = (int)loose_ends[k];
    if (s > max_subid) max_subid = s;
  }
  const int N = max_subid + 1;

  std::vector<int>  row_a(N, -1), row_b(N, -1);
  std::vector<char> side_a(N, 0), side_b(N, 0);

  for (int i = 0; i < nrow; i++) {
    int s = left_col[i];
    if (s >= 0 && s < N) {
      if (row_a[s] == -1) { row_a[s] = i; side_a[s] = 0; }
      else                { row_b[s] = i; side_b[s] = 0; }
    }
    s = right_col[i];
    if (s >= 0 && s < N) {
      if (row_a[s] == -1) { row_a[s] = i; side_a[s] = 1; }
      else                { row_b[s] = i; side_b[s] = 1; }
    }
  }

  std::vector<char> is_loose(N, 0);
  for (int k = 0; k < loose_ends.size(); k++) {
    int s = (int)loose_ends[k];
    if (s >= 0 && s < N) is_loose[s] = 1;
  }

  auto other_row = [&](int subid, int from_row) -> int {
    if (subid < 0 || subid >= N) return -1;
    if (row_a[subid] != -1 && row_a[subid] != from_row) return row_a[subid];
    if (row_b[subid] != -1 && row_b[subid] != from_row) return row_b[subid];
    return -1;
  };
  auto side_at = [&](int subid, int row) -> int {
    return (row_a[subid] == row) ? side_a[subid] : side_b[subid];
  };

  // ---------------- 2. Traversal ----------------
  std::vector<char> visited(nrow, 0);
  std::vector<char> loose_consumed(N, 0);
  List paths, cycles;

  // 2a. Linear walks from loose ends.
  for (int k = 0; k < loose_ends.size(); k++) {
    int start = (int)loose_ends[k];
    if (start < 0 || start >= N)   continue;
    if (loose_consumed[start])     continue;
    if (row_a[start] == -1)        continue;

    std::vector<int> nodepath;
    nodepath.reserve(nrow);

    int current_edge = start;
    int from_row     = -1;

    while (true) {
      int row = other_row(current_edge, from_row);
      if (row == -1 || visited[row]) break;
      visited[row] = 1;

      int side = side_at(current_edge, row);
      int next_edge;
      if (side == 0) {                          // entered left, exit right
        nodepath.push_back( n_col[row]);
        next_edge = right_col[row];
      } else {                                  // entered right, exit left
        nodepath.push_back(-n_col[row]);
        next_edge = left_col[row];
      }

      from_row     = row;
      current_edge = next_edge;

      if (current_edge >= 0 && current_edge < N && is_loose[current_edge]) {
        loose_consumed[current_edge] = 1;
        break;
      }
    }

    loose_consumed[start] = 1;
    if (!nodepath.empty()) {
      paths.push_back(IntegerVector(nodepath.begin(), nodepath.end()));
    }
  }

  // 2b. Cycles: any remaining unvisited rows form cycles.
  for (int start_row = 0; start_row < nrow; start_row++) {
    if (visited[start_row]) continue;

    std::vector<int> nodepath;
    nodepath.reserve(nrow);
    nodepath.push_back(n_col[start_row]);
    visited[start_row] = 1;

    int current_edge = right_col[start_row];
    int from_row     = start_row;

    while (true) {
      int row = other_row(current_edge, from_row);
      if (row == -1 || visited[row]) break;
      visited[row] = 1;

      int side = side_at(current_edge, row);
      int next_edge;
      if (side == 0) {
        nodepath.push_back( n_col[row]);
        next_edge = right_col[row];
      } else {
        nodepath.push_back(-n_col[row]);
        next_edge = left_col[row];
      }

      from_row     = row;
      current_edge = next_edge;
    }

    cycles.push_back(IntegerVector(nodepath.begin(), nodepath.end()));
  }

  // ---------------- 3. Combine into baseline return shape ----------------
  const int n_paths  = paths.size();
  const int n_cycles = cycles.size();
  const int n_total  = n_paths + n_cycles;

  List combined(n_total);
  LogicalVector circular(n_total);

  for (int i = 0; i < n_paths;  i++) { combined[i]           = paths[i];  circular[i]           = false; }
  for (int i = 0; i < n_cycles; i++) { combined[n_paths + i] = cycles[i]; circular[n_paths + i] = true;  }

  return List::create(
    _["snode.id"] = combined,
    _["circular"] = circular
  );
}
