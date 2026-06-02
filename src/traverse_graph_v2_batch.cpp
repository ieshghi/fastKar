#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Batched version of traverse_graph_v2_cpp: takes a List of N permutation
// vectors and returns a List of N walk-set Lists. Amortizes:
//   * the R<->C boundary cost (1 call instead of N)
//   * the allocation of scratch buffers (visited, loose_consumed, ...)
//
// Per-perm cost is the same algorithm as traverse_graph_v2_cpp:
//   rebuild attachment lookup (O(V)), reset visited bitmap, traverse.
//
// The lookup uses the same 2-slot model as traverse_graph_v2_cpp (row_a/row_b
// hold two rows per subid regardless of side) — necessary because inversion-
// type edges (where both endpoints are on the same side of their respective
// node-copies) cause the same subid to appear twice in left_col or right_col.
//
// API:
//   A          : DataFrame with columns n (int), left (int)
//                The 'right' column is ignored if present; perms supply it.
//   perms      : List of length N, each entry an IntegerVector of length nrow,
//                holding the per-sample 'right' assignment.
//   loose_ends : NumericVector of edge subids that are loose ends.
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
List traverse_graph_v2_batch_cpp(DataFrame A,
                                  List perms,
                                  NumericVector loose_ends) {
  IntegerVector n_col    = A["n"];
  IntegerVector left_col = A["left"];
  const int nrow = n_col.size();
  const int N    = perms.size();

  if (N == 0) return List();

  // ----- Determine max_subid across all inputs -----
  int max_subid = 0;
  for (int i = 0; i < nrow; i++) {
    if (left_col[i] > max_subid) max_subid = left_col[i];
  }
  for (int s = 0; s < N; s++) {
    IntegerVector perm = perms[s];
    if (perm.size() != nrow) {
      Rcpp::stop("perms[%d] has length %d; expected %d to match A",
                 s + 1, (int)perm.size(), nrow);
    }
    for (int i = 0; i < nrow; i++) {
      if (perm[i] > max_subid) max_subid = perm[i];
    }
  }
  for (int k = 0; k < loose_ends.size(); k++) {
    int s = (int)loose_ends[k];
    if (s > max_subid) max_subid = s;
  }
  const int N_subid = max_subid + 1;

  // ----- Loose-end bitmap (invariant) -----
  std::vector<char> is_loose(N_subid, 0);
  for (int k = 0; k < loose_ends.size(); k++) {
    int s = (int)loose_ends[k];
    if (s >= 0 && s < N_subid) is_loose[s] = 1;
  }

  // ----- Scratch buffers (reused per perm) -----
  // Two-slot lookup matching traverse_graph_v2_cpp's data model. Each subid
  // can attach to at most 2 rows; side: 0 = left of node-copy, 1 = right.
  std::vector<int>  row_a(N_subid),  row_b(N_subid);
  std::vector<char> side_a(N_subid), side_b(N_subid);
  std::vector<char> visited(nrow);
  std::vector<char> loose_consumed(N_subid);
  std::vector<int>  nodepath_buf;
  nodepath_buf.reserve(nrow);

  List result(N);

  // ----- Per-perm loop -----
  for (int s = 0; s < N; s++) {
    IntegerVector perm = perms[s];

    // Rebuild attachments: both left and right contribute. Reset row_a to -1
    // (row_b is set lazily once row_a fills, so we don't need to reset it
    // explicitly — but we DO need to reset row_a because it gates the if/else).
    std::fill(row_a.begin(), row_a.end(), -1);
    std::fill(row_b.begin(), row_b.end(), -1);
    for (int i = 0; i < nrow; i++) {
      int sub = left_col[i];
      if (sub >= 0 && sub < N_subid) {
        if (row_a[sub] == -1) { row_a[sub] = i; side_a[sub] = 0; }
        else                  { row_b[sub] = i; side_b[sub] = 0; }
      }
      sub = perm[i];
      if (sub >= 0 && sub < N_subid) {
        if (row_a[sub] == -1) { row_a[sub] = i; side_a[sub] = 1; }
        else                  { row_b[sub] = i; side_b[sub] = 1; }
      }
    }

    // Reset traversal state.
    std::fill(visited.begin(), visited.end(), 0);
    std::fill(loose_consumed.begin(), loose_consumed.end(), 0);

    auto other_row = [&](int subid, int from_row) -> int {
      if (subid < 0 || subid >= N_subid) return -1;
      if (row_a[subid] != -1 && row_a[subid] != from_row) return row_a[subid];
      if (row_b[subid] != -1 && row_b[subid] != from_row) return row_b[subid];
      return -1;
    };
    auto side_at = [&](int subid, int row) -> int {
      return (row_a[subid] == row) ? side_a[subid] : side_b[subid];
    };

    List paths, cycles;

    // 2a. Linear walks from loose ends.
    for (int k = 0; k < loose_ends.size(); k++) {
      int start = (int)loose_ends[k];
      if (start < 0 || start >= N_subid) continue;
      if (loose_consumed[start])         continue;
      if (row_a[start] == -1)            continue;

      nodepath_buf.clear();
      int current_edge = start;
      int from_row     = -1;

      while (true) {
        int row = other_row(current_edge, from_row);
        if (row == -1 || visited[row]) break;
        visited[row] = 1;

        int side = side_at(current_edge, row);
        int next_edge;
        if (side == 0) {                          // entered left, exit right
          nodepath_buf.push_back( n_col[row]);
          next_edge = perm[row];
        } else {                                  // entered right, exit left
          nodepath_buf.push_back(-n_col[row]);
          next_edge = left_col[row];
        }

        from_row     = row;
        current_edge = next_edge;

        if (current_edge >= 0 && current_edge < N_subid && is_loose[current_edge]) {
          loose_consumed[current_edge] = 1;
          break;
        }
      }

      loose_consumed[start] = 1;
      if (!nodepath_buf.empty()) {
        paths.push_back(IntegerVector(nodepath_buf.begin(), nodepath_buf.end()));
      }
    }

    // 2b. Cycles.
    for (int start_row = 0; start_row < nrow; start_row++) {
      if (visited[start_row]) continue;

      nodepath_buf.clear();
      nodepath_buf.push_back(n_col[start_row]);
      visited[start_row] = 1;

      int current_edge = perm[start_row];
      int from_row     = start_row;

      while (true) {
        int row = other_row(current_edge, from_row);
        if (row == -1 || visited[row]) break;
        visited[row] = 1;

        int side = side_at(current_edge, row);
        int next_edge;
        if (side == 0) {
          nodepath_buf.push_back( n_col[row]);
          next_edge = perm[row];
        } else {
          nodepath_buf.push_back(-n_col[row]);
          next_edge = left_col[row];
        }

        from_row     = row;
        current_edge = next_edge;
      }

      cycles.push_back(IntegerVector(nodepath_buf.begin(), nodepath_buf.end()));
    }

    // Combine into baseline return shape.
    const int n_paths  = paths.size();
    const int n_cycles = cycles.size();
    const int n_total  = n_paths + n_cycles;

    List combined(n_total);
    LogicalVector circular(n_total);

    for (int i = 0; i < n_paths;  i++) { combined[i]           = paths[i];  circular[i]           = false; }
    for (int i = 0; i < n_cycles; i++) { combined[n_paths + i] = cycles[i]; circular[n_paths + i] = true;  }

    result[s] = List::create(
      _["snode.id"] = combined,
      _["circular"] = circular
    );
  }

  return result;
}
