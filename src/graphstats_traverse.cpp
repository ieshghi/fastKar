#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Exact, linear-time re-implementation of the reference traverse_graph_cpp
// (src/traverse_graph.cpp at fastKar 0a5497e). The output, including walk
// order, walk orientation and the rotation of cycles, is identical to the
// reference for any input; only the lookups change:
//   * "first row (in row order) whose left or right equals edge e and whose id
//     is not yet visited" uses a per-edge CSR list of rows instead of a scan
//     over all rows (or over the remaining rows in the cycle phase);
//   * visited rows are tracked per distinct id value (the reference keys its
//     std::set on id, not on the row), visited edges per distinct edge value;
//   * the next loose end to start from is found with a monotone pointer into
//     the sorted distinct loose ends instead of rebuilding a std::set;
//   * the "next edge is a loose end not yet visited when this path started"
//     test uses the identity  in_snapshot(e) = is_loose(e) && (!visited(e) || e == start),
//     which holds because the only loose end a path visits before it stops is
//     its start (any other loose end met as next_edge either stops the path or
//     was already visited before the path began).
// Edge values are compressed to dense codes over the union of left, right and
// (int-truncated, as std::set<int>::insert(double) does) loose ends, so the
// code works for arbitrary integer values.

namespace {

struct Dense {
  std::vector<int> vals;  // sorted distinct values
  int code(int v) const {
    return (int)(std::lower_bound(vals.begin(), vals.end(), v) - vals.begin());
  }
};

}  // namespace

List traverse_graph_core(const int* n, const int* left, const int* right, const int* id,
                         int nrow, const NumericVector& loose_ends) {
  const int nl = loose_ends.size();
  std::vector<int> loose_int(nl);
  for (int k = 0; k < nl; k++) loose_int[k] = (int)loose_ends[k];

  // ---- dense codes for edge values ----
  Dense ed;
  ed.vals.reserve(2 * (size_t)nrow + nl);
  for (int i = 0; i < nrow; i++) { ed.vals.push_back(left[i]); ed.vals.push_back(right[i]); }
  for (int k = 0; k < nl; k++) ed.vals.push_back(loose_int[k]);
  std::sort(ed.vals.begin(), ed.vals.end());
  ed.vals.erase(std::unique(ed.vals.begin(), ed.vals.end()), ed.vals.end());
  const int ne = (int)ed.vals.size();
  std::vector<int> L(nrow), R(nrow);
  bool identity = true;
  for (int i = 0; i < nrow; i++) {
    L[i] = ed.code(left[i]);
    R[i] = ed.code(right[i]);
  }

  // ---- dense codes for row ids (visited is keyed on id, as in the reference) ----
  for (int i = 0; i < nrow; i++) if (id[i] != i + 1) { identity = false; break; }
  std::vector<int> idc(nrow);
  int nid = nrow;
  if (identity) {
    for (int i = 0; i < nrow; i++) idc[i] = i;
  } else {
    Dense idd;
    idd.vals.assign(id, id + nrow);
    std::sort(idd.vals.begin(), idd.vals.end());
    idd.vals.erase(std::unique(idd.vals.begin(), idd.vals.end()), idd.vals.end());
    nid = (int)idd.vals.size();
    for (int i = 0; i < nrow; i++) idc[i] = idd.code(id[i]);
  }

  // ---- CSR: rows containing each edge code, ascending row order, no repeats ----
  std::vector<int> start(ne + 1, 0);
  for (int i = 0; i < nrow; i++) {
    start[L[i] + 1]++;
    if (R[i] != L[i]) start[R[i] + 1]++;
  }
  for (int e = 0; e < ne; e++) start[e + 1] += start[e];
  std::vector<int> rows(start[ne]);
  {
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int i = 0; i < nrow; i++) {
      rows[fill[L[i]]++] = i;
      if (R[i] != L[i]) rows[fill[R[i]]++] = i;
    }
  }

  std::vector<char> vis_id(nid, 0), vis_edge(ne, 0), is_loose(ne, 0);
  for (int k = 0; k < nl; k++) is_loose[ed.code(loose_int[k])] = 1;

  auto find_row = [&](int e) -> int {
    for (int p = start[e]; p < start[e + 1]; p++) {
      int r = rows[p];
      if (!vis_id[idc[r]]) return r;
    }
    return -1;
  };

  std::vector<std::vector<int>> paths, cycles;

  // ---- part 1: paths from loose ends, smallest unvisited loose end first ----
  std::vector<int> loose_codes;
  loose_codes.reserve(nl);
  for (int k = 0; k < nl; k++) loose_codes.push_back(ed.code(loose_int[k]));
  std::sort(loose_codes.begin(), loose_codes.end());
  loose_codes.erase(std::unique(loose_codes.begin(), loose_codes.end()), loose_codes.end());
  size_t lp = 0;
  while (true) {
    while (lp < loose_codes.size() && vis_edge[loose_codes[lp]]) lp++;
    if (lp >= loose_codes.size()) break;
    const int s = loose_codes[lp];
    std::vector<int> nodepath;
    int cur = s;
    while (true) {
      vis_edge[cur] = 1;
      int r = find_row(cur);
      if (r < 0) break;
      vis_id[idc[r]] = 1;
      int nxt;
      if (L[r] == cur) { nodepath.push_back(n[r]); nxt = R[r]; }
      else             { nodepath.push_back(-n[r]); nxt = L[r]; }
      if (is_loose[nxt] && (!vis_edge[nxt] || nxt == s)) break;
      int r2 = find_row(nxt);
      if (r2 < 0) { vis_edge[nxt] = 1; break; }
      cur = nxt;
    }
    if (!nodepath.empty()) paths.push_back(std::move(nodepath));
  }

  // ---- part 2: cycles over the remaining rows, first unvisited row first ----
  int rp = 0;
  while (true) {
    while (rp < nrow && vis_id[idc[rp]]) rp++;
    if (rp >= nrow) break;
    const int s = rp;
    vis_id[idc[s]] = 1;
    std::vector<int> nodepath;
    nodepath.push_back(n[s]);
    int cur = R[s];
    while (true) {
      int r = find_row(cur);
      if (r < 0) { cycles.push_back(std::move(nodepath)); break; }
      vis_id[idc[r]] = 1;
      if (L[r] == cur) { nodepath.push_back(n[r]); cur = R[r]; }
      else             { nodepath.push_back(-n[r]); cur = L[r]; }
    }
  }

  const int np = (int)paths.size(), nc = (int)cycles.size();
  List combined(np + nc);
  LogicalVector circular(np + nc);
  for (int i = 0; i < np; i++) { combined[i] = wrap(paths[i]); circular[i] = false; }
  for (int i = 0; i < nc; i++) { combined[np + i] = wrap(cycles[i]); circular[np + i] = true; }
  return List::create(_["snode.id"] = combined, _["circular"] = circular);
}

// Same traversal as traverse_graph_cpp, taking the wiring columns as vectors so
// callers that only permute `right` need not rebuild a data.frame per step.
// [[Rcpp::export]]
List traverse_graph_fast_cpp(IntegerVector n, IntegerVector left, IntegerVector right,
                             IntegerVector id, NumericVector loose_ends) {
  const int nrow = n.size();
  if (left.size() != nrow || right.size() != nrow || id.size() != nrow)
    stop("n, left, right and id must have the same length");
  return traverse_graph_core(n.begin(), left.begin(), right.begin(), id.begin(), nrow, loose_ends);
}
