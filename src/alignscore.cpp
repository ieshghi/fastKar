#include <Rcpp.h>
#include <vector>
#include <map>
#include <cstdlib>
using namespace Rcpp;

// -----------------------------
// Core alignment (fast DP)
// -----------------------------
double alignscore_core(const IntegerVector& x,
                       const IntegerVector& y,
                       const NumericVector& penalty) {

  int n = x.size();
  int m = y.size();

  NumericVector gx(n), gy(m);

  for (int i = 0; i < n; i++) {
    gx[i] = penalty[ std::abs(x[i])-1 ];
  }
  for (int j = 0; j < m; j++) {
    gy[j] = penalty[ std::abs(y[j])-1 ];
  }

  NumericVector prev(m + 1), curr(m + 1);

  for (int j = 1; j <= m; j++) {
    prev[j] = prev[j - 1] + gy[j - 1];
  }

  for (int i = 1; i <= n; i++) {

    curr[0] = prev[0] + gx[i - 1];
    int xi = x[i - 1];
    double gxi = gx[i - 1];

    for (int j = 1; j <= m; j++) {

      double gyj = gy[j - 1];

      double match = prev[j - 1] +
        (xi == y[j - 1] ? 0.0 : (gxi + gyj));

      double del = prev[j] + gxi;
      double ins = curr[j - 1] + gyj;

      double best = match > del ? match : del;
      curr[j] = best > ins ? best : ins;
    }

    std::swap(prev, curr);
  }

  return prev[m];
}
// -----------------------------
// reverse-complement score
// -----------------------------
double alignscore_compl_core(const IntegerVector& x,
                             const IntegerVector& y,
                             const NumericVector& penalty) {

  double s1 = alignscore_core(x, y, penalty);

  if (y.size() == 0) return s1;

  int m = y.size();
  IntegerVector y2(m);

  for (int i = 0; i < m; i++) {
    y2[i] = -y[m - 1 - i];
  }

  double s2 = alignscore_core(x, y2, penalty);

  return s1 > s2 ? s1 : s2;
}

// -----------------------------
// COST MATRIX
// -----------------------------
// Same scores as calling alignscore_compl_core() for every (i, j) (the reference loop, kept above),
// computed without per-pair R allocations: each walk's penalties and each y's reverse complement
// are built once, and the DP rows are reused. The DP performs the same floating-point operations
// in the same order, so every entry is bit-identical.
namespace {
struct Seq {
  std::vector<int> v, rc;
  std::vector<double> g, grc;
};
inline double dp_score(const Seq& X, const std::vector<int>& y, const std::vector<double>& gy,
                       std::vector<double>& prev, std::vector<double>& curr) {
  const int n = (int)X.v.size(), m = (int)y.size();
  prev.assign(m + 1, 0.0);
  curr.assign(m + 1, 0.0);
  for (int j = 1; j <= m; j++) prev[j] = prev[j - 1] + gy[j - 1];
  for (int i = 1; i <= n; i++) {
    curr[0] = prev[0] + X.g[i - 1];
    const int xi = X.v[i - 1];
    const double gxi = X.g[i - 1];
    for (int j = 1; j <= m; j++) {
      const double gyj = gy[j - 1];
      const double match = prev[j - 1] + (xi == y[j - 1] ? 0.0 : (gxi + gyj));
      const double del = prev[j] + gxi;
      const double ins = curr[j - 1] + gyj;
      const double best = match > del ? match : del;
      curr[j] = best > ins ? best : ins;
    }
    std::swap(prev, curr);
  }
  return prev[m];
}
Seq make_seq(SEXP s, const NumericVector& penalty) {
  Seq out;
  IntegerVector v(s);
  const int n = v.size();
  out.v.assign(v.begin(), v.end());
  out.g.resize(n);
  out.rc.resize(n);
  out.grc.resize(n);
  for (int i = 0; i < n; i++) out.g[i] = penalty[std::abs(v[i]) - 1];
  for (int i = 0; i < n; i++) {
    out.rc[i] = -v[n - 1 - i];
    out.grc[i] = penalty[std::abs(out.rc[i]) - 1];
  }
  return out;
}
}  // namespace

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
NumericMatrix compute_cost_matrix_cpp(List sn_x,
                                      List sn_y,
                                      NumericVector penalty) {
  int n = sn_x.size();
  int m = sn_y.size();
  int dim = std::max(n, m);
  std::vector<Seq> X(dim), Y(dim);
  for (int i = 0; i < n; i++) X[i] = make_seq(sn_x[i], penalty);
  for (int j = 0; j < m; j++) Y[j] = make_seq(sn_y[j], penalty);
  // identical walks give identical rows / columns: score each distinct (x, y) content pair once
  auto uniq = [](const std::vector<Seq>& V, std::vector<int>& code) {
    std::map<std::vector<int>, int> seen;
    std::vector<int> rep;
    code.resize(V.size());
    for (size_t i = 0; i < V.size(); i++) {
      auto it = seen.find(V[i].v);
      if (it == seen.end()) { seen.emplace(V[i].v, (int)rep.size()); code[i] = (int)rep.size(); rep.push_back((int)i); }
      else code[i] = it->second;
    }
    return rep;
  };
  std::vector<int> cx, cy;
  std::vector<int> rx = uniq(X, cx), ry = uniq(Y, cy);
  const int ux = (int)rx.size(), uy = (int)ry.size();
  std::vector<double> U((size_t)ux * uy);
  std::vector<double> prev, curr;
  for (int a = 0; a < ux; a++) {
    const Seq& xa = X[rx[a]];
    for (int b = 0; b < uy; b++) {
      const Seq& yb = Y[ry[b]];
      double s1 = dp_score(xa, yb.v, yb.g, prev, curr);
      double score = s1;
      if (!yb.v.empty()) {
        double s2 = dp_score(xa, yb.rc, yb.grc, prev, curr);
        score = s1 > s2 ? s1 : s2;
      }
      U[(size_t)a * uy + b] = -score;
    }
  }
  NumericMatrix costmat(dim, dim);
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++) costmat(i, j) = U[(size_t)cx[i] * uy + cy[j]];
  return costmat;
}
