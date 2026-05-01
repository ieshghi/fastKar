#include <Rcpp.h>
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

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
NumericMatrix compute_cost_matrix_cpp(List sn_x,
                                      List sn_y,
                                      NumericVector penalty) {
  int n = sn_x.size();
  int m = sn_y.size();
  int dim = std::max(n, m);
  NumericMatrix costmat(dim, dim);
  for (int i = 0; i < dim; i++) {
    IntegerVector xi = (i < n) ? sn_x[i] : IntegerVector(0);
    for (int j = 0; j < dim; j++) {
      IntegerVector yj = (j < m) ? sn_y[j] : IntegerVector(0);
      double score = alignscore_compl_core(xi, yj, penalty);
      costmat(i, j) = -score;
    }
  }
  return costmat;
}
