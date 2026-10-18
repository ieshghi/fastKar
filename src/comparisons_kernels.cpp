// Small kernels for R/comparisons.R: template (upper-triangle pixel table) layout checks,
// construction and indexing, without materialising make_template_dat() tables.
//
// Template layout for n tiles (as make_template_dat_cpp): rows (i, j) for i = 1..n, j = i..n,
// in that order; id = row number; widthprod = width[i] * width[j].
#include <Rcpp.h>
using namespace Rcpp;

// TRUE iff i, j, id are plain integer vectors (no attributes) that are exactly the template
// columns for n tiles, i.e. identical() to make_template_dat(gr)$i / $j / $id.
// [[Rcpp::export]]
bool cmp_template_match(SEXP i, SEXP j, SEXP id, int n) {
    if (TYPEOF(i) != INTSXP || TYPEOF(j) != INTSXP || TYPEOF(id) != INTSXP) return false;
    if (ATTRIB(i) != R_NilValue || ATTRIB(j) != R_NilValue || ATTRIB(id) != R_NilValue) return false;
    R_xlen_t P = (R_xlen_t)n * (n + 1) / 2;
    if (XLENGTH(i) != P || XLENGTH(j) != P || XLENGTH(id) != P) return false;
    const int *pi = INTEGER(i), *pj = INTEGER(j), *pid = INTEGER(id);
    R_xlen_t r = 0;
    for (int a = 1; a <= n; ++a) {
        for (int b = a; b <= n; ++b, ++r) {
            if (pi[r] != a || pj[r] != b || pid[r] != r + 1) return false;
        }
    }
    return true;
}

// Template columns i, j (integer) and widthprod (double, width[i] * width[j]) for the given
// tile widths; the same values make_template_dat_cpp() produces, preallocated.
// [[Rcpp::export]]
List cmp_template_cols(NumericVector widths, bool with_widthprod) {
    int n = widths.size();
    R_xlen_t P = (R_xlen_t)n * (n + 1) / 2;
    IntegerVector ii(P), jj(P);
    NumericVector wp(with_widthprod ? P : 0);
    int *pi = INTEGER(ii), *pj = INTEGER(jj);
    double *pw = with_widthprod ? REAL(wp) : nullptr;
    const double *w = REAL(widths);
    R_xlen_t r = 0;
    for (int a = 0; a < n; ++a) {
        for (int b = a; b < n; ++b, ++r) {
            pi[r] = a + 1;
            pj[r] = b + 1;
            if (pw) pw[r] = w[a] * w[b];
        }
    }
    if (with_widthprod) return List::create(_["i"] = ii, _["j"] = jj, _["widthprod"] = wp);
    return List::create(_["i"] = ii, _["j"] = jj);
}

// Template row (id) of each (i, j) pair for n tiles: NA unless 1 <= i <= j <= n
// (rows that a join on the template's (i, j) would not match).
// [[Rcpp::export]]
IntegerVector cmp_pix_index(IntegerVector i, IntegerVector j, int n) {
    R_xlen_t m = i.size();
    if (j.size() != m) stop("i and j must have the same length");
    IntegerVector out(m);
    const int *pi = INTEGER(i), *pj = INTEGER(j);
    int *po = INTEGER(out);
    for (R_xlen_t r = 0; r < m; ++r) {
        int a = pi[r], b = pj[r];
        if (a == NA_INTEGER || b == NA_INTEGER || a < 1 || b < a || b > n) { po[r] = NA_INTEGER; continue; }
        // rows before tile a: sum_{k=1}^{a-1} (n - k + 1)
        double before = (double)(a - 1) * n - (double)(a - 1) * (a - 2) / 2.0;
        po[r] = (int)(before + (b - a) + 1);
    }
    return out;
}

// Per-tile coverage of a template-ordered value vector: cov[t] = sum over pixels with i == t
// plus pixels with j == t and i != j. Summed in the order hic_cov()'s reference grouping uses
// (rows with i == t in j order, then rows with j == t in i order), in double.
// [[Rcpp::export]]
NumericVector cmp_tile_cov(NumericVector v, int n) {
    R_xlen_t P = (R_xlen_t)n * (n + 1) / 2;
    if (v.size() != P) stop("value vector does not match the template");
    NumericVector out(n);
    double *acc = REAL(out);
    const double *pv = REAL(v);
    R_xlen_t r = 0;
    for (int a = 0; a < n; ++a)
        for (int b = a; b < n; ++b, ++r) acc[a] += pv[r];
    r = 0;
    for (int a = 0; a < n; ++a)
        for (int b = a; b < n; ++b, ++r) if (b != a) acc[b] += pv[r];
    return out;
}
