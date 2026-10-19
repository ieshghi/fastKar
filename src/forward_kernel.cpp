// Hi-C / long-read forward simulation kernels for forward.R.
//
// These reproduce simulate_map() + the sparse lift t(T) L T + symmetrize() + upper-triangle
// extraction of the reference simulate_walks(), the template builder make_template_dat_cpp(),
// and simulate_diploid_map(), to floating-point summation order.
//
// The Hi-C distance model is R's fmm splinefun (stats splines.c, spline_eval, method 3) on
// log(distance), with the power-law extension below the first knot used by distance_density().
// Spline coefficients are passed in from fastKar::splineobj at run time.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <climits>
#include <cstdint>
using namespace Rcpp;

namespace {

// The distance part of the model: knots and coefficients of the fmm spline plus the
// constants of distance_density() and simulate_map().
struct DistModel {
    const double *x, *y, *b, *c, *d;
    R_xlen_t n;
    double lmin;          // first knot, log scale
    double slope;         // spline derivative at the first knot
    double maxval;        // 1e8: above this, 1/d tail
    double maxval_density;// distance_density(maxval)
    int model_type;       // 0 = Hi-C spline, 1 = long-read exponential
    double model;         // read length (model_type 1)
    double lr_denom;      // 10 * min.res^2 (model_type 1)

    DistModel(List z, double lmin_, double slope_, double maxval_, double maxval_density_,
              double model_, double lr_denom_) {
        NumericVector zx = z["x"], zy = z["y"], zb = z["b"], zc = z["c"], zd = z["d"];
        int method = as<int>(z["method"]);
        if (method != 3) stop("forward kernel: only the fmm spline (method 3) is supported");
        x = zx.begin(); y = zy.begin(); b = zb.begin(); c = zc.begin(); d = zd.begin();
        n = zx.size();
        if (n < 1 || zy.size() != n || zb.size() != n || zc.size() != n || zd.size() != n)
            stop("forward kernel: malformed spline coefficients");
        lmin = lmin_; slope = slope_; maxval = maxval_; maxval_density = maxval_density_;
        model = model_;
        model_type = model_ > 0 ? 1 : 0;
        lr_denom = lr_denom_;
    }

    // R's spline_eval for one point, method != 2 (no natural-spline left extrapolation).
    // Interval: the binary search of spline_eval, i.e. the largest i with x[i] <= u (0 if u < x[0]).
    inline double spline(double u, R_xlen_t &hint) const {
        R_xlen_t i = hint;
        if (!(x[i] <= u && (i == n - 1 || u < x[i + 1]))) {
            // bracket [lo, hi) around the interval by galloping from the hint, then bisect as
            // spline_eval does; the result is the same canonical interval as a full search
            R_xlen_t lo = 0, hi = n;
            if (!(u < x[0]) && !std::isnan(u)) {
                if (x[i] <= u) {
                    R_xlen_t step = 1;
                    lo = i;
                    while (lo + step < n && x[lo + step] <= u) { lo += step; step *= 2; }
                    hi = std::min<R_xlen_t>(lo + step, n);
                } else {
                    R_xlen_t step = 1;
                    hi = i;
                    while (hi - step > 0 && u < x[hi - step]) { hi -= step; step *= 2; }
                    lo = std::max<R_xlen_t>(hi - step, 0);
                }
            }
            i = lo;
            R_xlen_t j = hi;
            while (j > i + 1) {
                R_xlen_t k = (i + j) / 2;
                if (u < x[k]) j = k; else i = k;
            }
            hint = i;
        }
        double dx = u - x[i];
        return y[i] + dx * (b[i] + dx * (c[i] + dx * d[i]));
    }

    // distance_density(dist) for dist > 0
    inline double density(double dist, R_xlen_t &hint) const {
        double ldist = std::log(dist);
        double out;
        if (ldist < lmin) {
            R_xlen_t h0 = 0;
            out = spline(lmin, h0) + slope * (ldist - lmin);
        } else {
            out = spline(ldist, hint);
        }
        return std::exp(out);
    }

    // value of one local pixel with width product wp at distance dist, as in simulate_map():
    // returns false when dist == 0 in the Hi-C model (diagonal: diag_spline, filled by the caller).
    // NA distances give 0.
    inline bool value(double wp, double dist, double &v, R_xlen_t &hint) const {
        if (std::isnan(dist)) { v = 0; return true; }
        if (model_type == 0) {
            if (dist > 0) {
                if (dist > maxval) v = wp * maxval_density * maxval / dist;
                else v = wp * density(dist, hint);
            } else if (dist == 0) {
                return false;
            } else {
                v = 0;  // negative distances do not occur (abs())
            }
        } else {
            v = wp / lr_denom * std::exp(-dist / model);
        }
        if (std::isnan(v)) v = 0;
        return true;
    }
};

// pixel id of (p, q), 1-based, p <= q, on an n-tile template (row of make_template_dat_cpp), 0-based
inline int64_t pix0(int64_t p, int64_t q, int64_t n) {
    return (p - 1) * n - (p - 1) * (p - 2) / 2 + (q - p);
}

// deferred dist == 0 contributions, evaluated through diag_spline in one R call
struct Deferred {
    std::vector<int64_t> id;
    std::vector<double> wp, mult;
};

void apply_deferred(const Deferred &df, Function diagf, double *acc) {
    if (df.id.empty()) return;
    NumericVector lwp(df.wp.size());
    for (size_t k = 0; k < df.wp.size(); ++k) lwp[k] = std::log(df.wp[k]);
    NumericVector dv = diagf(lwp);
    if ((size_t) dv.size() != df.wp.size()) stop("diag_spline returned the wrong length");
    for (size_t k = 0; k < df.id.size(); ++k) {
        double v = std::exp(dv[k]);
        if (std::isnan(v)) v = 0;
        acc[df.id[k]] += df.mult[k] * v;
    }
}

} // namespace

//' Template of all pixels (i <= j) of a tiling (moved from a load-time cppFunction)
//'
//' @param A data.frame with a "width" column
//' @return data.frame with columns i, j, value (0), widthprod, id
// [[Rcpp::export]]
DataFrame make_template_dat_cpp(DataFrame A) {
    NumericVector widths = A["width"];
    int64_t l = widths.size();
    int64_t P = l * (l + 1) / 2;
    if (P > INT_MAX) stop("too many tiles for an integer pixel id");
    IntegerVector iv = no_init(P), jv = no_init(P), id = no_init(P);
    NumericVector wp = no_init(P), val(P);
    int64_t k = 0;
    for (int64_t i = 0; i < l; ++i) {
        for (int64_t j = i; j < l; ++j, ++k) {
            iv[k] = i + 1;
            jv[k] = j + 1;
            wp[k] = widths[i] * widths[j];
            id[k] = k + 1;
        }
    }
    return DataFrame::create(_["i"] = iv, _["j"] = jv, _["value"] = val,
                             _["widthprod"] = wp, _["id"] = id);
}

//' distance_density() in C++ (kernel test hook)
// [[Rcpp::export]]
NumericVector forward_density_cpp(NumericVector dist, List z, double lmin, double slope) {
    DistModel m(z, lmin, slope, 1e8, 0, 0, 1);
    NumericVector out(dist.size());
    R_xlen_t hint = 0;
    for (R_xlen_t k = 0; k < dist.size(); ++k) out[k] = m.density(dist[k], hint);
    return out;
}

//' Simulate a set of walks into one reference-pixel accumulator
//'
//' Walk w occupies local rows walk_ptr[w] .. walk_ptr[w+1]-1 (0-based). For each walk the local
//' coordinates are the cumulative sum of the row widths (as the reference cumsum(), including its
//' integer-overflow / NA behaviour), and every pair of local rows a <= b that both map to a target
//' tile (orig non-NA) is simulated and lifted: a == b adds L_aa to (p, p); a < b adds L_ab to
//' (min(p,q), max(p,q)) if p != q and 2 L_ab to (p, p) if p == q. Each walk's map is multiplied by cn[w].
//'
//' @return numeric vector over the n_ref*(n_ref+1)/2 template pixels, in pixel-id order
// [[Rcpp::export]]
NumericVector forward_walks_cpp(IntegerVector walk_ptr, IntegerVector orig, NumericVector width,
                                bool int_width, LogicalVector circular, NumericVector cn, int n_ref,
                                List z, double lmin, double slope, double maxval, double maxval_density,
                                double model, double lr_denom, Function diagf) {
    DistModel m(z, lmin, slope, maxval, maxval_density, model, lr_denom);
    int nw = walk_ptr.size() - 1;
    if (nw < 0 || circular.size() != nw || cn.size() != nw) stop("forward_walks_cpp: bad walk index");
    int64_t P = (int64_t) n_ref * (n_ref + 1) / 2;
    if (P > INT_MAX) stop("too many tiles for an integer pixel id");
    NumericVector accv(P);
    double *acc = accv.begin();
    Deferred df;
    bool overflow = false;
    std::vector<double> mid;
    std::vector<int> inrow;
    for (int w = 0; w < nw; ++w) {
        int r0 = walk_ptr[w], r1 = walk_ptr[w + 1];
        int nr = r1 - r0;
        mid.assign(nr, NA_REAL);
        // cumsum(width) as R does it: integer overflow or NA ends the sum (rest NA)
        double s = 0, totl = NA_REAL;
        bool dead = false;
        for (int r = 0; r < nr; ++r) {
            double wd = width[r0 + r];
            double end = NA_REAL;
            if (!dead) {
                if (std::isnan(wd)) dead = true;
                else {
                    s += wd;
                    if (int_width && s > (double) INT_MAX) { dead = true; overflow = true; }
                    else end = s;
                }
            }
            double start = end - wd + 1;
            mid[r] = (start + end) / 2;
            if (r == nr - 1) totl = end;
        }
        inrow.clear();
        for (int r = 0; r < nr; ++r) if (orig[r0 + r] != NA_INTEGER) inrow.push_back(r);
        bool circ = circular[w] == TRUE;
        double cw = cn[w];
        int ni = inrow.size();
        for (int ia = 0; ia < ni; ++ia) {
            int a = inrow[ia];
            int p = orig[r0 + a];
            if (p < 1 || p > n_ref) stop("orig.id out of range");
            double wa = width[r0 + a], ma = mid[a];
            R_xlen_t hint = 0;
            for (int ib = ia; ib < ni; ++ib) {
                int b = inrow[ib];
                int q = orig[r0 + b];
                if (q < 1 || q > n_ref) stop("orig.id out of range");
                double wp = wa * width[r0 + b];
                double dist = std::fabs(mid[b] - ma);
                if (circ) {
                    double d2 = totl - dist;
                    if (std::isnan(d2)) dist = d2; else if (d2 < dist) dist = d2;
                }
                int64_t id;
                double mult;
                if (a == b || p == q) { id = pix0(p, p, n_ref); mult = (a == b) ? 1.0 : 2.0; }
                else if (p < q) { id = pix0(p, q, n_ref); mult = 1.0; }
                else { id = pix0(q, p, n_ref); mult = 1.0; }
                double v;
                if (m.value(wp, dist, v, hint)) {
                    acc[id] += cw * (mult * v);
                } else {
                    df.id.push_back(id); df.wp.push_back(wp); df.mult.push_back(cw * mult);
                }
            }
        }
    }
    if (overflow) Rf_warning("integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))'");
    apply_deferred(df, diagf, acc);
    return accv;
}

//' Assemble the summed map: i, j, id and value = depth * (acc [+ widthprod * rho * cn_i * cn_j])
// [[Rcpp::export]]
List forward_assemble_cpp(NumericVector width, NumericVector cn, NumericVector acc, double depth,
                          double rho, bool interchr) {
    int64_t l = width.size();
    int64_t P = l * (l + 1) / 2;
    if ((int64_t) acc.size() != P) stop("forward_assemble_cpp: accumulator length mismatch");
    if (interchr && cn.size() != l) stop("forward_assemble_cpp: cn length mismatch");
    IntegerVector iv = no_init(P), jv = no_init(P), id = no_init(P);
    NumericVector val = no_init(P);
    int64_t k = 0;
    for (int64_t i = 0; i < l; ++i) {
        for (int64_t j = i; j < l; ++j, ++k) {
            iv[k] = i + 1; jv[k] = j + 1; id[k] = k + 1;
            double v = acc[k];
            if (interchr) v = v + width[i] * width[j] * rho * (cn[i] * cn[j]);
            val[k] = depth * v;
        }
    }
    List out = List::create(_["i"] = iv, _["j"] = jv, _["id"] = id, _["value"] = val);
    out.attr("row.names") = IntegerVector::create(NA_INTEGER, -(int) P);
    out.attr("class") = "data.frame";
    return out;
}

//' simulate_diploid_map() values (Hi-C model) in pixel-id order, before is.na and depth:
//' 2 * the single-copy map on same-chromosome pixels, plus the interchromosomal term if asked.
// [[Rcpp::export]]
NumericVector forward_diploid_cpp(IntegerVector chrom, NumericVector mid, NumericVector width,
                                  List z, double lmin, double slope, double maxval, double maxval_density,
                                  bool interchr, double rho, double depth, Function diagf) {
    DistModel m(z, lmin, slope, maxval, maxval_density, 0, 1);
    int64_t l = width.size();
    int64_t P = l * (l + 1) / 2;
    if (P > INT_MAX) stop("too many tiles for an integer pixel id");
    NumericVector outv(P);
    double *out = outv.begin();
    Deferred df;
    int64_t k = 0;
    for (int64_t i = 0; i < l; ++i) {
        R_xlen_t hint = 0;
        for (int64_t j = i; j < l; ++j, ++k) {
            if (chrom[i] != chrom[j] || chrom[i] == NA_INTEGER) continue;
            double wp = width[i] * width[j];
            double dist = std::fabs(mid[j] - mid[i]);
            if (std::isnan(dist)) continue;
            if (dist > 0) {
                if (dist > maxval) out[k] = 2 * wp * maxval_density * maxval / dist;
                else out[k] = 2 * wp * m.density(dist, hint);
            } else if (dist == 0) {
                df.id.push_back(k); df.wp.push_back(wp); df.mult.push_back(2.0);
            }
        }
    }
    apply_deferred(df, diagf, out);
    k = 0;
    for (int64_t i = 0; i < l; ++i) {
        for (int64_t j = i; j < l; ++j, ++k) {
            double v = out[k];
            if (interchr) v = v + width[i] * width[j] * rho * 4;
            if (std::isnan(v)) v = 0;
            out[k] = depth * v;
        }
    }
    return outv;
}
