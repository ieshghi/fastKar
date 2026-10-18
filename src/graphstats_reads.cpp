#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Exact sweep replacement for the per-walk foverlaps() + paste() in reads_fromwalk()
// (R/graphstats.R, fastKar 0a5497e). For each walk the reference
//   1. lays the nodes end to end (start = cumsum(width) - width + 1); a circular walk is
//      extended by a second copy shifted by the walk length and clipped at
//      total + eff - 1, with eff = min(readL, total) (eff = readL for a linear walk);
//   2. cuts the read-start domain ([2 - eff, total] linear, [1, total] circular) at
//      every start - eff + 1 and end + 1, and per bin of starts takes the nodes that the
//      representative read [s, s + eff - 1] overlaps (closed intervals; integer
//      coordinates), in node order, dropping bins without hits;
//   3. emits, per bin, the word of the kept (non-tiny) node ids joined by '|' and the word
//      of the reverse complement, each with multiplicity = bin width;
//   4. sums multiplicities by word and orders by (-nums, word) (C-locale byte order).
// Across walks, words are concatenated in walk order and the first occurrence of each word
// keeps its place, with nums summed over walks.
// Here the hits of a bin are a contiguous node range [lo, hi] found by a two-pointer sweep
// (starts and ends both increase along the laid-out walk); bins are first reduced to
// distinct kept-node ranges, so each distinct word is built once.
// Arithmetic on positions follows the reference expression order exactly; the R wrapper only
// uses this path for integer-valued readL, where foverlaps' double end-adjustment
// ((end + eps) * (1 + eps)) cannot change any overlap.

namespace {

struct RangeHash {
  size_t operator()(const std::pair<int, int>& p) const {
    return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned int)p.second);
  }
};

}  // namespace

// widths, fwd_tok, rev_tok, keep: flattened over walks (lens gives each walk's length).
// fwd_tok[k] = as.character(x[k]), rev_tok[k] = as.character(-x[k]) for snode id x[k].
// [[Rcpp::export]]
List reads_fromwalk_cpp(IntegerVector widths, IntegerVector lens, LogicalVector circular,
                        double readL, CharacterVector fwd_tok, CharacterVector rev_tok,
                        LogicalVector keep) {
  const int nw = lens.size();
  if (circular.size() != nw) stop("circular and lens differ in length");

  std::vector<std::string> gwords;  // global, first-appearance order
  std::vector<long double> gnums;
  std::unordered_map<std::string, int> gidx;

  int off = 0;
  for (int w = 0; w < nw; w++) {
    const int m = lens[w];
    if (m <= 0) stop("reads_fromwalk_cpp: empty walk");
    // ---- 1. lay out the walk ----
    std::vector<double> st, en;   // laid-out (possibly extended) nodes
    std::vector<int> src;         // index into the flat arrays
    st.reserve(2 * m); en.reserve(2 * m); src.reserve(2 * m);
    double cum = 0;
    for (int k = 0; k < m; k++) {
      int wd = widths[off + k];
      if (wd == NA_INTEGER || wd < 1) stop("reads_fromwalk_cpp: node widths must be >= 1");
      cum += wd;
      st.push_back(cum - wd + 1);
      en.push_back(cum);
      src.push_back(off + k);
    }
    const double total = cum;
    double eff, dlo, dhi;
    const bool circ = circular[w] == TRUE;
    if (circ) {
      eff = std::min(readL, total);
      dlo = 1;
      dhi = total;
      const double max_read_end = total + eff - 1;
      for (int k = 0; k < m; k++) {
        double s2 = st[k] + total, e2 = en[k] + total;
        if (!(s2 <= max_read_end)) continue;
        if (e2 > max_read_end) e2 = max_read_end;
        st.push_back(s2); en.push_back(e2); src.push_back(off + k);
      }
    } else {
      eff = readL;
      dlo = st[0] - eff + 1;
      dhi = total;
    }
    const int mx = (int)st.size();

    // ---- 2. cuts ----
    std::vector<double> cuts;
    cuts.reserve(2 * mx + 2);
    const double dhi1 = dhi + 1;
    auto add_cut = [&](double c) { if (c >= dlo && c <= dhi1) cuts.push_back(c); };
    add_cut(dlo);
    add_cut(dhi1);
    for (int k = 0; k < mx; k++) { add_cut(st[k] - eff + 1); add_cut(en[k] + 1); }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

    // kept-node prefix counts over the laid-out nodes
    std::vector<int> kpre(mx + 1, 0), kidx;
    kidx.reserve(mx);
    for (int k = 0; k < mx; k++) {
      bool kp = keep[src[k]] == TRUE;
      kpre[k + 1] = kpre[k] + (kp ? 1 : 0);
      if (kp) kidx.push_back(src[k]);
    }

    // ---- sweep bins -> distinct kept ranges (a, b) half-open in kept space ----
    std::unordered_map<std::pair<int, int>, int, RangeHash> ridx;
    std::vector<std::pair<int, int>> ranges;
    std::vector<long double> rmult;
    int lo = 0, hi = -1;  // lo: first node with end >= rs; hi: last node with start <= re
    for (size_t b = 0; b + 1 < cuts.size(); b++) {
      const double rs = cuts[b];
      const double ns = cuts[b + 1];
      const double mult = ns - rs;
      const double re = rs + eff - 1;
      while (lo < mx && en[lo] < rs) lo++;
      while (hi + 1 < mx && st[hi + 1] <= re) hi++;
      if (lo > hi) continue;  // nomatch = 0: bin dropped
      std::pair<int, int> key(kpre[lo], kpre[hi + 1]);
      auto it = ridx.find(key);
      if (it == ridx.end()) {
        ridx.emplace(key, (int)ranges.size());
        ranges.push_back(key);
        rmult.push_back(mult);
      } else {
        rmult[it->second] += mult;
      }
    }

    // ---- 3. words (forward, then reverse complement) summed by word ----
    std::unordered_map<std::string, int> widx;
    std::vector<std::string> words;
    std::vector<long double> nums;
    auto add_word = [&](std::string&& s, long double v) {
      auto it = widx.find(s);
      if (it == widx.end()) {
        widx.emplace(s, (int)words.size());
        words.push_back(std::move(s));
        nums.push_back(v);
      } else {
        nums[it->second] += v;
      }
    };
    for (int pass = 0; pass < 2; pass++) {
      for (size_t r = 0; r < ranges.size(); r++) {
        const int a = ranges[r].first, bnd = ranges[r].second;
        std::string s;
        if (pass == 0) {
          for (int q = a; q < bnd; q++) {
            if (q > a) s.push_back('|');
            SEXP t = STRING_ELT(fwd_tok, kidx[q]);
            s.append(CHAR(t), LENGTH(t));
          }
        } else {
          for (int q = bnd - 1; q >= a; q--) {
            if (q < bnd - 1) s.push_back('|');
            SEXP t = STRING_ELT(rev_tok, kidx[q]);
            s.append(CHAR(t), LENGTH(t));
          }
        }
        add_word(std::move(s), rmult[r]);
      }
    }

    // ---- 4. order by (-nums, word) ----
    std::vector<int> ord(words.size());
    for (size_t i = 0; i < ord.size(); i++) ord[i] = (int)i;
    std::sort(ord.begin(), ord.end(), [&](int x, int y) {
      if (nums[x] != nums[y]) return nums[x] > nums[y];
      return words[x] < words[y];
    });
    for (int i : ord) {
      double v = (double)nums[i];
      auto it = gidx.find(words[i]);
      if (it == gidx.end()) {
        gidx.emplace(words[i], (int)gwords.size());
        gwords.push_back(words[i]);
        gnums.push_back(v);
      } else {
        gnums[it->second] += v;
      }
    }
    off += m;
  }
  if (off != widths.size()) stop("reads_fromwalk_cpp: lens do not sum to the number of nodes");

  const int nout = (int)gwords.size();
  CharacterVector ow(nout);
  NumericVector on(nout);
  for (int i = 0; i < nout; i++) {
    ow[i] = Rf_mkCharLenCE(gwords[i].data(), (int)gwords[i].size(), CE_NATIVE);
    on[i] = (double)gnums[i];
  }
  return List::create(_["words"] = ow, _["nums"] = on);
}
