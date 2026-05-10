#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <cstdint>
#include <cstdio>
using namespace Rcpp;

// Canonical-karyotype hash, C++ replacement for hash_snodelist (R/graphstats.R).
// Behavior matches v1 for the partition relation it induces:
//   1. For each circular walk: rotate to its lex-min rotation (Booth).
//   2. For each walk: replace with min(walk, rc(walk)) lex-wise, where
//      rc(w)[i] = -w[n-1-i].   (This matches sort_snodes; see note below.)
//   3. Build a multiset including each canonicalized walk AND its rc (the
//      "doubling" from hash_snodelist) annotated with its circular flag.
//   4. Sort the multiset lex.
//   5. FNV-1a 64-bit hash over the bytes; return as a 16-char hex string.
//
// Note on canonicalization: for circular walks, R's sort_snodes computes
// rc(Booth(W)) and compares against Booth(W), rather than Booth(rc(W))
// vs Booth(W). For a strict canonical form under cyclic+RC equivalence
// you'd want the latter. We deliberately replicate the former here so
// the v1/v2 partition equivalence check in the benchmark passes; tightening
// the canonicalization is a separate (correctness) PR.

// ----------------------------------------------------------------------------
// Booth's algorithm: index of the lex-min cyclic rotation of s.
// Linear time, no auxiliary suffix array.
// ----------------------------------------------------------------------------
static int booth_rotation_start(const std::vector<int>& s) {
  const int n = (int)s.size();
  if (n <= 1) return 0;
  std::vector<int> s2(2 * n);
  for (int i = 0; i < n; i++) { s2[i] = s[i]; s2[i + n] = s[i]; }
  int i = 0, j = 1, k = 0;
  while (i < n && j < n && k < n) {
    int a = s2[i + k], b = s2[j + k];
    if (a == b) {
      k++;
    } else if (a > b) {
      i = i + k + 1;
      if (i <= j) i = j + 1;
      k = 0;
    } else {
      j = j + k + 1;
      if (j <= i) j = i + 1;
      k = 0;
    }
  }
  int start = std::min(i, j);
  if (start >= n) start -= n;
  return start;
}

static void apply_booth_rotation(std::vector<int>& s) {
  const int n = (int)s.size();
  if (n <= 1) return;
  int start = booth_rotation_start(s);
  if (start > 0 && start < n) {
    std::rotate(s.begin(), s.begin() + start, s.end());
  }
}

// rc(s)[i] = -s[n-1-i]
static std::vector<int> rc_walk(const std::vector<int>& s) {
  std::vector<int> r(s.size());
  for (size_t i = 0; i < s.size(); i++) r[i] = -s[s.size() - 1 - i];
  return r;
}

// Replace s with min(s, rc(s)) lex-wise.
static void rc_canonicalize_inplace(std::vector<int>& s) {
  std::vector<int> r = rc_walk(s);
  if (r < s) s = std::move(r);
}

// ----------------------------------------------------------------------------
// FNV-1a 64-bit byte stream hash.
// ----------------------------------------------------------------------------
struct FNV1a64 {
  uint64_t h;
  FNV1a64() : h(0xcbf29ce484222325ULL) {}
  void update_byte(uint8_t b) { h = (h ^ (uint64_t)b) * 0x100000001b3ULL; }
  void update_int(int v) {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(&v);
    for (size_t i = 0; i < sizeof(v); i++) update_byte(p[i]);
  }
};

// ----------------------------------------------------------------------------
// Public entry point.
// ----------------------------------------------------------------------------

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
std::string hash_karyotype_cpp(List snode_id, LogicalVector circular) {
  const int K = snode_id.size();
  if (K == 0) return std::string("0000000000000000");
  if (K != circular.size())
    Rcpp::stop("hash_karyotype_cpp: snode_id and circular have different lengths");

  // 1. Pull walks; canonicalize each (Booth on circular, then RC-min).
  std::vector<std::vector<int>> walks(K);
  std::vector<bool> circ(K);
  for (int k = 0; k < K; k++) {
    IntegerVector iv = snode_id[k];
    walks[k].assign(iv.begin(), iv.end());
    circ[k] = (bool)circular[k];
    if (circ[k]) apply_booth_rotation(walks[k]);
    rc_canonicalize_inplace(walks[k]);
  }

  // 2. Build doubled multiset: (walk, circ) + (rc(walk), circ).
  std::vector<std::pair<std::vector<int>, bool>> entries;
  entries.reserve(2 * K);
  for (int k = 0; k < K; k++) {
    entries.emplace_back(walks[k], circ[k]);
    entries.emplace_back(rc_walk(walks[k]), circ[k]);
  }

  // 3. Sort lex.
  std::sort(entries.begin(), entries.end(),
    [](const std::pair<std::vector<int>, bool>& a,
       const std::pair<std::vector<int>, bool>& b) {
      if (a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });

  // 4. Hash.
  FNV1a64 hasher;
  for (const auto& e : entries) {
    hasher.update_byte(0xFE);                       // entry separator
    hasher.update_byte(e.second ? (uint8_t)'C' : (uint8_t)'L');
    for (int v : e.first) hasher.update_int(v);
  }

  // 5. Format as 16-char lowercase hex.
  char buf[17];
  std::snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)hasher.h);
  return std::string(buf);
}
