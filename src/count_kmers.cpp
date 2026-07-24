#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
using namespace Rcpp;
// Hash function for a vector<int>
struct VectorHash {
  std::size_t operator()(const std::vector<int>& x) const noexcept {
    std::size_t seed = x.size();

    for (int value : x) {
      seed ^= std::hash<int>{}(value)
            + 0x9e3779b9
            + (seed << 6)
            + (seed >> 2);
    }

    return seed;
  }
};

struct KmerCount {
  std::vector<int> kmer;
  int count;
};
// [[Rcpp::export]]
DataFrame count_kmers_cpp(
    List snode_ids,
    LogicalVector circular,
    int k
) {
  if (k < 1)
    stop("k must be a positive integer");

  int nwalks = snode_ids.size();

  if (circular.size() != nwalks)
    stop("snode_ids and circular must have the same length");

  // Estimate the number of k-mer occurrences for hash-table allocation
  std::size_t estimated_occurrences = 0;

  for (int i = 0; i < nwalks; ++i) {
    IntegerVector x = snode_ids[i];
    int n = x.size();

    if (n == 0)
      continue;

    if (circular[i] == TRUE) {
      if (k <= n + 1)
        estimated_occurrences += 2ULL * n;
    } else {
      if (n >= k)
        estimated_occurrences += 2ULL * (n - k + 1);
    }
  }

  std::unordered_map<
    std::vector<int>,
    int,
    VectorHash
  > counts;

  counts.reserve(estimated_occurrences);

  std::vector<int> forward(k);
  std::vector<int> reverse_complement(k);

  for (int i = 0; i < nwalks; ++i) {
    IntegerVector x = snode_ids[i];
    int n = x.size();

    if (n == 0)
      continue;

    bool is_circular = circular[i] == TRUE;
    int nstarts;

    if (is_circular) {
      // ABC may produce ABCA, but not ABCAB
      if (k > n + 1)
        continue;

      nstarts = n;
    } else {
      if (n < k)
        continue;

      nstarts = n - k + 1;
    }

    for (int start = 0; start < nstarts; ++start) {
      for (int j = 0; j < k; ++j) {
        int index = start + j;

        // Circular walks repeat as needed
        if (is_circular)
          index %= n;

        int value = x[index];

        forward[j] = value;
        reverse_complement[k - j - 1] = -value;
      }

      ++counts[forward];
      ++counts[reverse_complement];
    }
  }

  std::vector<KmerCount> result;
  result.reserve(counts.size());

  for (const auto& entry : counts) {
    result.push_back({entry.first, entry.second});
  }

  std::sort(
    result.begin(),
    result.end(),
    [](const KmerCount& a, const KmerCount& b) {
      return a.count > b.count;
    }
  );

  int n_unique = result.size();

  List output(k + 2);
  CharacterVector output_names(k + 2);

  // One integer output column per k-mer position
  for (int j = 0; j < k; ++j) {
    IntegerVector column(n_unique);

    for (int i = 0; i < n_unique; ++i)
      column[i] = result[i].kmer[j];

    output[j] = column;
    output_names[j] = "node" + std::to_string(j + 1);
  }

  IntegerVector count_column(n_unique);
  CharacterVector kmer_column(n_unique);

  for (int i = 0; i < n_unique; ++i) {
    count_column[i] = result[i].count;

    std::ostringstream stream;

    for (int j = 0; j < k; ++j) {
      if (j > 0)
        stream << "|";

      stream << result[i].kmer[j];
    }

    kmer_column[i] = stream.str();
  }

  output[k] = count_column;
  output_names[k] = "count";

  output[k + 1] = kmer_column;
  output_names[k + 1] = "kmer";

  output.attr("names") = output_names;
  output.attr("class") = "data.frame";
  output.attr("row.names") = IntegerVector::create(NA_INTEGER, -n_unique);

  return DataFrame(output);
}
