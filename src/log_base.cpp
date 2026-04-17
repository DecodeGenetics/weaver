#include "log_base.hpp"

#include <array>
#include <cmath>

#include "logging.hpp"
#include "options.hpp"

namespace weaver
{
double alignment_score_partition_function(double const lambda,
                                          double const match,
                                          double const mismatch,
                                          std::array<double, 4> const nt_freqs)
{
  double partition = 0.0;

  for (int i{0}; i < 4; i++)
  {
    double const i_freq = nt_freqs[i];
    partition += i_freq * i_freq * std::exp(lambda * match); // covers the match case

    // the two following loops cover all mismatch cases
    for (int j{0}; j < i; j++)
      partition += i_freq * nt_freqs[j] * std::exp(lambda * mismatch);

    for (int j{i + 1}; j < 4; j++)
      partition += i_freq * nt_freqs[j] * std::exp(lambda * mismatch);
  }

  if (std::isnan(partition))
  {
    weaver::print_error("overflow error finding the partition");
    std::exit(1);
  }

  return partition;
}

double calculate_log_base()
{
  // Inspired from "recover_log_base" from
  // https://github.com/vgteam/vg/blob/e8b8e0a493a60884acc4f3728f832620852af541/src/aligner.cpp
  Options const & copts = *(weaver::Options::const_instance());
  double const gc_content = copts.gc_content;
  double const match = copts.match;
  double const mismatch = -copts.mismatch;
  double constexpr tol = 1e-12;

  // convert GC content into nucleotide frequencies
  std::array<double, 4> nt_freqs;
  nt_freqs[0] = 0.5 * (1 - gc_content); // A
  nt_freqs[1] = 0.5 * gc_content;       // C
  nt_freqs[2] = 0.5 * gc_content;       // G
  nt_freqs[3] = 0.5 * (1 - gc_content); // T

  // searching for a positive value (because it's a base of a logarithm)
  double lower_bound{};
  double upper_bound{};

  // arbitrary starting point greater than zero
  double lambda{1.0};

  // exponential search for a window containing lambda where total probability is 1
  double partition = alignment_score_partition_function(lambda, match, mismatch, nt_freqs);

  if (partition < 1.0)
  {
    do
    {
      lower_bound = lambda;
      lambda *= 2.0;
      partition = alignment_score_partition_function(lambda, match, mismatch, nt_freqs);
    } while (partition <= 1.0);

    upper_bound = lambda;
  }
  else
  {
    do
    {
      upper_bound = lambda;
      lambda /= 2.0;
      partition = alignment_score_partition_function(lambda, match, mismatch, nt_freqs);
    } while (partition >= 1.0);

    lower_bound = lambda;
  }

  // bisect to find a log base where total probability is 1
  while ((upper_bound / lower_bound) - 1.0 > tol)
  {
    lambda = (lower_bound + upper_bound) / 2.0;

    if (alignment_score_partition_function(lambda, match, mismatch, nt_freqs) < 1.0)
      lower_bound = lambda;
    else
      upper_bound = lambda;
  }

  return (lower_bound + upper_bound) / 2.0;
}

} // namespace weaver
