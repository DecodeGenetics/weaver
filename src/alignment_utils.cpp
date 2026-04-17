#include <algorithm> // std::max_element
#include <cassert>
#include <cmath> // std::exp, std::log1p
#include <cstdint>
#include <limits>
#include <vector>

namespace weaver
{
double add_log(double const x, double const y)
{
  assert(not std::isinf(x));
  assert(not std::isinf(y));

  return x > y ? x + std::log1p(std::exp(y - x)) : y + std::log1p(std::exp(x - y));
}

void add_log_ref(double & x, double const y)
{
  assert(not std::isinf(x));
  assert(not std::isinf(y));

  x = x > y ? x + std::log1p(std::exp(y - x)) : y + std::log1p(std::exp(x - y));
}

template <typename T>
double log_sum(std::vector<T> const & vals)
{
  assert(vals.size() > 0);
  auto max_it = std::max_element(vals.begin(), vals.end());
  double sum_exp{0.0};

  for (auto it = vals.begin(); it != vals.end(); ++it)
  {
    if (it == max_it)
      continue;

    sum_exp += std::exp(*it - *max_it);
  }

  return *max_it + std::log1p(sum_exp);
}

template <typename T>
double log_sum_base(std::vector<T> const & vals, double const base)
{
  assert(vals.size() > 0);

  if (base >= 0)
  {
    auto max_it = std::max_element(vals.begin(), vals.end());
    double sum_exp{0.0};

    for (auto it = vals.begin(); it != vals.end(); ++it)
    {
      if (it == max_it)
        continue;

      sum_exp += std::exp(static_cast<double>(*it) * base - static_cast<double>(*max_it) * base);
    }

    return static_cast<double>(*max_it) * base + std::log1p(sum_exp);
  }
  else
  {
    auto max_it = std::min_element(vals.begin(), vals.end());
    double sum_exp{0.0};

    for (auto it = vals.begin(); it != vals.end(); ++it)
    {
      if (it == max_it)
        continue;

      sum_exp += std::exp(static_cast<double>(*it) * base - static_cast<double>(*max_it) * base);
    }

    return static_cast<double>(*max_it) * base + std::log1p(sum_exp);
  }
}

template <typename T>
double log_sum_max_first(std::vector<T> const & vals)
{
  assert(vals.size() > 0);
  assert(vals.begin() == std::max_element(vals.begin(), vals.end()));

  double sum_exp{0.0};

  for (int i{1}; i < static_cast<int>(vals.size()); ++i)
    sum_exp += std::exp(vals[i] - vals[0]);

  return vals[0] + std::log1p(sum_exp);
}

double subtract_log(double const x, double const y)
{
  assert(y <= x);
  assert(not std::isinf(y));
  return x + std::log1p(-std::exp(y - x));
}

double phred_to_prob(double const phred)
{
  return std::pow(10.0, -phred / 10.0);
}

double prob_to_phred(double prob_error)
{
  return -10.0 * std::log10(prob_error);
}

template double log_sum(std::vector<double> const & vals);
template double log_sum(std::vector<float> const & vals);
template double log_sum(std::vector<int> const & vals);
template double log_sum(std::vector<int64_t> const & vals);

template double log_sum_max_first(std::vector<double> const & vals);

template double log_sum_base(std::vector<double> const & vals, double base);
template double log_sum_base(std::vector<uint8_t> const & vals, double base);
} // namespace weaver
