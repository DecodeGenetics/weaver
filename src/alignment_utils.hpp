#pragma once
/*!
 * @file alignment_utils.hpp
 *
 * @brief Utility functions for alignments
 */

#include <vector>

namespace weaver
{
/*!
 * @brief Add two values in log space.
 *
 * @details
 * Return the log of the sum of two log-transformed values without taking them out of log-space.
 *
 * Source: https://stackoverflow.com/questions/778047/we-know-log-add-but-how-to-do-log-subtract
 *
 * @param[in] x One of the values to add.
 * @param[in] y Other value to add.
 *
 * @returns Sum of the two values in log space.
 */
double add_log(double const x, double const y);

/*!
 * @brief Add value in log space to another referenced value.
 *
 * @details
 * No value is returned, only the referenced value is updated instead.
 * add_log_ref(x,y) has the same effect as x = add_log(x,y)
 *
 * @param[in,out] x Add to this referenced value.
 * @param[in]     y Add this value.
 *
 * @see add_log(x,y)
 */
void add_log_ref(double & x, double const y);

/*!
 * @brief Calculate the sum of multiple values
 *
 * @details
 * Inspired from: https://stackoverflow.com/questions/65233445/how-to-calculate-sums-in-log-space-without-underflow
 */
template <typename T>
double log_sum(std::vector<T> const & vals);

/*!
 * @brief Calculate the sum of multiple values
 *
 * @details
 * base will be added to each value first.
 *
 * @details
 * Inspired from: https://stackoverflow.com/questions/65233445/how-to-calculate-sums-in-log-space-without-underflow
 */
template <typename T>
double log_sum_base(std::vector<T> const & vals, double const base);

/*!
 * @brief Calculate the sum of multiple values with largest value first
 *
 * @details
 * Inspired from: https://stackoverflow.com/questions/65233445/how-to-calculate-sums-in-log-space-without-underflow
 */
template <typename T>
double log_sum_max_first(std::vector<T> const & vals);

/*!
 * @brief Subtract two values in log-space.
 *
 * @details
 * Calculate log of the difference of two log-transformed values without taking them out of log-space.
 *
 * Source: https://stackoverflow.com/questions/778047/we-know-log-add-but-how-to-do-log-subtract
 *
 * @param[in] x value to subtract.
 * @param[in] y value to subtract.
 *
 * @results The result of subtracting the two values in log-space.
 */
double subtract_log(double const x, double const y);

/*!
 * @brief Convert a PHRED quality score to probability of error.
 *
 * @param[in] phred PHRED quality score.
 *
 * @returns The Probability of error.
 */
double phred_to_prob(double const phred);

/*!
 * @brief Convert probability of error to a PHRED quality score.
 *
 * The function returns a float, while PHREDs are normally represented as integers.
 * Use std::round(prob_to_phred(prob)) to get an integer.
 *
 * @param[in] prob_error Probability of an error.
 *
 * @return PHRED scaled value of the probability.
 */
double prob_to_phred(double prob_error);

template <typename T>
inline T & get_element_reference(int index, std::vector<T> & vec, T & empty)
{
  return (index < 0 || index >= static_cast<int>(vec.size())) ? empty : vec[index];
}

} // namespace weaver
