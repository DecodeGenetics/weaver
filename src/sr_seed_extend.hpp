#pragma once

#include <cstdint>

#include "gfa.hpp"
#include "icu.hpp"
#include "sr_seed.hpp"

namespace weaver
{
/*!
 * @brief Extend a seed if \a to is exactly \a distance from its end, otherwise do nothing and return false.
 *
 * @details
 * This function essentially wraps is_exact_distance() and Seed::extend(). Since the logic of extending can be
 * a bit tricky, the main motivation of this function is to make extending a little simpler.
 *
 * @param[in] gfa GFA graph
 * @param[in,out] seed Extend from this seed, will be modified iff match is found.
 * @param[in] icu The "I see you" index.
 * @param[in] to To this point
 * @param[in] read_value The read minimizer index value at \a to.
 * @param[in] distance The exact distance to check for.
 *
 * @return True iff the seed was extended.
 *
 * @see seed_extend_if_approximate_distance() if an approximate distance should be enough to extend.
 * @see is_extact_distance().
 */
bool seed_extend_if_exact_distance(
  GFA const & gfa, SRSeed & seed, T_icu const & icu, uint64_t to, uint64_t read_value, int distance);

/*!
 * @brief Extend a seed if \a to is approximately \a distance from its end, otherwise do nothing and return false.
 *
 * @param[in] gfa GFA graph
 * @param[in,out] seed Extend from this seed, will be modified iff match is found.
 * @param[in] icu The "I see you" index.
 * @param[in] to To this point
 * @param[in] read_value The read minimizer index value at \a to.
 * @param[in] distance The approximate distance to check for.
 *
 * @returns True iff the seed was extended.
 *
 * @see seed_extend_if_exact_distance() if the distance should be exactly \a distance .
 */
bool seed_extend_if_approximate_distance(
  GFA const & gfa, SRSeed & seed, T_icu const & icu, uint64_t to, uint64_t read_value, int distance);

} // namespace weaver
