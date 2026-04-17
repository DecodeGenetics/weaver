#pragma once
/*!
 * @file icu.hpp
 * @brief Defines functions to create an "I see you" index.
 */

#include <cstdint> // uint64_t
#include <vector>  // std::vector

#include <parallel_hashmap/phmap_fwd_decl.h>

#include "gfa.hpp" // GFA

#include "gfa.h" // gfa_arc_t

namespace weaver
{
/*!
 * @brief Type definition for the ICU ("I see you") index.
 *
 * The index keys are 64 bits composed of:<br>
 *   from_vertex (32 bits) | to_vertex (32 bits)<br>
 * which is the same as<br>
 *   from_rid (31 bits) | from_strand (1 bit) | to_rid (31 bits) | to_strand (1 bit)<br>
 * and the values of the index are the shortest distance between those positions
 * in the graph.
 */
using T_icu = phmap::flat_hash_map<uint64_t, int>;

/*!
 * @brief Estimates the shortest distance between \a from and \a to in the \a gfa graph.
 *
 * @param[in] gfa graph
 * @param[in] icu The I see you index.
 * @param[in] from position to start at.
 * @param[in] to position to end at.
 *
 * @returns shortest distance. -1 if unreachable or far (not defined what far is, but at least 1501 bp).
 */
int estimate_shortest_distance(GFA const & gfa, //
                               T_icu const & icu,
                               uint64_t const from,
                               uint64_t const to);

/*!
 * @brief Finds if an exact distance is between two (minimizer) values in a graph.
 *
 * Any passed arcs will be added to the \a new_arcs vector.
 *
 * @param[in] gfa graph.
 * @param[in] icu The I see you index. See T_icu
 * @param[in] from from this position value.
 * @param[in] to to this position value.
 * @param[out] new_arcs Passed/crossed arcs will be added to this vector.
 * @param[in] distance the exact distance to check for.
 *
 * @returns true if the exact distance is found, otherwise false.
 */
bool is_exact_distance(GFA const & gfa, //
                       T_icu const & icu,
                       uint64_t from,
                       uint64_t to,
                       std::vector<gfa_arc_t const *> & new_arcs,
                       int distance);

/*!
 * @brief Finds if an approximate distance is between two (minimizer) values in a graph.
 *
 * @param[in] gfa graph.
 * @param[in] icu The "I see you" index.
 * @param[in] from position to start at.
 * @param[in] to position to end at.
 * @param[out] new_arcs Passed/crossed arcs will be added to this vector.
 * @param[in] distance distance to check for.
 * @param[in] t threshold of the distance to allow. It is expected that \a t is less or equal to \a distance .
 *
 * @returns The distance between \a from and \a to.
 * @retval std::numeric_limits<int>::min() if the distance is outside of the approximate range
 *
 * @see is_exact_distance() if the distance between \a from and \a to needs to be exact.
 * @see is_within_distance()
 */
int approximate_gfa_distance(GFA const & gfa,
                             T_icu const & icu,
                             uint64_t from,
                             uint64_t to,
                             std::vector<gfa_arc_t const *> & new_arcs,
                             int distance,
                             int t);

/*!
 * @brief Checks if the distance in the graph between positions \a from and \a to is \a max_distance or shorter.
 *
 * The strand of \a from and \a to are checked, the two positions are only considered to be within a distance
 * if they are reachable to each other without changing the strand.
 *
 * @param[in] gfa graph.
 * @param[in] icu The "I see you" index.
 * @param[in] from position to start at.
 * @param[in] to position to end at.
 * @param[in] max_distance maxmimum distance between the points.
 *
 * @returns true iff the distance between from and to is \a max_distance or less.
 *
 * @see is_exact_distance() if the distance between \a from and \a to needs to be exact.
 * @see approximate_gfa_distance() if some distance threshold can be allowed.
 */
bool is_within_distance(GFA const & gfa, T_icu const & icu, uint64_t from, uint64_t to, int max_distance);

/*!
 * @brief Creates an "I see you"/icu index.
 *
 * @param[in] gfa graph.
 * @param[in] max_distance Maxmimum distance between vertices that will be added to the index.
 *
 * @returns The ICU index.
 *
 * @see T_icu has description of the ICU index.
 */
T_icu make_icu_index(GFA const & gfa, int const max_distance);

} // namespace weaver
