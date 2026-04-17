#pragma once

#include <gfa.h>

#include "gfa.hpp"

namespace weaver
{
/*!
 * @brief Given a vertex \a v, get a pointer to the first arc element.
 *
 * @param[in] gfa graph.
 * @param[in] v vertex index.
 *
 * @returns Pointer to the first arc element. If there are no arcs, then arc_begin() and arc_end() will return the same
 * pointer.
 *
 * @see arc_end(GFA const & gfa, uint32_t const v)
 * @see arc_end(GFA const & gfa, uint32_t const v, gfa_arc_t const * arc_begin)
 */
gfa_arc_t const * arc_begin(GFA const & gfa, uint32_t const v);

/*!
 * @brief Given a vertex \a v, get a pointer behind the last arc element.
 *
 * @param[in] gfa graph.
 * @param[in] v vertex index.
 *
 * @returns Pointer behind the last arc element. If there are no arcs,
 * then arc_begin() and arc_end(GFA const & gfa, uint32_t const v) will return the same pointer.
 *
 * @see arc_begin()
 * @see arc_end(GFA const & gfa, uint32_t const v, gfa_arc_t const * arc_begin)
 */
gfa_arc_t const * arc_end(GFA const & gfa, uint32_t const v);

/*!
 * @brief Given a vertex \a v and a pointer to the first arc, get a pointer behind the last arc element.
 *
 * @see arc_begin() for getting the first arc.
 * @see arc_end(GFA const & gfa, uint32_t const v) .
 */
gfa_arc_t const * arc_end(GFA const & gfa, uint32_t const v, gfa_arc_t const * arc_begin);

/*!
 * @brief Checks if arcs \a a and \a b are complement of each other.
 *
 * @details
 * Two arcs are complements of each other if you can follow them and end up on the same segment as you started.
 *
 * @note
 * Neither of the two parameters can be a nullptr.
 *
 * @param[in] a first arc to compare.
 * @param[in] b second arc to compare.
 *
 * @see is_same_arc() checks if two arcs are the same.
 */
bool is_complement_arc(gfa_arc_t const * a, gfa_arc_t const * b);

/*!
 * @brief Checks if arcs \a a and \a b are the same.
 *
 * @note
 * Neither of the two parameters can be a nullptr.
 *
 * @param[in] a first arc to compare.
 * @param[in] b second arc to compare.
 *
 * @see is_complement_arc() checks if two arcs are complement of each other.
 */
bool is_same_arc(gfa_arc_t const * a, gfa_arc_t const * b);

/*!
 * @brief Returns a pretty string showing the contents of an arc.
 *
 * @details
 * Useful for debugging purposes.
 *
 * @param[in] arc Arc to show the contents of.
 *
 * @returns A pretty string with the contents of \a arc .
 */
std::string arc_to_string(gfa_arc_t const & arc);

} // namespace weaver
