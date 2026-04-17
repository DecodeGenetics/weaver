#pragma once
/*!
 * @file sketch_value.hpp
 * @brief Defines and implements functions to extract from sketch values.
 */
#include <cstdint>

namespace weaver
{
class GFA;

//! Get the reference ID from a sketch value.
inline int sketch_value_rid(uint64_t const sketch_value)
{
  return sketch_value >> 32;
}

//! Get the position from a sketch value.
inline int sketch_value_pos(uint64_t const sketch_value)
{
  return static_cast<uint32_t>(sketch_value) >> 1;
}

//! Get the strand from a sketch value. 0=forward, 1=reverse.
inline bool sketch_value_strand(uint64_t const sketch_value)
{
  return sketch_value & 1;
}

/*!
 * @brief Checks if the sketch value is a valid position in the graph \a gfa
 *
 * @param[in] gfa graph to check validity.
 * @param[in] sketch_value The sketch value to check.
 *
 * @returns True iff the sketch value is valid.
 */
bool is_sketch_value_valid(GFA const & gfa, uint64_t const sketch_value);

} // namespace weaver
