#pragma once
/*!
 * @file sketch_to_string.hpp
 * @brief Defines the functions to convert sketches to strings.
 */

#include <cstdint>
#include <gfa.h> // gfa_arc_t
#include <string>

namespace weaver
{
//! Returns a pretty string showing the key (bits) of a sketch.
std::string sketch_key_to_string(uint64_t const sketch_first);

//! Returns a pretty string showing the value (bits) of a sketch.
std::string sketch_value_to_string(uint64_t const sketch_value);

//! Returns a pretty string showing both the key (bits) and values of a sketch.
template <typename u64_pair>
std::string sketch_to_string(u64_pair const & sketch);

} // namespace weaver
