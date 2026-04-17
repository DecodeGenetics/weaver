#include "sketch_to_string.hpp"

#include <bitset>
#include <cstdint>
#include <sstream>
#include <string>

#include <weaver/constants.hpp>

#include "read_sketch.hpp"
#include "sketch_value.hpp"

namespace weaver
{
//! How much the current minimizer hash should be shifted up in the minimizer key.
int constexpr K_SKETCH_TO_STRING = 23; // assuming k=23
static_assert(K_SKETCH_TO_STRING >= 5);
static_assert(K_SKETCH_TO_STRING <= 31);

std::string sketch_key_to_string(uint64_t const sketch_first)
{
  std::ostringstream ss;

  ss << std::bitset<2 * K_SKETCH_TO_STRING>{sketch_first};

  return ss.str();
}

std::string sketch_value_to_string(uint64_t const sketch_value)
{
  std::ostringstream ss;

  ss << sketch_value_rid(sketch_value)            // rid
     << " " << sketch_value_pos(sketch_value)     // pos
     << " " << sketch_value_strand(sketch_value); // strand

  return ss.str();
}

template <typename u64_pair>
std::string sketch_to_string(u64_pair const & sketch)
{
  std::ostringstream ss;

  ss << std::bitset<2 * K_SKETCH_TO_STRING>{sketch.first} // minimizer hash
     << " " << sketch_value_rid(sketch.second)            // rid
     << " " << sketch_value_pos(sketch.second)            // pos
     << " " << sketch_value_strand(sketch.second);        // strand

  return ss.str();
}

// explicit instantiation
//
template std::string sketch_to_string(std::pair<uint64_t, uint64_t> const & sketch);
template std::string sketch_to_string(ReadSketch const & sketch);

} // namespace weaver
