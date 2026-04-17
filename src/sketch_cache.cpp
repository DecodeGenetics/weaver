/*!
 * @file sketch_cache.cpp
 * @brief Implementation of the SketchCache class.
 */

#include "sketch_cache.hpp"

#include <algorithm> // std::fill
#include <cassert>   // assert
#include <cstdint>   // uint64_t
#include <limits>    // std::numeric_limits
#include <utility>   // std::pair

#include "read_sketch.hpp"

namespace weaver
{
template <typename u64_pair>
SketchCache<u64_pair>::SketchCache()
{
  clear_buffer(MAX_W); // Clear all buffer
}

template <typename u64_pair>
void SketchCache<u64_pair>::clear()
{
  kmer = {0, 0};
  last_key = std::numeric_limits<uint64_t>::max();
  min = {std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max()};
  clear_buffer(MAX_W);
  // prev_rid, prev_rid_size unchanged
  l = 0;
  buf_i = 0;
  min_i = 0;
}

template <typename u64_pair>
void SketchCache<u64_pair>::clear_buffer(int buf_len)
{
  assert(buf_len >= 0);
  assert(buf_len <= static_cast<int>(MAX_W));

  u64_pair const empty(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max());
  std::fill(buf.begin(), buf.begin() + buf_len, empty);
}

// explicit instantiation
//
template SketchCache<std::pair<uint64_t, uint64_t>>::SketchCache();
template void SketchCache<std::pair<uint64_t, uint64_t>>::clear();
template void SketchCache<std::pair<uint64_t, uint64_t>>::clear_buffer(int buf_len);

template SketchCache<ReadSketch>::SketchCache();
template void SketchCache<ReadSketch>::clear();
template void SketchCache<ReadSketch>::clear_buffer(int buf_len);

} // namespace weaver
