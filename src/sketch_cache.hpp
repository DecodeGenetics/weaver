#pragma once
/*!
 * @file sketch_cache.hpp
 * @brief Defines the SketchCache class.
 */

#include <array>
#include <cstdint>
#include <limits>
#include <utility>

#include <weaver/constants.hpp>

namespace weaver
{
//! Caches containing sketches.
/*!
 * Sketches in cache are defined as follows:
 *
 * (64 bits) Key:
 * minimizer hash: 46 bits,
 * partial hash of next kmers: 18 bits, first four bits of the next hash are skipped (it is so often 0000)
 *
 * (64 bits) Value:
 * segment ID: 32 bits
 * segment offset: 31 bits
 * strand: 1 bit
 *
 * The SketchCache stores the internal state while sketching. By storing everything in the cache we can copy the cache
 * when needed, i.e. when the graph branches.
 */
template <typename u64_pair>
class SketchCache
{
public:
  // using u64_pair = std::pair<uint64_t, uint64_t>;

  //! current kmer. First is forward, second is reverse.
  u64_pair kmer{0, 0};
  uint64_t last_key{std::numeric_limits<uint64_t>::max()};                                  //!< last written key
  u64_pair min{std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max()}; //!< current minimizer
  std::array<u64_pair, MAX_W> buf;                                                          //!< kmer buffer.

  //! Previous minimizer reference ID. -1 if no contig was previous to this one.
  int prev_rid{-1};

  //! Size of previous minimizer reference ID contig. The value is negative iff that was from read in reverse.
  int prev_rid_size{0};

  //! Stores the length of the \a kmer.
  int l{0};

  //! The current index in the \a buf array. Updates after reading a base the string to sketch.
  int buf_i{0};

  //! min_i is the index of minimizer in the \a buf array.
  int min_i{0};

  SketchCache();                                              //!< Construct an empty cache.
  SketchCache(SketchCache const &) = default;                 //!< Default copy constructor.
  SketchCache(SketchCache &&) noexcept = default;             //!< Default move constructor.
  SketchCache & operator=(SketchCache const &) = default;     //!< Default copy assignment.
  SketchCache & operator=(SketchCache &&) noexcept = default; //!< Default move assignment.
  ~SketchCache() = default;

  //! True iff minimizer was read on the forward strand.
  inline bool was_prev_forward()
  {
    return prev_rid_size >= 0;
  }

  //! Get the size of the minimizer reference ID contig.
  inline int get_prev_str_size()
  {
    return was_prev_forward() ? prev_rid_size : -prev_rid_size;
  }

  void clear();

  //! Clears the kmer buffer \a buf.
  void clear_buffer(int buf_len);
};
} // namespace weaver
