/*!
 * @file sketch.cpp
 * @brief Implements the functions to sketch a sequence.
 */

#include "sketch.hpp"

#include <bitset>
#include <cassert>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "logging.hpp"
#include "read_sketch.hpp"
#include "sequence_utils.hpp"
#include "sketch_cache.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp"

namespace
{
/*
inline uint64_t hash64(uint64_t key, uint64_t mask)
{
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ (key >> 14);
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ (key >> 28);
  key = (key + (key << 31)) & mask;
  return key;
}
*/

// Alternative hashing
inline uint64_t hash64_alt(uint64_t key, uint64_t mask)
{
  key = (~key + (key << 21)); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = ((key + (key << 3)) + (key << 8)); // key * 265
  key = key ^ (key >> 14);
  key = ((key + (key << 2)) + (key << 4)); // key * 21
  key = key ^ (key >> 28);
  key = (key + (key << 31));
  return key & mask;
}

} // namespace

namespace weaver
{
template <typename u64_pair>
void sketch_boundary(std::vector<u64_pair> & sketches,
                     SketchCache<u64_pair> & cache,
                     std::string_view str,
                     int const w,
                     int const k,
                     int const rid,
                     bool is_reading_forward,
                     int const max_read_size,
                     int const read_init,
                     int const add_pos)
{
  int const boundary_size = k + w - 1;

  if (static_cast<int>(str.size()) <= (3 * boundary_size) || !is_reading_forward || max_read_size >= 0 || add_pos > 0)
  {
    // string is small enough, or a special case where boundary sketching has not been implemented is requested
    sketch(sketches, cache, str, w, k, rid, is_reading_forward, max_read_size, read_init, add_pos);
    return;
  }

  int const pos_offset = static_cast<int>(str.size()) - boundary_size;
  assert(pos_offset > 0);
  assert(is_reading_forward); // reverse not implemented yet!
  assert(max_read_size < 0);  // boundary sketching does not make sense otherwise
  assert(add_pos == 0);

  // print_warning("str=", str, " (size=", str.size(), ")");
  // print_warning("prefix_str=", prefix_str);
  // print_warning("suffix_str=", suffix_str);
  // print_warning("pos_offset=", pos_offset);
  sketch(sketches, cache, str, w, k, rid, is_reading_forward, /*max_read_size=*/boundary_size, read_init, add_pos);
  push_last_sketch(sketches, cache, /*max_read_size=*/-1, k);
  cache.clear();
  sketch(sketches, cache, str, w, k, rid, is_reading_forward, /*max_read_size=*/-1, read_init + pos_offset, add_pos);
}

template <typename u64_pair>
void sketch(std::vector<u64_pair> & sketches,
            SketchCache<u64_pair> & cache,
            std::string_view str,
            int const w,
            int const k,
            int const rid,
            bool is_reading_forward,
            int const max_read_size,
            int const read_init,
            int const add_pos)
{
  uint64_t const rev_shift = 2ull * static_cast<uint64_t>(k - 1); // how much reverse kmer should be shifted
  uint64_t const mask = (1ull << (2 * k)) - 1ull;                 // kmer mask
  int const str_size{static_cast<int>(str.size())};               // size of the input string
  int str_read_size{0};                                           // how far into the input str should be read

  if (max_read_size >= 0)
  {
    if (is_reading_forward)
      str_read_size = max_read_size - k / 2;
    else
      str_read_size = str_size - max_read_size - 1 + k / 2;
  }

  // loop over str
  for (int i{read_init}; i < str_size; ++i)
  {
    int str_i;  // index in str to use
    int idx_i;  // index to put in minimizer index
    uint64_t c; // 0=a,1=c,2=g,3=t,4=n (or other)

    if (is_reading_forward)
    {
      str_i = i;
      c = nt5_char_to_ull(str[str_i]);
      idx_i = str_i - k / 2 + add_pos;
    }
    else
    {
      str_i = str_size - i - 1;
      c = complement_nt5_char_to_ull(str[str_i]);
      idx_i = str_i + k / 2 - add_pos;
    }

    auto is_stopping = [&](u64_pair const & min) -> bool
    {
      if (max_read_size < 0)
        return false;

      if (min.second == std::numeric_limits<uint64_t>::max())
        return false; // no minimizer saved

      // TODO test if it is worth it to use
      // if ((is_reading_forward && str_i <= str_read_size) || (str_i >= str_read_size))
      //  return false;

      if (sketch_value_rid(min.second) != rid)
        return false; // different rids

      int const sketch_pos = sketch_value_pos(min.second);

      if (is_reading_forward)
        return sketch_pos > str_read_size;
      else
        return sketch_pos < str_read_size;
    };

    if (c > 3)
    {
      cache.l = 0;           // reset kmer length if we see an ambigous base
      cache.clear_buffer(w); // clear first w elements in buffer
      continue;
    }

    cache.kmer.first = ((cache.kmer.first << 2ull) | c) & mask;                  // forward kmer
    cache.kmer.second = (cache.kmer.second >> 2ull) | ((3ull ^ c) << rev_shift); // reverse kmer

    // skip kmers with same forward and reverse complement sequence
    // the problem is that we cannot pick which one to use
    if (cache.kmer.first == cache.kmer.second)
      continue;

    ++cache.l;
    u64_pair current_kv{std::numeric_limits<uint64_t>::max(), // current
                                                              // key-value pair
                        std::numeric_limits<uint64_t>::max()};

    if (cache.l >= k)
    {
      bool const is_reverse_smaller = cache.kmer.second < cache.kmer.first;

      if (is_reverse_smaller)
        current_kv.first = hash64_alt(cache.kmer.second, mask);
      else
        current_kv.first = hash64_alt(cache.kmer.first, mask);

      if ((is_reading_forward && idx_i >= 0) || (!is_reading_forward && idx_i < str_size))
      {
        // zero or positive value of idx_i
        assert(rid >= 0);
        current_kv.second = static_cast<uint64_t>(rid) << 32 |                             // rid 32 bits
                            static_cast<uint32_t>(idx_i) << 1 |                            // rid index 31 bits
                            static_cast<uint8_t>(is_reverse_smaller ^ is_reading_forward); // strand 1 bit
      }
      else
      {
        // when idx_i is negative, we should use the previous rid rather than the current one
        assert(cache.prev_rid >= -1); // another RID should have been read
        int const rem_bases = is_reading_forward ? -idx_i - 1 : idx_i - 1 - (str_size - 1);

        // print_debug(_HERE_,
        //             " ",
        //             is_reading_forward,
        //             " rem_bases=",
        //             rem_bases,
        //             " idx_i=",
        //             idx_i,
        //             " str_size=",
        //             str_size);

        assert(rem_bases >= 0);
        int min_rid_pos;

        if (cache.prev_rid_size >= 0)
        {
          min_rid_pos = cache.prev_rid_size - 1 - rem_bases;

          if (min_rid_pos < 0)
          {
            print_debug(_HERE_, " skipping. cache.min_rid_size=", cache.prev_rid_size, " rem_bases=", rem_bases);
            continue; // we cannot use this
          }
        }
        else
        {
          min_rid_pos = rem_bases;

          if (rem_bases > -cache.prev_rid_size)
          {
            print_debug(_HERE_, " skipping. cache.prev_rid_size=", cache.prev_rid_size, " rem_bases=", rem_bases);
            continue;
          }
        }

        current_kv.second = static_cast<uint64_t>(cache.prev_rid) << 32 |                        // rid 32 bits
                            static_cast<uint32_t>(min_rid_pos) << 1 |                            // rid index 31 bits
                            static_cast<uint8_t>(is_reverse_smaller ^ cache.was_prev_forward()); // strand 1 bit
      }
    } // ENDs if (cache.l >= k)

    cache.buf[cache.buf_i] = current_kv;

    // special case for the first window, it is needed because identical k-mers are not stored yet
    if (cache.l == w + k - 1 && cache.min.first != std::numeric_limits<uint64_t>::max())
    {
      for (int j{cache.buf_i + 1}; j < w; ++j)
      {
        if (cache.min.first == cache.buf[j].first && cache.buf[j].second != cache.min.second)
          sketches.push_back(cache.buf[j]);
      }

      for (int j{0}; j < (cache.buf_i - 1); ++j)
      {
        if (cache.min.first == cache.buf[j].first && cache.buf[j].second != cache.min.second)
          sketches.push_back(cache.buf[j]);
      }
    }

    if (current_kv.first <= cache.min.first)
    {
      // a new minimum. Write the old min if there is one
      if (cache.l >= w + k - 1 && cache.min.first != std::numeric_limits<uint64_t>::max())
        sketches.push_back(cache.min);

      cache.min = current_kv;    // set new min
      cache.min_i = cache.buf_i; // set new index to min
    }
    else if (cache.buf_i == cache.min_i) // old min has moved outside the window
    {
      // Push minimizer
      if (cache.l >= w + k - 1 && cache.min.first != std::numeric_limits<uint64_t>::max())
        sketches.push_back(cache.min);

      // Find new minimizer in buffer. cache.min_i will be set as the index of the new minimizer
      {
        uint64_t new_min_key{std::numeric_limits<uint64_t>::max()}; // hash of new minimizer

        // Check values in buf_i+1..w-1
        for (int j{cache.buf_i + 1}; j < w; ++j)
        {
          if (cache.buf[j].first <= new_min_key)
          {
            new_min_key = cache.buf[j].first;
            cache.min_i = j;
          }
        }

        // Check values in 0..buf_i. The two loops are necessary when there are identical k-mers.
        for (int j{0}; j <= cache.buf_i; ++j)
        {
          if (cache.buf[j].first <= new_min_key)
          {
            new_min_key = cache.buf[j].first;
            cache.min_i = j;
          }
        }

        // Set new minimizer in cache
        cache.min = cache.buf[cache.min_i];
      }

      // write identical k-mers in the new window
      if (cache.l >= w + k - 1 && cache.min.first != std::numeric_limits<uint64_t>::max())
      {
        // we need to have these two loops to make sure the output is sorted
        for (int j{cache.buf_i + 1}; j < w; ++j)
        {
          if (cache.min.first == cache.buf[j].first && cache.min.second != cache.buf[j].second)
            sketches.push_back(cache.buf[j]);
        }

        for (int j{0}; j <= cache.buf_i; ++j)
        {
          if (cache.min.first == cache.buf[j].first && cache.min.second != cache.buf[j].second)
            sketches.push_back(cache.buf[j]);
        }
      }
    }

    // Update buf_i for next iteration
    if (++cache.buf_i == w)
      cache.buf_i = 0;

    if (is_stopping(current_kv))
    {
      // print_debug(_HERE_,
      //             " Stopping at ",
      //             sketch_to_string(current_kv),
      //             " rid=",
      //             rid,
      //             " str_read_size=",
      //             str_read_size);
      break;
    }
  }

  // Set previous RID on cache
  cache.prev_rid = rid;
  cache.prev_rid_size = is_reading_forward ? str_size : -str_size;
}

template <typename u64_pair>
void push_last_sketch(std::vector<u64_pair> & sketches, SketchCache<u64_pair> & cache, int max_read_size, int const k)
{
  if (cache.min.first != std::numeric_limits<uint64_t>::max()) // min not missing
  {
    assert(cache.prev_rid >= 0); // some rid has been read
    bool const is_reading_forward = cache.was_prev_forward();

    // how far into the input str should be read
    int str_read_size;

    if (max_read_size >= 0)
    {
      if (is_reading_forward)
        str_read_size = max_read_size - k / 2;
      else
        str_read_size = cache.get_prev_str_size() - max_read_size - 1 + k / 2;
    }
    else
    {
      if (is_reading_forward)
        str_read_size = std::numeric_limits<int>::max(); // set maximum limit
      else
        str_read_size = 0;
    }

    int const sketch_pos = sketch_value_pos(cache.min.second);

    if (sketch_value_rid(cache.min.second) != cache.prev_rid || // different rid
        (is_reading_forward && sketch_pos <= str_read_size) ||  // forward pos ok
        (!is_reading_forward && sketch_pos >= str_read_size))   // reverse pos ok
    {
      // print_debug(_HERE_, " sketch_pos=", sketch_pos, " str_read_size=",
      // str_read_size);
      sketches.push_back(cache.min);
    }
  }
}

int get_position_of_the_first_sketch(std::string_view str, //
                                     int const rid,        //
                                     int const w,
                                     int const k,
                                     bool is_reading_forward)
{
  int max_read_size{k + k / 2 + 1 + 2 * w};

  // Sketch next segment and find the first sketch with an empty cache
  SketchCache<std::pair<uint64_t, uint64_t>> empty_cache;
  std::vector<std::pair<uint64_t, uint64_t>> test_sketches;

  sketch(/*sketches=*/test_sketches,
         /*cache=*/empty_cache,
         /*str=*/str,
         /*w=*/w,
         /*k=*/k,
         /*rid=*/rid,
         /*is_reading_forward=*/is_reading_forward,
         /*max_read_size=*/max_read_size,
         /*read_init=*/0,
         /*add_pos=*/0);

  // test_sketches.size() == 0 can happen if sequence is short or contains ambigous bases
  if (test_sketches.size() == 0 && empty_cache.min.first != std::numeric_limits<uint64_t>::max())
    push_last_sketch(test_sketches, empty_cache, max_read_size, k);

  if (test_sketches.size() > 0)
  {
    if (is_reading_forward)
      max_read_size = sketch_value_pos(test_sketches[0].second) - 1 + k / 2;
    else
      max_read_size = empty_cache.get_prev_str_size() - sketch_value_pos(test_sketches[0].second) + k / 2;
  }
  else
  {
    print_debug(_HERE_, " W no key for str: ", str);
  }

  return max_read_size;
}

// explicit instantiation
//
template void sketch(std::vector<std::pair<uint64_t, uint64_t>> & sketches,
                     SketchCache<std::pair<uint64_t, uint64_t>> & cache,
                     std::string_view str,
                     int const w,
                     int const k,
                     int const rid,
                     bool const is_reading_forward,
                     int const max_read_size = -1,
                     int const read_init = 0,
                     int const add_pos = 0);

template void sketch_boundary(std::vector<std::pair<uint64_t, uint64_t>> & sketches,
                              SketchCache<std::pair<uint64_t, uint64_t>> & cache,
                              std::string_view str,
                              int const w,
                              int const k,
                              int const rid,
                              bool const is_reading_forward,
                              int const max_read_size = -1,
                              int const read_init = 0,
                              int const add_pos = 0);

template void push_last_sketch(std::vector<std::pair<uint64_t, uint64_t>> & sketches,
                               SketchCache<std::pair<uint64_t, uint64_t>> & cache,
                               int max_read_size,
                               int const k);

template void sketch(std::vector<ReadSketch> & sketches,
                     SketchCache<ReadSketch> & cache,
                     std::string_view str,
                     int const w,
                     int const k,
                     int const rid,
                     bool const is_reading_forward,
                     int const max_read_size = -1,
                     int const read_init = 0,
                     int const add_pos = 0);

template void push_last_sketch(std::vector<ReadSketch> & sketches,
                               SketchCache<ReadSketch> & cache,
                               int max_read_size,
                               int const k);

} // namespace weaver
