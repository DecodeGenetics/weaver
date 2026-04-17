#pragma once
/*!
 * @file sketch.hpp
 * @brief Defines the functions to sketch the minimizers of a sequence.
 *
 * @details
 * The most important function here is simple called \a sketch() which pushes back minimizer sketches of a sequence to a
 * vector.
 */

#include <string_view> // std::string_view
#include <utility>     // std::pair
#include <vector>      // std::vector

#include "sketch_cache.hpp" // weaver::SketchCache

namespace weaver
{
/*!
 * @brief Make minimizer sketches of a string.
 *
 * @details
 * Sketches are the kmers with minimum local hash in a window of size w (minimizers) and their locations on a segment.
 * This function is inspired by the \a mm_sketch function in minimap2.
 *
 * @param[out] sketches key-value sketches of a string. Key is 46 bits kmer, 18 bits least significant hash of the next
 *                      minimizer. The value contains information about where the minimizer is located at. It is split
 *                      into rid (32 bits), center position of the minimizer on the rid (31 bits) and the minimizer
 *                      strand/orienation (1 bit). Forward strand is 0 and the reverse strand is 1.
 * @param[in,out] cache Cache of sketches.
 * @param[in] str Input string to sketch.
 * @param[in] w Window size.
 * @param[in] k kmer size.
 * @param[in] rid Reference/contig ID. Stored in value.
 * @param[in] is_reading_forward False iff reading should start at end. Useful for partial sketching.
 * @param[in] max_read_size Indicates how many bases should be read in \a str. Useful for partial sketching.
 *                          Use \a max_read_size < 0 if \a str should all be read.
 * @param[in] read_init First position to read from str. Default 0 (which is first position).
 * @param[in] add_pos Add some number to the positions if reading forward, or subtract if reading in reverse.
 *                    Default is 0.
 */
template <typename u64_pair>
void sketch(std::vector<u64_pair> & sketches,
            SketchCache<u64_pair> & cache,
            std::string_view str,
            int const w,
            int const k,
            int const rid,
            bool const is_reading_forward,
            int const max_read_size = -1,
            int const read_init = 0,
            int const add_pos = 0);

/*!
 * @brief Same as sketch() except sketches only each boundary of the input \a str .
 *
 * @details
 * Useful if the input string has been already sketched with other sequences surrounding it so the center is not
 * adding any new sketches. If the input string is small enough it might be fully sketched.
 */
template <typename u64_pair>
void sketch_boundary(std::vector<u64_pair> & sketches,
                     SketchCache<u64_pair> & cache,
                     std::string_view str,
                     int const w,
                     int const k,
                     int const rid,
                     bool const is_reading_forward = true,
                     int const max_read_size = -1,
                     int const read_init = 0,
                     int const add_pos = 0);

/*!
 * @brief Pushes the last sketch to \a sketches from \a cache.
 *
 * @details
 * When a string is sketched the current minimizer is stored in the cache. After sketching,
 * it is therefore needed to push the last sketch using this function when the string/segment has no continuation.
 *
 * @param[out]    sketches Key-value sketches of a string. Key is 46 bits kmer, and the remaining 18 bits
 *                         are the most significant hash of the next minimizer.
 * @param[in,out] cache Cache sketch data.
 * @param[in]     max_read_size How many bases should be read in the input string.
 * @param[in]     k kmer size.
 */
template <typename u64_pair>
void push_last_sketch(std::vector<u64_pair> & sketches, SketchCache<u64_pair> & cache, int max_read_size, int const k);

/*!
 * @brief Get the position of the first sketch of an input \a str .
 *
 * @details
 * If the sequence is too short such that no minimizer was sketched, a default value is returned (determined from
 * k and w), but is always larger than the input string and indicates the highest possible position of the first
 * sketch.
 *
 * @param[in] str String to sketch.
 * @param[in] rid Reference/contig ID.
 * @param[in] w Window size.
 * @param[in] k kmer size.
 * @param[in] is_reading_forward False iff reading should start at end. Useful for partial sketching.
 *
 * @returns The position of the first sketch.
 */
int get_position_of_the_first_sketch(std::string_view str, //
                                     int const rid,
                                     int const w,
                                     int const k,
                                     bool is_reading_forward);

} // namespace weaver
