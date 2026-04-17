#pragma once
/*!
 * @file segment.hpp
 * @brief Contains free functions for segment and vertices.
 */

#include <gfa.h> // gfa_seg_t

namespace weaver
{
class GFA;
class GFALocation;

//! Print debug messages that describe a segment.
/*!
 * Does nothing in release mode.
 *
 * @param[in] segment segment to describe.
 */
void print_segment(gfa_seg_t const & segment);

//! Determines the ID of a vertex from the segment ID going forward (out of the back).
/*!
 * @param[in] segment_id from this segment id.
 *
 * @returns vertex ID that comes out of the back of segment_id.
 *
 * @see reverse_segment_id2vertex_id(uint32_t segment_id)
 */
inline uint32_t forward_segment_id2vertex_id(uint32_t segment_id)
{
  return segment_id << 1;
}

//! Determines the ID of a vertex from the segment ID going reverse (out of the front).
/*!
 * @param[in] segment_id from this segment id.
 *
 * @returns vertex ID that comes out of the front of segment_id.
 *
 * @see forward_segment_id2vertex_id(uint32_t segment_id)
 */
inline uint32_t reverse_segment_id2vertex_id(uint32_t segment_id)
{
  return (segment_id << 1) | 1u;
}

//! Get the sequence prefix of a segment sequence.
/*!
 * If the length of the segment is shorter or equal than \a size then the entire sequence is returned.
 *
 * @param[in] segment segment to get prefix sequence of.
 * @param[in] size maxmimum size of the prefix.
 * @returns The prefix sequence. The type is deduced from the only template parameter.
 * @see get_segment_suffix(gfa_seg_t const & segment, int size)
 *
 * \par Example
 * get_segment_prefix(gfa_seg_t(seq="ACA", len=3), 2) == "AC"<br>
 * get_segment_prefix(gfa_seg_t(seq="ACA", len=3), 3) == "ACA"<br>
 * get_segment_prefix(gfa_seg_t(seq="ACA", len=3), 4) == "ACA"
 */
template <typename Tstring>
inline Tstring get_segment_prefix(gfa_seg_t const & segment, int const size)
{
  if (segment.len > size)
    return Tstring(segment.seq, segment.seq + size); // prefix of size "size"

  return Tstring(segment.seq, segment.len);
}

//! Get the suffix of a sequence that is associated with \a segment.
/*!
 * If the length of the segment is shorter than \a size then the entire sequence is returned.
 *
 * @param[in] segment segment to get suffix sequence of.
 * @param[in] size maxmimum size of the suffix.
 * @returns The suffix sequence. The type is deduced from the only template parameter.
 * @see get_segment_prefix(gfa_seg_t const & segment, int size)
 *
 * \par Example
 * get_segment_suffix(gfa_seg_t(seq="ACT", len=3), 2) == "CT"<br>
 * get_segment_suffix(gfa_seg_t(seq="ACT", len=3), 3) == "ACT"<br>
 * get_segment_suffix(gfa_seg_t(seq="ACT", len=3), 4) == "ACT"
 */
template <typename Tstring>
inline Tstring get_segment_suffix(gfa_seg_t const & segment, int const size)
{
  if (segment.len > size)
    return Tstring(segment.seq + segment.len - size, size);

  return Tstring(segment.seq, segment.len);
}

//! Return the number of bases that are remaining on a segment, given a \a rid, \a pos and \a strand.
/*!
 * @param[in] gfa graph.
 * @param[in] rid segment id.
 * @param[in] pos segment position.
 * @param[in] strand segment strand.
 *
 * @returns The number of bases remaining.
 *
 * @see vertex_bases_remaining(GFA const & gfa, uint64_t value)
 * @see vertex_bases_passed(GFA const & gfa, uint32_t rid, int pos, bool strand)
 * @see vertex_bases_passed(GFA const & gfa, uint64_t value)
 */
int vertex_bases_remaining(GFA const & gfa, uint32_t rid, int pos, bool strand);

//! Return the number of bases that are remaining on a segment, given a vertex position.
/*!
 * @param[in] gfa graph.
 * @param[in] value segment minimizer value.
 *
 * @returns The number of bases remaining.
 *
 * @see vertex_bases_remaining(GFA const & gfa, uint32_t rid, int pos, bool strand)
 * @see vertex_bases_passed(GFA const & gfa, uint32_t rid, int pos, bool strand)
 * @see vertex_bases_passed(GFA const & gfa, uint64_t value)
 */
int vertex_bases_remaining(GFA const & gfa, uint64_t value);

//! Return the number of bases that are remaining on a segment, given a GFA location.
/*!
 * @param[in] gfa graph.
 * @param[in] location the location on the graph.
 *
 * @returns The number of bases remaining.
 *
 * @see vertex_bases_remaining(GFA const & gfa, uint32_t rid, int pos, bool strand)
 */
int vertex_bases_remaining(GFA const & gfa, GFALocation const & location);

//! Return the number of bases that are passed on a segment, given a rid+pos+strand.
/*!
 * @param[in] gfa graph.
 * @param[in] rid segment id.
 * @param[in] pos segment position.
 * @param[in] strand segment strand.
 *
 * @returns The number of bases passed.
 *
 * @see vertex_bases_remaining(GFA const & gfa, uint32_t rid, int pos, bool strand)
 * @see vertex_bases_remaining(GFA const & gfa, uint64_t value)
 * @see vertex_bases_passed(GFA const & gfa, uint64_t value)
 */
int vertex_bases_passed(GFA const & gfa, uint32_t rid, int pos, bool strand);

//! Return the number of bases that are passed on a segment, given a vertex position.
/*!
 * @param[in] gfa graph.
 * @param[in] value segment minimizer value.
 *
 * @returns How many bases are passed on a segment.
 *
 * @see vertex_bases_remaining(GFA const & gfa, uint32_t rid, int pos, bool strand)
 * @see vertex_bases_remaining(GFA const & gfa, uint64_t value)
 * @see vertex_bases_passed(GFA const & gfa, uint32_t rid, int pos, bool strand)
 */
int vertex_bases_passed(GFA const & gfa, uint64_t value);

//! Check if the following position is on a different segment or not present in the graph.
/*!
 * @param[in] gfa graph.
 * @param[in] ref_value graph reference position value.
 *
 * @returns true if the position value is at the end of a segment, false otherwise.
 */
bool is_at_end_of_segment(GFA const & gfa, uint64_t ref_value);

} // namespace weaver
