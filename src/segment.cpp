/*!
 * @file segment.cpp
 * @brief Implements the functions to handle for segments and vertices.
 */
#include "segment.hpp"

#include <cstdint> // uint32_t uint64_t
#include <gfa.h>   // gfa_seg_t

#include "gfa.hpp"
#include "gfa_location.hpp"
#include "logging.hpp"
#include "sketch_value.hpp"

namespace weaver
{
void print_segment(gfa_seg_t const & segment)
{
#ifndef NDEBUG
  if (segment.len > 50)
  {
    char const * end_of_str = segment.seq + segment.len - 20;
    std::string end_str(end_of_str, 20);
    print_info(segment.name, '\t', segment.len, '\t', std::string_view(segment.seq, 20), "...", end_str);
  }
  else
  {
    print_info(segment.name, '\t', segment.len, '\t', segment.seq);
  }
#else
  print_info(segment.name, '\t', segment.len);
#endif // NDEBUG
}

int vertex_bases_remaining(GFA const & gfa, uint32_t rid, int pos, bool strand)
{
#ifndef NDEBUG
  {
    gfa_seg_t const & segment = gfa.get_segment(rid);
    assert(pos < segment.len);
  }
#endif // NDEBUG

  if (strand == 1)
    return pos;

  gfa_seg_t const & segment = gfa.get_segment(rid);
  return segment.len - 1 - pos;
}

int vertex_bases_remaining(GFA const & gfa, uint64_t value)
{
  return vertex_bases_remaining(gfa, sketch_value_rid(value), sketch_value_pos(value), sketch_value_strand(value));
}

int vertex_bases_remaining(GFA const & gfa, GFALocation const & location)
{
  return vertex_bases_remaining(gfa, location.rid, location.pos, location.strand);
}

int vertex_bases_passed(GFA const & gfa, uint32_t rid, int pos, bool strand)
{
  return vertex_bases_remaining(gfa, rid, pos, !strand);
}

int vertex_bases_passed(GFA const & gfa, uint64_t value)
{
  return vertex_bases_remaining(gfa, sketch_value_rid(value), sketch_value_pos(value), !sketch_value_strand(value));
}

int vertex_bases_passed(GFA const & gfa, GFALocation const & location)
{
  return vertex_bases_remaining(gfa, location.rid, location.pos, !location.strand);
}

bool is_at_end_of_segment(GFA const & gfa, uint64_t ref_value)
{
  int const ref_pos = sketch_value_pos(ref_value);
  bool const ref_strand = sketch_value_strand(ref_value);

  if (ref_strand)
    return ref_pos == 0; // reverse strand

  return (ref_pos + 1) == gfa.get_segment(sketch_value_rid(ref_value)).len; // forward strand
}

} // namespace weaver
