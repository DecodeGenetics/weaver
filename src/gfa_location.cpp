/*!
 * @file gfa_location.cpp
 * @brief Implements the GFALocation class
 */

#include "gfa_location.hpp"

#include <iostream>
#include <string>

#include "gfa.hpp"
#include "gfa_arc.hpp"
#include "logging.hpp"
#include "sequence_utils.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp"

namespace weaver
{
GFA const * GFALocation::gfa = nullptr;

GFALocation::GFALocation(uint64_t begin_value) :
  rid(sketch_value_rid(begin_value)), //
  pos(sketch_value_pos(begin_value)), //
  strand(sketch_value_strand(begin_value))
{
}

uint64_t GFALocation::get_value() const
{
  return static_cast<uint64_t>(rid) << 32ull | // rid 32 bits
         static_cast<uint64_t>(pos) << 1ull |  // pos 31 bits
         static_cast<uint64_t>(strand);        // strand 1 bit
}

std::string GFALocation::to_string() const
{
  return sketch_value_to_string(get_value());
}

bool GFALocation::is_valid() const
{
#ifndef NDEBUG
  bool valid{true};

  if (gfa == nullptr)
  {
    print_warning(_HERE_, " gfa is null.");
    valid = false;
  }

  if (static_cast<int>(rid) >= gfa->get_num_segments())
  {
    print_warning(_HERE_, " rid >= segments (", rid, " >= ", gfa->get_num_segments(), ").");
    valid = false;
  }

  if (pos >= gfa->get_segment(rid).len)
  {
    print_warning(_HERE_, " pos >= segment length (", pos, " >= ", gfa->get_segment(rid).len, ").");
    valid = false;
  }

  if (pos < 0)
  {
    print_warning(_HERE_, " negative pos, pos=", pos);
    valid = false;
  }

  return valid;
#else
  return gfa &&                                             // gfa not null
         static_cast<int>(rid) < gfa->get_num_segments() && // rid in range
         pos < gfa->get_segment(rid).len &&                 // pos in range
         pos >= 0 &&                                        // pos non-negative
         (strand == 0 || strand == 1);                      // strand is valid
#endif
}

bool GFALocation::is_same_location(GFALocation const & o) const
{
  return rid == o.rid && pos == o.pos && strand == o.strand;
}

void GFALocation::get_base(std::string & seq) const
{
  assert(is_valid());

  gfa_seg_t const & segment = gfa->get_segment(rid);
  char const base = segment.seq[pos];

  if (is_strand_forward())
    seq += base;
  else
    seq += complement(base);
}

bool GFALocation::is_strand_forward() const
{
  return !strand;
}

void GFALocation::flip_strand()
{
  strand = !strand;
}

void GFALocation::advance_arc(gfa_arc_t const * arc_ptr)
{
  assert(arc_ptr != nullptr);
  assert(static_cast<int>(arc_ptr->v_lv >> 33ull) == rid);
  assert(is_valid());

  gfa_arc_t const & arc = *arc_ptr;
  uint32_t const next_rid = static_cast<uint32_t>(arc.w) >> 1;
  gfa_seg_t const & next_segment = gfa->get_segment(next_rid);

  rid = next_rid;
  strand = (arc.w & 1u);

  if (is_strand_forward())
  {
    pos = 0; // forward strand
  }
  else
  {
    // reverse strand
    int const segment_length = next_segment.len;
    assert(segment_length > 0);
    pos = static_cast<uint32_t>(segment_length) - 1u;
  }

  assert(is_valid());
}

int GFALocation::advance_on_segment(int by)
{
  assert(is_valid());
  assert(by >= 0);
  int safe_advance{by};

  if (is_strand_forward())
  {
    assert(gfa);
    gfa_seg_t const & segment = gfa->get_segment(rid);
    assert(pos < segment.len);
    int const max_advance = segment.len - (pos + 1);

    if (max_advance < safe_advance)
      safe_advance = max_advance;

    pos += safe_advance;
  }
  else
  {
    if (pos < safe_advance)
      safe_advance = pos;

    pos -= safe_advance;
  }

  by -= safe_advance;
  assert(pos < gfa->get_segment(rid).len);
  return by;
}

int GFALocation::advance_on_segment_and_get_sequence(int by, std::string & seq)
{
  assert(gfa);
  assert(pos < gfa->get_segment(rid).len);
  assert(by >= 0);

  int safe_advance{by};
  gfa_seg_t const & segment = gfa->get_segment(rid);

  if (is_strand_forward())
  {
    assert(pos < segment.len);
    int const max_advance = segment.len - (pos + 1);

    if (max_advance < safe_advance)
      safe_advance = max_advance;

    assert(segment.len - pos - safe_advance >= 0);
    seq += std::string(segment.seq + pos, safe_advance);
    pos += safe_advance;
  }
  else
  {
    if (pos < safe_advance)
      safe_advance = pos;

    std::string_view new_seq(segment.seq + pos - safe_advance + 1, safe_advance);
    seq += get_reverse_complement(new_seq);
    pos -= safe_advance;
  }

  by -= safe_advance;
  assert(pos < gfa->get_segment(rid).len);
  return by;
}

int GFALocation::advance_until_and_get_sequence(GFALocation const & end,
                                                std::vector<gfa_arc_t const *> const & arcs,
                                                std::string & seq)
{
  assert(is_valid());
  assert(gfa);
  assert(pos < gfa->get_segment(rid).len);

  if (seq.size() > 10000)
    print_warning(_HERE_, " Unexpectedly large sequence has been extracted: ", seq.size());

  for (gfa_arc_t const * arc_ptr : arcs)
  {
    int constexpr by = std::numeric_limits<int>::max();
    advance_on_segment_and_get_sequence(by, seq);
    get_base(seq);
    gfa_arc_t const & arc = *(arc_ptr);

    /*
    print_debug(_HERE_,
                " rid=",
                rid,
                " strand=",
                static_cast<uint32_t>(strand),
                " ",
                (arc.v_lv >> 33ull),
                "|",
                (arc.v_lv >> 32ull) & 1,
                " -> ",
                arc.w >> 1,
                "|",
                arc.w & 1);
    */

    assert(static_cast<int>(arc.v_lv >> 33ull) == rid);
    assert(((arc.v_lv >> 32ull) & 1) == strand);

    rid = (arc.w >> 1);
    strand = (arc.w & 1);

    if (is_strand_forward())
    {
      pos = 0;
    }
    else
    {
      int const segment_length = gfa->get_segment(rid).len;
      assert(segment_length > 0);
      pos = segment_length - 1;
    }

    assert(is_valid());
  }

  assert(rid == end.rid);
  assert(strand == end.strand);
  int diff = end.is_strand_forward() ? end.pos - pos : pos - end.pos;
  assert(diff >= 0);
  int by = advance_on_segment_and_get_sequence(diff, seq);
  assert(static_cast<int>(pos) < gfa->get_segment(rid).len);
  get_base(seq);
  return by;
}

int GFALocation::advance_when_same_contig(int by, std::vector<gfa_arc_t const *> & new_arcs)
{
  assert(is_valid());
  assert(by >= 0);

  if (by == 0)
    return 0;

  by = advance_on_segment(by);
  assert(is_valid());

  while (by > 0)
  {
    // Try to go onto a new segment
    gfa_seg_t const & segment = gfa->get_segment(rid);
    uint32_t const vertex_id{rid << 1 | static_cast<uint32_t>(strand)};
    bool is_advanced{false};
    assert(gfa != nullptr);
    gfa_arc_t const * arc_end_ptr = arc_end(*gfa, vertex_id);

    for (gfa_arc_t const * arc_ptr = arc_begin(*gfa, vertex_id); arc_ptr != arc_end_ptr; ++arc_ptr)
    {
      gfa_arc_t const & arc = *arc_ptr;

      if (arc.rank != segment.rank)
        continue;

      uint32_t const next_rid = static_cast<uint32_t>(arc.w) >> 1;
      gfa_seg_t const & next_segment = gfa->get_segment(next_rid);

      if (next_segment.rank != segment.rank || next_segment.snid != segment.snid)
        continue;

      // same contig, advance
      --by; // the last base of the previous segment
      new_arcs.push_back(arc_ptr);
      is_advanced = true;
      assert(is_valid());
      rid = next_rid;
      strand = (arc.w & 1u);

      if (is_strand_forward())
      {
        pos = 0; // forward strand
      }
      else
      {
        int const segment_length = next_segment.len;
        assert(segment_length > 0);
        pos = static_cast<uint32_t>(segment_length) - 1u;
      }

      by = advance_on_segment(by);
      break; // ok because it is always an unique arc that advances on the same
             // contig
    }

    if (!is_advanced)
      break;
  }

  assert(is_valid());
  return by;
}

int GFALocation::advance_when_same_contig_and_get_sequence(int by,
                                                           std::vector<gfa_arc_t const *> & new_arcs,
                                                           std::string & seq)
{
  if (by == 0)
    return 0;

#ifndef NDEBUG
  int const old_by_seq = by + static_cast<int>(seq.size()); // for debugging
#endif                                                      // NDEBUG
  assert(is_valid());
  assert(by > 0);

  //--by; // reduce by now, and then add a base in the end without changing by
  by = advance_on_segment_and_get_sequence(by, seq);
  assert(by + static_cast<int>(seq.size()) == old_by_seq);
  assert(is_valid());

  while (by > 0)
  {
    // Try to go onto a new segment
    gfa_seg_t const & segment = gfa->get_segment(rid);
    uint32_t const vertex_id{static_cast<uint32_t>(rid) << 1 | strand};
    bool is_advanced{false};
    assert(gfa != nullptr);
    gfa_arc_t const * arc_end_ptr = arc_end(*gfa, vertex_id);

    for (gfa_arc_t const * arc_ptr = arc_begin(*gfa, vertex_id); arc_ptr != arc_end_ptr; ++arc_ptr)
    {
      gfa_arc_t const & arc = *arc_ptr;

      if (arc.rank != segment.rank)
        continue;

      uint32_t const next_rid = static_cast<uint32_t>(arc.w) >> 1;
      gfa_seg_t const & next_segment = gfa->get_segment(next_rid);

      if (next_segment.rank != segment.rank || next_segment.snid != segment.snid)
        continue;

      // same contig, advance
      get_base(seq);
      --by;
      assert(by + static_cast<int>(seq.size()) == old_by_seq);

      new_arcs.push_back(arc_ptr);
      is_advanced = true;
      assert(is_valid());
      rid = next_rid;
      strand = (arc.w & 1u);

      if (is_strand_forward())
      {
        pos = 0; // forward strand
      }
      else
      {
        int const segment_length = next_segment.len;
        assert(segment_length > 0);
        pos = static_cast<uint32_t>(segment_length) - 1u;
      }

      by = advance_on_segment_and_get_sequence(by, seq);
      assert(by + static_cast<int>(seq.size()) == old_by_seq);
      break; // break is ok because there is always an unique arc that advances on the same contig
    }

    // below condition is met if no arc is on the same contig
    if (!is_advanced)
      break;
  }

  assert(is_valid());
  assert(by + static_cast<int>(seq.size()) == old_by_seq);
  return by;
}

} // namespace weaver
