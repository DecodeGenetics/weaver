/*!
 * @file icu.cpp
 * @brief Implements functions to create an "I see you" index.
 */
#include "icu.hpp"

#include <cstdint>
#include <gfa.h>

#include "constants.hpp"
#include "gfa.hpp"
#include "gfa_arc.hpp"
#include "hashmap.hpp"
#include "logging.hpp" // print_debug
#include "segment.hpp"
#include "sketch_value.hpp" // sketch_value_rid, sketch_value_pos

namespace
{
//! Recursively add values to the ICU index.
void icu_recursive(weaver::T_icu & icu,      // icu index being created
                   weaver::GFA const & gfa,  // graph
                   uint64_t start_vertex_id, // set at the start of the recursion
                   uint64_t vertex_id,       // current vertex id
                   int distance,             // current distance covered, starts at 1
                   int const max_distance)   // maximum distance to cover
{
  using namespace weaver;

  if (distance > max_distance)
    return;

  gfa_arc_t const * arc_end_ptr = arc_end(gfa, vertex_id);

  for (gfa_arc_t const * arc_ptr = arc_begin(gfa, vertex_id); arc_ptr != arc_end_ptr; ++arc_ptr)
  {
    gfa_arc_t const & arc = *arc_ptr;

    {
      uint64_t const key = (start_vertex_id << 32) | arc.w;
      T_icu::iterator find_it = icu.find(key);

      if (find_it == icu.end())
        icu[key] = distance;
      else if (distance < find_it->second)
        find_it->second = distance; // set to the smallest value in case of multiple entries
    }

    uint64_t const next_rid = static_cast<uint32_t>(arc.w) >> 1; // rid of segment on the other side of the arc
    gfa_seg_t const & next_segment = gfa.get_segment(next_rid);
    assert(next_segment.len > 0);

    // print_debug(_HERE_, " searching with icu on arc ", arc_to_string(*arc_ptr));
    icu_recursive(icu, gfa, start_vertex_id, arc.w, distance + next_segment.len, max_distance);
  }
}

} // namespace

namespace weaver
{
bool is_exact_distance(GFA const & gfa,   // graph
                       T_icu const & icu, // The "I see you" index
                       uint64_t from,
                       uint64_t to,
                       std::vector<gfa_arc_t const *> & new_arcs,
                       int distance)
{
  // Negative distance not allowed
  if (distance < 0)
    return false;

  int const from_rid = sketch_value_rid(from);
  int const from_pos = sketch_value_pos(from);
  bool const from_strand = sketch_value_strand(from);

  int const to_rid = sketch_value_rid(to);
  int const to_pos = sketch_value_pos(to);
  bool const to_strand = sketch_value_strand(to);

  // trivial case if RIDs are the same
  if (from_rid == to_rid)
  {
    if (from_strand != to_strand)
    {
      print_debug(_HERE_, " different strands.");
      return false; // different strands
    }

    if (from_strand == 0)
      return to_pos - from_pos == distance; // forward
    else
      return from_pos - to_pos == distance; // reverse
  }

  // reduce distance by the amount of the "to" segment left
  if (to_strand == 0)
  {
    // "to" is going forward, so reduce distance by its position
    to -= to_pos << 1;
    distance -= to_pos;
    /* NOTE: to_pos should be changed here as well, however it is never used again. */
    // to_pos = 0;
  }
  else
  {
    // "to" is moving in reverse, so extend it to its end
    gfa_seg_t const & to_segment = gfa.get_segment(to_rid);
    int to_segment_remaining = to_segment.len - 1 - to_pos;

    assert(to_segment_remaining >= 0);

    to += to_segment_remaining << 1;
    distance -= to_segment_remaining;
    // NOTE: to_pos should be changed here as well, however it is never used again
    // to_pos = static_cast<int>(to_segment.len) - 1l;
  }

  // check distance again after adjusting it
  if (distance < 0)
    return false;

  /* NOTE: it is tempting to just check for exact distance using the icu index but
   * there are at least two problems:
   *   1) the index only holds the shortest distance.
   *   2) the index does not contain the arcs crossed/path.
   * if the index is extended to keep also this information it might grow in
   * size very fast which would slow down
   * everything.
   *
   * TLDR; Only use icu to check if the two position can (under any
   * circumstance) have this exact distance between them.
   */
  {
    // make sure that "from" can see "to" using the icu index.
    uint64_t const icu_key = static_cast<uint64_t>(from_rid) << 33 |    // from
                             static_cast<uint64_t>(from_strand) << 32 | //
                             static_cast<uint64_t>(to_rid) << 1 |       // to
                             static_cast<uint64_t>(to_strand);          //

    auto find_it = icu.find(icu_key);

    if (find_it == icu.end() || find_it->second > distance)
      return false;
  }

  // go forward onto a new segment.
  gfa_seg_t const & from_segment = gfa.get_segment(from_rid);

  if (from_strand == 0)
  {
    distance -= (from_segment.len - from_pos); // go forward (out of the back)
  }
  else
  {
    distance -= (from_pos + 1); // go backward
  }

  // Negative distance indicates that the distance must be further
  if (distance < 0)
    return false;

  uint32_t const from_vertex{static_cast<uint32_t>(from_rid) << 1 | from_strand};
  gfa_arc_t const * arc_end_ptr = arc_end(gfa, from_vertex);

  for (gfa_arc_t const * arc_ptr = arc_begin(gfa, from_vertex); arc_ptr != arc_end_ptr; ++arc_ptr)
  {
    gfa_arc_t const & arc = *arc_ptr;
    bool const is_new_from_forward = (arc.w & 1u) == 0u;
    uint64_t const new_from_rid = static_cast<uint32_t>(arc.w) >> 1;
    gfa_seg_t const & new_from_segment = gfa.get_segment(new_from_rid);
    uint64_t new_from = static_cast<uint64_t>(new_from_rid) << 32; // new from RID with pos=0,strand=0

    if (!is_new_from_forward)
      new_from |= static_cast<uint64_t>(new_from_segment.len - 1) << 1 | 1; // new from position/strand

    bool const is_extending = is_exact_distance(gfa, icu, new_from, to, new_arcs, distance);

    if (is_extending)
    {
      new_arcs.push_back(arc_ptr);
      return true;
    }
  }

  return false;
}

int approximate_gfa_distance(GFA const & gfa,   // graph
                             T_icu const & icu, // The "I see you" index
                             uint64_t from,
                             uint64_t to,
                             std::vector<gfa_arc_t const *> & new_arcs,
                             int distance,
                             int const t) // threshold
{
  assert(t >= 0);
  int constexpr NO_MATCH{std::numeric_limits<int>::lowest() / 2};

  // Negative distance not allowed
  if ((distance + t) < 0)
    return NO_MATCH;

  auto dist_in_range = [=](int const d, int const distance) -> int
  {
    if (d >= (distance - t) && d <= (distance + t))
      return d - distance;
    else
      return NO_MATCH;
  };

  int const from_rid = sketch_value_rid(from);
  int const from_pos = sketch_value_pos(from);
  bool const from_strand = sketch_value_strand(from);

  int const to_rid = sketch_value_rid(to);
  int const to_pos = sketch_value_pos(to);
  bool const to_strand = sketch_value_strand(to);

  // trivial case if RIDs are the same
  if (from_rid == to_rid)
  {
    if (from_strand != to_strand)
      return NO_MATCH; // different strands

    if (from_strand == 0)
      return dist_in_range(to_pos - from_pos, distance); // forward
    else
      return dist_in_range(from_pos - to_pos, distance); // reverse
  }

  // reduce distance by the amount of the "to" segment left
  if (to_strand == 0)
  {
    // "to" is going forward, so reduce distance by its position
    to -= to_pos << 1;
    distance -= to_pos;
    /* NOTE: to_pos should be changed here as well, however it is never used again. */
    // to_pos = 0;
  }
  else
  {
    // "to" is moving in reverse, so extend it to its end
    gfa_seg_t const & to_segment = gfa.get_segment(to_rid);
    int const to_segment_remaining = to_segment.len - 1 - to_pos;
    assert(to_segment_remaining >= 0);

    to += to_segment_remaining << 1;
    distance -= to_segment_remaining;
    // NOTE: to_pos should be changed here as well, however it is never used again
    // to_pos = static_cast<int>(to_segment.len) - 1l;
  }

  // check distance again after adjusting it
  if ((distance + t) < 0)
    return NO_MATCH;

  /* NOTE: it is tempting to just check for exact distance using the icu index but
   * there are at least two problems:
   *   1) the index only holds the shortest distance.
   *   2) the index does not contain the arcs crossed/path.
   * if the index is extended to keep also this information it might grow in
   * size very fast which would slow down
   * everything.
   *
   * TLDR; Only use icu to check if the two position can (under any
   * circumstance) have this exact distance between them. */
  {
    // make sure that "from" can see "to" using the icu index.
    uint64_t const icu_key = static_cast<uint64_t>(from_rid) << 33 |    // from
                             static_cast<uint64_t>(from_strand) << 32 | //
                             static_cast<uint64_t>(to_rid) << 1 |       // to
                             static_cast<uint64_t>(to_strand);          //

    auto find_it = icu.find(icu_key);

    if (find_it == icu.end() || find_it->second > (distance + t))
      return NO_MATCH;
  }

  // go forward onto a new segment.
  gfa_seg_t const & from_segment = gfa.get_segment(from_rid);

  if (from_strand == 0)
    distance -= (from_segment.len - from_pos); // go forward (out of the back)
  else
    distance -= (from_pos + 1); // go backward

  // Negative distance indicates that the distance must be further
  if ((distance + t) < 0)
    return NO_MATCH;

  uint32_t const from_vertex{static_cast<uint32_t>(from_rid) << 1 | from_strand};
  gfa_arc_t const * arc_end_ptr = arc_end(gfa, from_vertex);

  for (gfa_arc_t const * arc_ptr = arc_begin(gfa, from_vertex); arc_ptr != arc_end_ptr; ++arc_ptr)
  {
    gfa_arc_t const & arc = *arc_ptr;
    bool const is_new_from_forward = (arc.w & 1u) == 0u;
    uint64_t const new_from_rid = static_cast<uint32_t>(arc.w) >> 1;
    gfa_seg_t const & new_from_segment = gfa.get_segment(new_from_rid);
    uint64_t new_from = static_cast<uint64_t>(new_from_rid) << 32; // new from RID with pos=0,strand=0

    if (!is_new_from_forward)
      new_from |= static_cast<uint64_t>(new_from_segment.len - 1) << 1 | 1; // new from position/strand

    // recursion
    int const is_extending = approximate_gfa_distance(gfa, icu, new_from, to, new_arcs, distance, t);

    if (is_extending != NO_MATCH)
    {
      new_arcs.push_back(arc_ptr);
      return is_extending;
    }
  }

  return NO_MATCH;
}

bool is_within_distance(GFA const & gfa, //
                        T_icu const & icu,
                        uint64_t const from,
                        uint64_t const to,
                        int const max_distance)
{
  assert(max_distance >= 0);
  int const from_rid = sketch_value_rid(from);
  int const from_pos = sketch_value_pos(from);
  bool const from_strand = sketch_value_strand(from);

  int const to_rid = sketch_value_rid(to);
  int const to_pos = sketch_value_pos(to);
  bool const to_strand = sketch_value_strand(to);

  // trivial case if RIDs are the same
  if (from_rid == to_rid)
  {
    if (from_strand != to_strand)
      return false; // different strands on the same contig don't see each other

    if (from_strand == 0)
      return to_pos >= from_pos && (to_pos - from_pos) <= max_distance; // forward
    else
      return from_pos >= to_pos && (from_pos - to_pos) <= max_distance; // reverse
  }

  // get how many bases are remaining on the vertices
  int const from_remaining = vertex_bases_remaining(gfa, from_rid, from_pos, from_strand);
  int const to_passed = vertex_bases_passed(gfa, to_rid, to_pos, to_strand);

  if ((1 + from_remaining + to_passed) > max_distance)
    return false;

  // generate the key for search the icu index
  uint64_t const icu_key = static_cast<uint64_t>(from_rid) << 33 |    // from
                           static_cast<uint64_t>(from_strand) << 32 | //
                           static_cast<uint64_t>(to_rid) << 1 |       // to
                           static_cast<uint64_t>(to_strand);          //

  // search in icu
  auto find_it = icu.find(icu_key);
  return find_it != icu.end() && (find_it->second + from_remaining + to_passed) <= max_distance;
}

int estimate_shortest_distance(GFA const & gfa,     // graph
                               T_icu const & icu,   // ICU index
                               uint64_t const from, // from this vertex position
                               uint64_t const to)   // to this vertex position
{
  int const from_rid = sketch_value_rid(from);
  int const from_pos = sketch_value_pos(from);
  int const from_strand = sketch_value_strand(from);

  int const to_rid = sketch_value_rid(to);
  int const to_pos = sketch_value_pos(to);
  int const to_strand = sketch_value_strand(to);

  // trivial case if RIDs are the same
  if (from_rid == to_rid)
  {
    if (from_strand != to_strand)
    {
      // print_debug(_HERE_, " different strand.");
      return -1;
    }

    if (from_strand == 0)
      return to_pos - from_pos; // forward
    else
      return from_pos - to_pos; // reverse
  }

  // get how many bases are remaining on the vertices
  int const from_remaining = vertex_bases_remaining(gfa, from_rid, from_pos, from_strand);
  int const to_passed = vertex_bases_passed(gfa, to_rid, to_pos, to_strand);

  // generate the key for search the icu index
  uint64_t const icu_key = static_cast<uint64_t>(from_rid) << 33 |    // from
                           static_cast<uint64_t>(from_strand) << 32 | //
                           static_cast<uint64_t>(to_rid) << 1 |       // to
                           static_cast<uint64_t>(to_strand);          //

  // search in icu
  auto find_it = icu.find(icu_key);

  print_debug(_HERE_, " from_remaining=", from_remaining, " to_passed=", to_passed);

  if (find_it == icu.end())
  {
    print_debug(_HERE_, " no icu key ", icu_key >> 32, "|", (icu_key << 32) >> 32);
    return -1;
  }

  print_debug(_HERE_, " icu value=", find_it->second);
  return find_it->second + from_remaining + to_passed;
}

T_icu make_icu_index(GFA const & gfa, int const max_distance)
{
  T_icu icu;
  int const num_segments = gfa.get_num_segments();

  for (int rid{0}; rid < num_segments; ++rid)
  {
    // forward vertex
    {
      auto const vertex_id_forward = forward_segment_id2vertex_id(rid);
      icu_recursive(icu, gfa, vertex_id_forward, vertex_id_forward, 1, max_distance);
    }

    // reverse vertex
    {
      auto const vertex_id_reverse = reverse_segment_id2vertex_id(rid);
      icu_recursive(icu, gfa, vertex_id_reverse, vertex_id_reverse, 1, max_distance);
    }
  }

  return icu;
}

} // namespace weaver
