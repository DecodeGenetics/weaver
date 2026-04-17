/*!
 * @file gfa.cpp
 * @brief Implements the GFA class.
 */

#include "gfa.hpp"

#include <cassert>
#include <gfa.h>

#include <weaver/constants.hpp>

#include "gfa_arc.hpp"      // arc_begin, arc_end
#include "gfa_location.hpp" // GFALocation
#include "hashmap.hpp"
#include "logging.hpp"        // print_error
#include "sketch_value.hpp"   // sketch_value_rid
#include "stable_contigs.hpp" // add_all_stable_contigs_from_graph

#include <gfa-priv.h>

namespace weaver
{
GFA::GFA(std::string const & _fn)
{
  g = gfa_read(_fn.c_str());

  if (!g)
  {
    print_error(_HERE_, " Could not read GFA:", _fn);
    std::exit(1);
  }

  fn = _fn;
  GFALocation::gfa = this; // Set the GFALocation pointer to this newly opened graph
  set_stable_contigs_from_graph(*this);
}

GFA::GFA(GFA && other) noexcept
{
  g = other.g;
  other.g = nullptr;
  GFALocation::gfa = this;
}

GFA::~GFA()
{
  gfa_destroy(g);
  GFALocation::gfa = nullptr;
}

bool GFA::is_arc_with_same_rank_as_segments(gfa_arc_t const * arc)
{
  return (arc->rank == get_segment(arc->w >> 1ull).rank) && (arc->rank == get_segment(arc->v_lv >> 33ull).rank);
}

bool GFA::are_arcs_with_same_rank_as_segments(std::vector<gfa_arc_t const *> const & arcs) const
{
  for (gfa_arc_t const * arc : arcs)
  {
    assert(arc != nullptr);

    if ((arc->rank != get_segment(arc->w >> 1ull).rank) || (arc->rank != get_segment(arc->v_lv >> 33ull).rank))
      return false;
  }

  return true;
}

gfa_t const & GFA::get_graph() const
{
  assert(g != nullptr);
  return *g;
}

gfa_t const * GFA::get_graph_ptr() const
{
  assert(g != nullptr);
  return g;
}

gfa_arc_t const * GFA::get_arc_to_next_stable_sequence(uint32_t vertex_id) const
{
  gfa_arc_t const * arc_end_ptr = arc_end(*this, vertex_id);

  for (gfa_arc_t const * arc_ptr = arc_begin(*this, vertex_id); arc_ptr != arc_end_ptr; ++arc_ptr)
  {
    if (arc_ptr->rank == 0)
      return arc_ptr;
  }

  return nullptr;
}

int64_t GFA::get_approximate_stable_position(uint64_t const value) const
{
  assert(g != nullptr);
  assert(value != std::numeric_limits<uint64_t>::max());

  gfa_seg_t const & segment = get_segment(sketch_value_rid(value));

  return (static_cast<uint64_t>(segment.snid) << 32ull) | //
         (static_cast<uint64_t>(sketch_value_pos(value) + segment.soff) >> 10ull);
}

uint64_t GFA::get_lowest_possible_stable_position(uint64_t const value) const
{
  assert(g != nullptr);
  assert(value != std::numeric_limits<uint64_t>::max());

  gfa_seg_t const & segment = get_segment(sketch_value_rid(value));
  uint64_t ret = static_cast<uint64_t>(segment.snid) << 32ull;

  // if (segment.rank == 0)
  return ret | static_cast<uint64_t>(sketch_value_pos(value) + segment.soff);
  // else
  //  return ret | static_cast<uint64_t>(segment.soff);
}

int GFA::get_num_segments() const
{
  assert(g != nullptr);
  return g->n_seg;
}

int GFA::get_num_stable_segments() const
{
  assert(g != nullptr);
  return g->n_sseq;
}

int GFA::get_num_arcs() const
{
  assert(g != nullptr);
  return g->n_arc;
}

gfa_seg_t const & GFA::get_segment(uint32_t const rid) const
{
#ifndef NDEBUG
  assert(g != nullptr);

  if (static_cast<int>(rid) >= get_num_segments())
  {
    print_warning(_HERE_, " No RID=", rid);
    log_singleton->flush_stream();
    assert(static_cast<int>(rid) < get_num_segments()); // will trigger
  }

#endif // NDEBUG

  return g->seg[rid];
}

std::string_view GFA::view_segment(uint32_t const rid) const
{
#ifndef NDEBUG
  assert(g);

  if (static_cast<int64_t>(rid) >= get_num_segments())
  {
    print_warning(_HERE_, " No RID=", rid);
    log_singleton->flush_stream();
  }

  assert(static_cast<int64_t>(rid) < get_num_segments());
#endif // NDEBUG

  gfa_seg_t const & segment = g->seg[rid];

#ifndef NDEBUG
  assert(static_cast<int>(strlen(segment.seq)) == segment.len);
#endif // NDEBUG

  return std::string_view(segment.seq, segment.len);
}

gfa_sseq_t const & GFA::get_stable_sequence(gfa_seg_t const & segment) const
{
  assert(segment.snid < static_cast<int>(g->n_sseq));
  return g->sseq[segment.snid];
}

int GFA::get_num_arcs_from_vertex(uint32_t const v) const
{
  return gfa_arc_n(g, v);
}

bool GFA::mmi_value_order(uint64_t const a, uint64_t const b) const
{
  return get_lowest_possible_stable_position(a) < get_lowest_possible_stable_position(b);
}

} // namespace weaver
