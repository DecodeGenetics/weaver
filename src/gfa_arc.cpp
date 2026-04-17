#include "gfa_arc.hpp"

#include <cassert>  // assert
#include <gfa.h>    // gfa_arc_a, gfa_arc_n
#include <iterator> // std::next
#include <sstream>  // std::ostringstream

#include "gfa.hpp" // GFA

namespace weaver
{
bool is_complement_arc(gfa_arc_t const * a, gfa_arc_t const * b)
{
  assert(a != nullptr);
  assert(b != nullptr);
  return (a->w == ((b->v_lv >> 32) ^ 1ull)) && (((a->v_lv >> 32) ^ 1ull) == b->w);
}

bool is_same_arc(gfa_arc_t const * a, gfa_arc_t const * b)
{
  assert(a != nullptr);
  assert(b != nullptr);
  return a->w == b->w && a->v_lv == b->v_lv;
}

gfa_arc_t const * arc_begin(GFA const & gfa, uint32_t const v)
{
  gfa_t const * g = gfa.get_graph_ptr();
  return gfa_arc_a(g, v);
}

gfa_arc_t const * arc_end(GFA const & gfa, uint32_t const v)
{
  gfa_t const * g = gfa.get_graph_ptr();
  return std::next(gfa_arc_a(g, v), gfa_arc_n(g, v));
}

gfa_arc_t const * arc_end(GFA const & gfa, uint32_t const v, gfa_arc_t const * arc_begin)
{
  gfa_t const * g = gfa.get_graph_ptr();
  return std::next(arc_begin, gfa_arc_n(g, v));
}

std::string arc_to_string(gfa_arc_t const & arc)
{
  std::ostringstream ss;
  ss << (arc.v_lv >> 33) << " " << ((arc.v_lv >> 32) & 1) << " " << (arc.w >> 1) << " " << (arc.w & 1);
  return ss.str();
}

} // namespace weaver
