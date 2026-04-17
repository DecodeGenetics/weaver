#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>

#include "gfa.hpp"
#include "icu.hpp"
#include "logging.hpp"
#include "sketch_value.hpp"
#include "sr_seed.hpp"

namespace weaver
{
bool seed_extend_if_exact_distance(GFA const & gfa,     // graph
                                   SRSeed & seed,       // from this seed, will be modified iff match is found
                                   T_icu const & icu,   // the "I see you" index
                                   uint64_t to,         // to this point
                                   uint64_t read_value, //
                                   int distance)        // distance to check
{
  assert(distance > 0);
  assert(sketch_value_rid(seed.end_ref_value) < gfa.get_num_segments());

  uint64_t const from = seed.end_ref_value;
  std::vector<gfa_arc_t const *> new_arcs;
  bool const found_exact_distance = is_exact_distance(gfa, icu, from, to, new_arcs, distance);

  if (found_exact_distance && gfa.are_arcs_with_same_rank_as_segments(new_arcs))
  {
    seed.arcs.insert(seed.arcs.end(), new_arcs.rbegin(), new_arcs.rend());
    seed.extend_end(distance, to, read_value & ~1ull); // Unset strand of read_value

    if (seed.is_extending_cut)
    {
      ++seed.num_cuts;
      seed.is_extending_cut = false; // continue from here
    }

    return true;
  }

  return false;
}

bool seed_extend_if_approximate_distance(GFA const & gfa,     // graph
                                         SRSeed & seed,       // from this seed, will be modified iff match is found
                                         T_icu const & icu,   // the "I see you" index
                                         uint64_t to,         // to this point
                                         uint64_t read_value, //
                                         int distance)        // distance to check
{
  assert(distance > 0);
  assert(sketch_value_rid(seed.end_ref_value) < gfa.get_num_segments());

  uint64_t const from = seed.end_ref_value;
  std::vector<gfa_arc_t const *> new_arcs;
  int const t = 1 + distance / 6; // threshold
  assert(t <= distance);
  int const off = approximate_gfa_distance(gfa, icu, from, to, new_arcs, distance, t);

  if (off > (std::numeric_limits<int>::lowest() / 2) && gfa.are_arcs_with_same_rank_as_segments(new_arcs))
  {
    assert(off != 0);
    seed.arcs.insert(seed.arcs.end(), new_arcs.rbegin(), new_arcs.rend());
    seed.extend_end(distance + off, to, read_value & ~1ull); // Unset strand of read_value
    assert(seed.is_extending_cut);
    seed.num_cuts += 2;            // approximate match counts as one extra cut
    seed.is_extending_cut = false; // continue from here
    return true;
  }

  return false;
}
} // namespace weaver
