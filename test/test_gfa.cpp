#include <cstdio>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include <weaver/gfa.hpp>
#include <weaver/io.hpp> // Note that it is important that this file is included, not kstring_t
#include <weaver/segment.hpp>
#include <weaver/sequence_utils.hpp> // get_reverse_complement(seq);
#include <weaver/sr_seed.hpp>

#include "test.hpp"

namespace weaver::test
{
//! Test for the GFA::GFA(std::string const &) constructor.
static void test_gfa_constructor()
{
  {
    GFA gfa = read_graph("test_small.gfa");
    REQUIRE(gfa.get_num_segments() == 6);
    REQUIRE(gfa.get_num_arcs() == 14);

    REQUIRE(gfa.get_segment(0).len == 50);
    REQUIRE(gfa.get_segment(1).len == 60);
    REQUIRE(gfa.get_segment(2).len == 40);
    REQUIRE(gfa.get_segment(3).len == 50);
    REQUIRE(gfa.get_segment(4).len == 10);
    REQUIRE(gfa.get_segment(5).len == 60);
  }

  {
    GFA gfa = read_graph("test_human_10k.gfa.gz");
    REQUIRE(gfa.get_num_segments() == 6);
    REQUIRE(gfa.get_num_arcs() == 14);

    REQUIRE(gfa.get_segment(0).len == 5180);
    REQUIRE(gfa.get_segment(1).len == 188);
    REQUIRE(gfa.get_segment(2).len == 634);
    REQUIRE(gfa.get_segment(3).len == 50);
    REQUIRE(gfa.get_segment(4).len == 3948);
    REQUIRE(gfa.get_segment(5).len == 187);
  }
}

//! Tests for GFA::get_num_arcs_from_vertex()
static void test_gfa_get_num_arcs_from_vertex()
{
  GFA gfa = read_graph("test_small.gfa");

  REQUIRE(gfa.get_num_arcs_from_vertex(0) == 2);
  REQUIRE(gfa.get_num_arcs_from_vertex(1) == 0);
  REQUIRE(gfa.get_num_arcs_from_vertex(1 << 1u | 0) == 1);
  REQUIRE(gfa.get_num_arcs_from_vertex(1 << 1u | 1) == 1);
}

//! \cond TESTS
TEST_CASE("Tests for gfa.cpp", "[gfa]")
{
  test_gfa_constructor();
  test_gfa_get_num_arcs_from_vertex();
}
//! \endcond
} // namespace weaver::test
