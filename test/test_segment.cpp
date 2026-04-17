#include <string> // std::string

#include <weaver/constants.hpp> // weaver_SOURCE_DIRECTORY
#include <weaver/gfa.hpp>       // weaver::GFA
#include <weaver/segment.hpp>   // weaver::vertex_bases_remaining

#include <catch2/catch.hpp> // TEST_CASE

namespace weaver::test
{
//! Tests for the vertex_bases_remaining() and vertex_bases_passed() functions.
static void test_vertex_bases_remaining_or_passed()
{
  std::string graph_path = std::string(weaver_SOURCE_DIRECTORY) + "/test/data/test_small.gfa";

  GFA gfa(graph_path);
  REQUIRE(gfa.get_num_segments() == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  SECTION("segment at index 0")
  {
    REQUIRE(vertex_bases_passed(gfa, 0) == 0);
    REQUIRE(vertex_bases_remaining(gfa, 0) == 49);
    REQUIRE(vertex_bases_passed(gfa, 1ull << 1) == 1);
    REQUIRE(vertex_bases_remaining(gfa, 1ull << 1) == 48);
    REQUIRE(vertex_bases_passed(gfa, 49ull << 1) == 49);
    REQUIRE(vertex_bases_remaining(gfa, 49ull << 1) == 0);
    REQUIRE(vertex_bases_passed(gfa, 49ull << 1 | 1) == 0);
    REQUIRE(vertex_bases_remaining(gfa, 49ull << 1 | 1) == 49);
    REQUIRE(vertex_bases_passed(gfa, 1) == 49);
    REQUIRE(vertex_bases_remaining(gfa, 1) == 0);
  }

  SECTION("segment at index 1")
  {
    REQUIRE(vertex_bases_passed(gfa, 1ull << 32 | 1ull << 1) == 1);
    REQUIRE(vertex_bases_remaining(gfa, 1ull << 32 | 1ull << 1) == 58);
    REQUIRE(vertex_bases_passed(gfa, 1ull << 32 | 1ull << 1 | 1ull) == 58);
    REQUIRE(vertex_bases_remaining(gfa, 1ull << 32 | 1ull << 1 | 1ull) == 1);
    REQUIRE(vertex_bases_passed(gfa, 1ull << 32 | 59ull << 1 | 1ull) == 0);
    REQUIRE(vertex_bases_remaining(gfa, 1ull << 32 | 59ull << 1 | 1ull) == 59);
  }

  SECTION("segment at index 1, second versions of the functions")
  {
    REQUIRE(vertex_bases_passed(gfa, 1, 1, 0) == 1);
    REQUIRE(vertex_bases_remaining(gfa, 1, 1, 0) == 58);
    REQUIRE(vertex_bases_passed(gfa, 1, 1, 1) == 58);
    REQUIRE(vertex_bases_remaining(gfa, 1, 1, 1) == 1);
    REQUIRE(vertex_bases_passed(gfa, 1, 59, 1) == 0);
    REQUIRE(vertex_bases_remaining(gfa, 1, 59, 1) == 59);
  }
}

//! \cond TESTS
TEST_CASE("Tests for segment.cpp", "[gfa]")
{
  test_vertex_bases_remaining_or_passed();
}
//! \endcond
} // namespace weaver::test
