#include <string>

#include <parallel_hashmap/phmap.h> // T_icu

#include <weaver/constants.hpp>
#include <weaver/gfa.hpp>
#include <weaver/icu.hpp>
#include <weaver/logging.hpp>

#include <catch2/catch.hpp>

namespace weaver::test
{
//! Tests the \c is_within_distance() function.
static void test_is_within_distance()
{
  std::ostringstream graph_path;
  graph_path << weaver_SOURCE_DIRECTORY << "/test/data/test_small.gfa";

  GFA gfa(graph_path.str());
  REQUIRE(gfa.get_num_segments() == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  T_icu icu = make_icu_index(gfa, 1500);

  {
    REQUIRE(is_within_distance(gfa, icu, 0, 0, 1500));
    REQUIRE(is_within_distance(gfa, icu, 0, 0, 0));
    REQUIRE(is_within_distance(gfa, icu, 1ull << 32 | 2ull << 1, 1ull << 32 | 2 << 1, 0));
    REQUIRE(is_within_distance(gfa, icu, 1ull << 32 | 2ull << 1 | 1, 1ull << 32 | 2 << 1 | 1, 0));
  }

  {
    REQUIRE(is_within_distance(gfa, icu, 0, 2 << 1, 1500));
    REQUIRE(is_within_distance(gfa, icu, 0, 2 << 1, 2));
    REQUIRE(!is_within_distance(gfa, icu, 0, 2 << 1, 1));
    REQUIRE(!is_within_distance(gfa, icu, 1, 2 << 1 | 1, 1));
    REQUIRE(!is_within_distance(gfa, icu, 1, 2 << 1 | 1, 1500));
  }

  {
    // shortest path s1++s2
    {
      // to pos s2:0
      REQUIRE(is_within_distance(gfa, icu, 0, 1ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 0, 1ull << 32, 50));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32, 49));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32, 0));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 1, 1500));

      // to pos s2:10
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 10ull << 1, 50));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 10ull << 1, 59));
      REQUIRE(is_within_distance(gfa, icu, 0, 1ull << 32 | 10ull << 1, 60));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 10ull << 1 | 1, 60));

      // to pos s2:59
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 59ull << 1, 50));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 59ull << 1, 108));
      REQUIRE(is_within_distance(gfa, icu, 0, 1ull << 32 | 59ull << 1, 109));
      REQUIRE(is_within_distance(gfa, icu, 0, 1ull << 32 | 59ull << 1, 110));
      REQUIRE(!is_within_distance(gfa, icu, 0, 1ull << 32 | 59ull << 1 | 1, 110));

      // from pos s1:40
      REQUIRE(!is_within_distance(gfa, icu, 40ull << 1, 1ull << 32 | 59ull << 1, 10));
      REQUIRE(!is_within_distance(gfa, icu, 40ull << 1, 1ull << 32 | 59ull << 1, 68));
      REQUIRE(is_within_distance(gfa, icu, 40ull << 1, 1ull << 32 | 59ull << 1, 69));
      REQUIRE(is_within_distance(gfa, icu, 40ull << 1, 1ull << 32 | 59ull << 1, 70));
      REQUIRE(!is_within_distance(gfa, icu, 40ull << 1, 1ull << 32 | 59ull << 1 | 1, 70));

      // switch from and to
      REQUIRE(!is_within_distance(gfa, icu, 1ull << 32, 0, 1500));
      REQUIRE(is_within_distance(gfa, icu, 1ull << 32 | 1, 1, 1500));
      REQUIRE(is_within_distance(gfa, icu, 1ull << 32 | 1, 1, 50));
      REQUIRE(!is_within_distance(gfa, icu, 1ull << 32 | 1, 1, 49));
    }

    // shortest path s1++s2++s3
    {
      // to pos s3:0
      REQUIRE(is_within_distance(gfa, icu, 0, 2ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 0, 2ull << 32, 110));
      REQUIRE(!is_within_distance(gfa, icu, 0, 2ull << 32, 109));
      REQUIRE(!is_within_distance(gfa, icu, 0, 2ull << 32, 0));
      REQUIRE(!is_within_distance(gfa, icu, 0, 2ull << 32 | 1, 1500));

      // from pos s1:49
      REQUIRE(is_within_distance(gfa, icu, 49ull << 1, 2ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 49ull << 1, 2ull << 32, 61));
      REQUIRE(!is_within_distance(gfa, icu, 49ull << 1, 2ull << 32, 60));
      REQUIRE(!is_within_distance(gfa, icu, 49ull << 1, 2ull << 32, 0));
      REQUIRE(!is_within_distance(gfa, icu, 49ull << 1, 2ull << 32 | 1, 1500));
    }

    // shortest path s1++s2++s3++s4
    {
      REQUIRE(is_within_distance(gfa, icu, 0, 3ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 0, 3ull << 32, 150));
      REQUIRE(!is_within_distance(gfa, icu, 0, 3ull << 32, 149));
    }

    // shorest path s1++s2++s3++s5
    {
      REQUIRE(is_within_distance(gfa, icu, 0, 4ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 0, 4ull << 32, 150));
      REQUIRE(!is_within_distance(gfa, icu, 0, 4ull << 32, 149));
    }

    // shortest path s1+-s6
    {
      REQUIRE(!is_within_distance(gfa, icu, 0, 5ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 0, 5ull << 32 | 1, 1500));
      REQUIRE(is_within_distance(gfa, icu, 0, 5ull << 32 | 1, 110));
      REQUIRE(is_within_distance(gfa, icu, 0, 5ull << 32 | 1, 109));
      REQUIRE(!is_within_distance(gfa, icu, 0, 5ull << 32 | 1, 108));
      REQUIRE(is_within_distance(gfa, icu, 0, 5ull << 32 | 59ull << 1 | 1, 50));
      REQUIRE(!is_within_distance(gfa, icu, 0, 5ull << 32 | 59ull << 1 | 1, 49));
    }

    // shortest path s3++s5
    {
      REQUIRE(is_within_distance(gfa, icu, 2ull << 32, 4ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 2ull << 32, 4ull << 32, 40));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32, 4ull << 32, 39));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32, 4ull << 32 | 1, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32 | 1, 4ull << 32 | 1, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32 | 1, 4ull << 32, 1500));
    }

    // shortest path s4++s5
    {
      REQUIRE(is_within_distance(gfa, icu, 0, 4ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 3ull << 32, 4ull << 32, 50));
      REQUIRE(!is_within_distance(gfa, icu, 3ull << 32, 4ull << 32, 49));
      REQUIRE(!is_within_distance(gfa, icu, 1, 4ull << 32, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 1, 4ull << 32 | 1, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 0, 4ull << 32 | 1, 1500));
    }

    // shortest path s6-+s3
    {
      REQUIRE(is_within_distance(gfa, icu, 5ull << 32 | 1, 2ull << 32, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 5ull << 32 | 1, 2ull << 32 | 1, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 5ull << 32, 2ull << 32 | 1, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 5ull << 32, 2ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 5ull << 32 | 1, 2ull << 32, 1));
      REQUIRE(!is_within_distance(gfa, icu, 5ull << 32 | 1, 2ull << 32, 0));

      // swap from and to
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32, 5ull << 32 | 1, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32, 5ull << 32, 1500));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32 | 1, 5ull << 32 | 1, 1500));
      REQUIRE(is_within_distance(gfa, icu, 2ull << 32 | 1, 5ull << 32, 1500));
      REQUIRE(is_within_distance(gfa, icu, 2ull << 32 | 1, 5ull << 32, 1));
      REQUIRE(!is_within_distance(gfa, icu, 2ull << 32 | 1, 5ull << 32, 0));
    }
  }
}

//! Tests for the function \c make_icu_index()
static void test_make_icu_index()
{
  std::string graph_path = static_cast<std::string>(weaver_source_dir) + std::string("/test/data/test_small.gfa");
  GFA gfa(graph_path);

  // Maximum distance 1500
  {
    T_icu icu = make_icu_index(gfa, 1500);
    REQUIRE(icu.size() == 28);

    for (auto it = icu.begin(); it != icu.end(); ++it)
    {
      REQUIRE(it->second <= 1500);
    }
  }

  // Maximum distance 50
  {
    T_icu icu = make_icu_index(gfa, 50);
    REQUIRE(icu.size() == 22);

    for (auto it = icu.begin(); it != icu.end(); ++it)
    {
      REQUIRE(it->second <= 50);
    }
  }
}

//! \cond TESTS
TEST_CASE("Tests for icu.cpp", "[gfa][icu]")
{
  test_is_within_distance();
  test_make_icu_index();
}
//! \endcond
} // namespace weaver::test
