#include <cstdio>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include <weaver/region.hpp>

#include "test.hpp"

namespace weaver::test
{
//! Tests for Region constructor
static void test_region_constructor()
{
  // Empty constructor
  {
    Region region;
    REQUIRE(region.chr.size() == 0);
    REQUIRE(region.begin == 0);
    REQUIRE(region.end == std::numeric_limits<int>::max());
    REQUIRE(!region.check()); // because missing chr
  }

  // chrN case
  {
    Region region("chr42");
    REQUIRE(region.chr == "chr42");
    REQUIRE(region.begin == 0);
    REQUIRE(region.end == std::numeric_limits<int>::max());
    REQUIRE(region.check());
  }

  // chrN:A case
  {
    Region region("chr42:1337");
    REQUIRE(region.chr == "chr42");
    REQUIRE(region.begin == 1336);
    REQUIRE(region.end == 1337);
    REQUIRE(region.size() == 1);
    REQUIRE(region.check());
  }

  // chrN:A-B case
  {
    Region region("chr42:1337-2674");
    REQUIRE(region.chr == "chr42");
    REQUIRE(region.begin == 1336);
    REQUIRE(region.end == 2674);
    REQUIRE(region.size() == 1338);
    REQUIRE(region.check());
  }

  // bad chrN:B-A case
  {
    Region region("chr42:2674-1337");
    REQUIRE(region.chr == "chr42");
    REQUIRE(region.begin == 2673);
    REQUIRE(region.end == 1337);
    REQUIRE(!region.check());
  }
}

//! \cond TESTS
TEST_CASE("Tests for region.cpp")
{
  test_region_constructor();
}
//! \endcond

} // namespace weaver::test
