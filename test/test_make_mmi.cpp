#include <cstdio>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include <weaver/constants.hpp>
#include <weaver/make_mmi.hpp>

#include "test.hpp"

namespace weaver::test
{
/*!
 * @brief Tests for make_mmi_index()
 *
 * @details
 * The functionality is already tested for the common cases with the black box test but here a VCF file that covers more
 * edge cases is tested.
 */
static void test_make_mmi_index()
{
  std::string gfa_path = std::string(weaver_SOURCE_DIRECTORY) + std::string("/test/data/test_human_10k.gfa.gz");
  std::string vcf_path = std::string(weaver_SOURCE_DIRECTORY) + std::string("/test/data/unit_test.vcf.gz");
  GFA gfa(gfa_path);

  weaver::MMI mmi = make_mmi_index(gfa, vcf_path, /*k=*/-1, /*w=*/-1);

  REQUIRE(mmi.haplotypes.size() == 1);

  auto const & vars = mmi.haplotypes[0];

  REQUIRE(vars.size() == 41);
  REQUIRE(vars[0].calls.size() == 9);
}

//! \cond TESTS
TEST_CASE("Tests for make_mmi.cpp", "[mmi]")
{
  test_make_mmi_index();
}
//! \endcond
} // namespace weaver::test
