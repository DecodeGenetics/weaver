#include <cstdio>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include <weaver/constants.hpp>
#include <weaver/io.hpp>
#include <weaver/variant.hpp>

#include "test.hpp"

namespace weaver::test
{
//! Tests for read_small_variants()
static void test_read_small_variants()
{
  /// Input streams
  weaver::hts_file_ptr in_vcf(nullptr, weaver::close_hts_file);
  weaver::tbx_t_ptr in_tbx(nullptr, weaver::close_tbx_t);
  weaver::hts_itr_t_ptr in_it(nullptr, weaver::close_hts_itr_t);

  std::string vcf_path = std::string(weaver_SOURCE_DIRECTORY) + std::string("/test/data/unit_test.vcf.gz");

  in_vcf = weaver::open_hts_file(vcf_path.c_str(), "r"); // open vcf.gz
  in_tbx = weaver::open_tbx_t(vcf_path.c_str());         // open vcf.gz.tbi
  in_it = open_hts_itr_t(in_tbx.get(), "chr20", 1, 100000);

  std::vector<Variant> small_variants = get_variants_in_a_region(in_vcf, in_tbx, in_it);

  REQUIRE(small_variants.size() == 32);

  REQUIRE(small_variants[0].pos == 499);
  REQUIRE(small_variants[1].pos == 500);
  REQUIRE(small_variants[2].pos == 4965);
  REQUIRE(small_variants[3].pos == 5165);

  REQUIRE(small_variants[0].calls.size() == 9);
  REQUIRE(small_variants[1].calls.size() == 9);

  {
    auto const & cls = small_variants[0].calls;
    REQUIRE(std::count(cls.begin(), cls.end(), 0) == 1);                             // 1 reference call
    REQUIRE(std::count(cls.begin(), cls.end(), 1) == 6);                             // 6 alt calls
    REQUIRE(std::count(cls.begin(), cls.end(), weaver::Variant::MISSING_CALL) == 2); // 2 missing calls
  }

  {
    auto const & cls = small_variants[2].calls;
    REQUIRE(cls[0] == 0);
    REQUIRE(cls[1] == 0);
    REQUIRE(cls[2] == 0);
    REQUIRE(cls[3] == 1);
    REQUIRE(cls[4] == 1);
    REQUIRE(cls[5] == 1);
    REQUIRE(cls[6] == 0);
    REQUIRE(cls[7] == 0);
  }
}

//! \cond TESTS
TEST_CASE("Tests for small_variants.cpp", "[gfa]")
{
  test_read_small_variants();
}
//! \endcond

} // namespace weaver::test
