#include <sstream> // std::ostringstream

#include <parallel_hashmap/phmap.h> // T_icu

#include <weaver/constants.hpp>
#include <weaver/gfa.hpp>
#include <weaver/gfa_location.hpp>
#include <weaver/make_mmi.hpp>
#include <weaver/read_sketch.hpp>
#include <weaver/sequence_utils.hpp> // get_reverse_complement(seq);
#include <weaver/sr_seed.hpp>
#include <weaver/sr_seed_pair.hpp>

#include <catch2/catch.hpp>

namespace weaver
{
//! Tests for the function \c get_best_seed_pairs() and \c filter_sr_seeds() .
static void test_get_best_seed_pairs_and_filter()
{
  std::string graph_path = std::string(weaver_SOURCE_DIRECTORY) + "/test/data/test_small.gfa";

  GFA gfa(graph_path);
  T_icu icu = make_icu_index(gfa, 1500);

  std::vector<SRSeed> seeds1(2);
  std::vector<SRSeed> seeds2(0);

  seeds1[0].begin_read_value = 1 << 1ull;
  seeds1[0].end_read_value = 21 << 1ull;

  {
    std::vector<SRSeedPair> const & best_seed_pairs = get_best_seed_pairs(gfa, icu, seeds1, seeds2);
    REQUIRE(best_seed_pairs.size() == 2);
    REQUIRE(best_seed_pairs[0].seed_i1 == 0);
    REQUIRE(best_seed_pairs[0].seed_i2 == -1);
    REQUIRE(best_seed_pairs[1].seed_i1 == 1);
    REQUIRE(best_seed_pairs[1].seed_i2 == -1);

    std::vector<SRSeed> seeds1_cp(seeds1);
    std::vector<SRSeed> seeds2_cp(seeds2);

    filter_sr_seeds(seeds1, seeds2, best_seed_pairs);

    REQUIRE(seeds1.size() == 2);
    REQUIRE(seeds2.size() == 0);
    REQUIRE(seeds1[0].end_read_value == (21 << 1ull));
    REQUIRE(seeds1[1].end_read_value == 0);
    seeds1 = std::move(seeds1_cp);
    seeds2 = std::move(seeds2_cp);
  }

  seeds1[1].begin_read_value = 1 << 1ull;
  seeds1[1].end_read_value = 151 << 1ull;

  /*
  {
    std::vector<SRSeedPair> const & best_seed_pairs = get_best_seed_pairs(gfa, icu, seeds1, seeds2);
    REQUIRE(best_seed_pairs.size() == 2);
    REQUIRE(best_seed_pairs[0].seed_i1 == 1);
    REQUIRE(best_seed_pairs[0].seed_i2 == -1);

    std::vector<SRSeed> seeds1_cp(seeds1);
    std::vector<SRSeed> seeds2_cp(seeds2);

    filter_sr_seeds(seeds1, seeds2, best_seed_pairs);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 0);
    REQUIRE(seeds1[0].end_read_value == (151 << 1ull));
    seeds1 = std::move(seeds1_cp);
    seeds2 = std::move(seeds2_cp);
  }

  std::swap(seeds1[0], seeds1[1]);

  {
    std::vector<SRSeedPair> const & best_seed_pairs = get_best_seed_pairs(gfa, icu, seeds1, seeds2);
    REQUIRE(best_seed_pairs.size() == 2);
    REQUIRE(best_seed_pairs[0].seed_i1 == 0);
    REQUIRE(best_seed_pairs[0].seed_i2 == -1);

    std::vector<SRSeed> seeds1_cp(seeds1);
    std::vector<SRSeed> seeds2_cp(seeds2);

    filter_sr_seeds(seeds1, seeds2, best_seed_pairs);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 0);
    REQUIRE(seeds1[0].end_read_value == (151 << 1ull));
    seeds1 = std::move(seeds1_cp);
    seeds2 = std::move(seeds2_cp);
  }

  seeds2.resize(1);
  seeds2[0].begin_read_value = 1 << 1ull;
  seeds2[0].end_read_value = 81 << 1ull;

  REQUIRE(!seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
  REQUIRE(!seeds1[1].do_you_see_me(gfa, icu, seeds2[0]));

  {
    std::vector<SRSeedPair> const & best_seed_pairs = get_best_seed_pairs(gfa, icu, seeds1, seeds2);
    REQUIRE(best_seed_pairs.size() == 1);
    REQUIRE(best_seed_pairs[0].seed_i1 == 0);
    REQUIRE(best_seed_pairs[0].seed_i2 == 0);

    std::vector<SRSeed> seeds1_cp(seeds1);
    std::vector<SRSeed> seeds2_cp(seeds2);

    filter_sr_seeds(seeds1, seeds2, best_seed_pairs);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);
    REQUIRE(seeds1[0].end_read_value == (151 << 1ull));
    seeds1 = std::move(seeds1_cp);
    seeds2 = std::move(seeds2_cp);
  }

  seeds2[0].begin_ref_value = 41 << 1ull;
  seeds2[0].end_ref_value = 11 << 1ull | 1ull;
  seeds1[1].begin_ref_value = 1 << 1ull | 1ull;
  seeds1[1].end_ref_value = 5 << 1ull;

  REQUIRE(!seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
  REQUIRE(seeds1[1].do_you_see_me(gfa, icu, seeds2[0]));

  {
    std::vector<SRSeedPair> const & best_seed_pairs = get_best_seed_pairs(gfa, icu, seeds1, seeds2);
    REQUIRE(best_seed_pairs.size() == 2);
    REQUIRE(best_seed_pairs[0].seed_i1 == 0);
    REQUIRE(best_seed_pairs[0].seed_i2 == 0);
    REQUIRE(best_seed_pairs[1].seed_i1 == 1);
    REQUIRE(best_seed_pairs[1].seed_i2 == 0);
    REQUIRE(best_seed_pairs[0].score > best_seed_pairs[1].score);

    std::vector<SRSeed> seeds1_cp(seeds1);
    std::vector<SRSeed> seeds2_cp(seeds2);

    filter_sr_seeds(seeds1, seeds2, best_seed_pairs);

    REQUIRE(seeds1.size() == 2);
    REQUIRE(seeds2.size() == 1);
    seeds1 = std::move(seeds1_cp);
    seeds2 = std::move(seeds2_cp);
  }
  */
}

//! \cond TESTS
TEST_CASE("Tests for sr_seed_pair.cpp.", "[gfa]")
{
  test_get_best_seed_pairs_and_filter();
}
//! \endcond

} // namespace weaver
