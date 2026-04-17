#include <sstream> // std::ostringstream

#include <parallel_hashmap/phmap.h> // T_icu

#include <weaver/constants.hpp>
#include <weaver/gfa.hpp>
#include <weaver/gfa_location.hpp>
#include <weaver/make_mmi.hpp>
#include <weaver/read_sketch.hpp>
#include <weaver/sequence_utils.hpp> // get_reverse_complement(seq);
#include <weaver/sr_seed.hpp>

#include <catch2/catch.hpp>

namespace weaver
{
//! Tests for the function \c get_seeds_with_pair()
static void test_get_seeds_with_pair()
{
  std::string graph_path = std::string(weaver_SOURCE_DIRECTORY) + "/test/data/test_small.gfa";

  GFA gfa(graph_path);
  auto const num_segments = gfa.get_num_segments();
  REQUIRE(num_segments == 6);
  REQUIRE(gfa.get_num_arcs() == 14);
  int const k{23};
  int const w{5};

  MMI mmi = make_mmi_index(gfa, "", k, w);
  REQUIRE(mmi.num_keys() == 86);
  T_icu icu = make_icu_index(gfa, 1500);

  std::string s1 = "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC";
  std::string s2 = "ATCCTTGAAAATAAACACTAAAAATACATCCAAATGTTTAACCCAGTTTGGTCATTTTTT";
  std::string s3 = "TAAAAGTCAAAGATCCCTGGAGGGACAGGTGGGGGTGAGG";
  std::string s4 = "CTTAGACTTATTGTGAAATATGTAAGTGCTTATCTTGTAAAAATAGTATT";
  std::string s5 = "GAGCACAGGG";
  std::string s6 = "CCCCCCCAGGTATTAGGAAATAAAGCACAGAGAAAAATAAGAATGCTGAGCCGGGATCAT"; // length 60

  std::string s5_rev = get_reverse_complement(s5);
  std::string s6_rev = get_reverse_complement(s6);

  // Sequence and reverse sequence on same rid
  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s6, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s6_rev, num_segments, mmi, k, w);

    // s6 forward
    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds1[0].get_read_length() == 36);
    REQUIRE(seeds1[0].begin_ref_value == (5ull << 32 | 11ull << 1 | 1ull));
    REQUIRE(seeds1[0].end_ref_value == (5ull << 32 | 46ull << 1));
    REQUIRE(seeds1[0].begin_read_value == (6ull << 32 | 11ull << 1 | 1ull));
    REQUIRE(seeds1[0].end_read_value == (6ull << 32 | 46ull << 1));

    std::string ref_seq = seeds1[0].get_ref_sequence();
    std::string read_seq = seeds1[0].get_read_sequence(s6);
    REQUIRE(ref_seq == read_seq);
    REQUIRE(ref_seq == "ATTAGGAAATAAAGCACAGAGAAAAATAAGAATGCT");

    // reverse
    REQUIRE(seeds2.size() == 1);
    REQUIRE(seeds2[0].get_read_length() == 36);
    REQUIRE(seeds2[0].begin_ref_value == (5ull << 32 | 46ull << 1));
    REQUIRE(seeds2[0].end_ref_value == (5ull << 32 | 11ull << 1 | 1ull));
    REQUIRE(seeds2[0].begin_read_value == (6ull << 32 | 13ull << 1 | 1ull));
    REQUIRE(seeds2[0].end_read_value == (6ull << 32 | 48ull << 1));

    std::string ref_seq_rev = seeds2[0].get_ref_sequence();
    std::string read_seq_rev = seeds2[0].get_read_sequence(s6_rev);

    REQUIRE(ref_seq_rev == read_seq_rev);
    REQUIRE(ref_seq_rev.size() == 36);
    REQUIRE(ref_seq == get_reverse_complement(ref_seq_rev));
  }

  {
    // s5 is too small
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s5, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s5_rev, num_segments, mmi, k, w);

    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 0);
    REQUIRE(seeds2.size() == 0);
  }

  {
    // s4>s5
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s4 + s5, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 =
      get_read_sketches(get_reverse_complement(s4 + s5), num_segments, mmi, k, w);

    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds1[0].get_read_length() == 38);
    REQUIRE(seeds2.size() == 1);
    REQUIRE(seeds2[0].get_read_length() == 38);
  }
}

static void test_do_you_see_me()
{
  std::string graph_path = std::string(weaver_SOURCE_DIRECTORY) + "/test/data/test_small.gfa";

  GFA gfa(graph_path);
  long const num_segments = gfa.get_num_segments();
  REQUIRE(num_segments == 6);
  REQUIRE(gfa.get_num_arcs() == 14);
  int const k{23};
  int const w{5};

  MMI mmi = make_mmi_index(gfa, "", k, w);
  REQUIRE(mmi.num_keys() == 86);
  T_icu icu = make_icu_index(gfa, 1500);

  std::string s1 = "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC";
  std::string s2 = "ATCCTTGAAAATAAACACTAAAAATACATCCAAATGTTTAACCCAGTTTGGTCATTTTTT";
  std::string s3 = "TAAAAGTCAAAGATCCCTGGAGGGACAGGTGGGGGTGAGG";
  std::string s4 = "CTTAGACTTATTGTGAAATATGTAAGTGCTTATCTTGTAAAAATAGTATT";
  std::string s5 = "GAGCACAGGG";
  std::string s6 = "CCCCCCCAGGTATTAGGAAATAAAGCACAGAGAAAAATAAGAATGCTGAGCCGGGATCAT"; // length 60

  std::string s1_rev = get_reverse_complement(s1);
  std::string s2_rev = get_reverse_complement(s2);

  std::string s1_sub = s1.substr(0, 40);
  std::string s1_rev_sub = s1_rev.substr(0, 40);

  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s1_sub, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s1_rev_sub, num_segments, mmi, k, w);

    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);

    REQUIRE(seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
    REQUIRE(seeds2[0].do_you_see_me(gfa, icu, seeds1[0]));
  }

  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s1_rev_sub, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s1_sub, num_segments, mmi, k, w);

    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);

    REQUIRE(seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
    REQUIRE(seeds2[0].do_you_see_me(gfa, icu, seeds1[0]));
  }

  // Sequence and reverse sequence on the next rid, they see each other
  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s1, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s2_rev, num_segments, mmi, k, w);

    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);

    REQUIRE(seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
    REQUIRE(seeds2[0].do_you_see_me(gfa, icu, seeds1[0]));
  }

  // Reversed sequence and forward sequence on the next rid, they don't see each other
  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s1_rev, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s2, num_segments, mmi, k, w);

    // s6 forward
    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);

    REQUIRE(!seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
    REQUIRE(!seeds2[0].do_you_see_me(gfa, icu, seeds1[0]));
  }

  // Both forward, they don't see each other
  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s1, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s2, num_segments, mmi, k, w);

    // s6 forward
    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);

    REQUIRE(!seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
    REQUIRE(!seeds2[0].do_you_see_me(gfa, icu, seeds1[0]));
  }

  // Both reverse, they don't see each other
  {
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(s1_rev, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(s2_rev, num_segments, mmi, k, w);

    // s6 forward
    std::vector<SRSeed> seeds1;
    std::vector<SRSeed> seeds2;

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu);

    REQUIRE(seeds1.size() == 1);
    REQUIRE(seeds2.size() == 1);

    REQUIRE(!seeds1[0].do_you_see_me(gfa, icu, seeds2[0]));
    REQUIRE(!seeds2[0].do_you_see_me(gfa, icu, seeds1[0]));
  }
}

static void test_get_seed_estimated_score()
{
  std::vector<SRSeed> seeds(3);
  seeds[0].begin_read_value = 3 << 1ull;
  seeds[0].end_read_value = 11 << 1ull; // length 11-3+1=9
  seeds[0].num_cuts = 0;
  REQUIRE(seeds[0].get_est_score() == 9); // score = 9-0=9

  seeds[1].begin_read_value = 3 << 1ull;
  seeds[1].end_read_value = 12 << 1ull; // length 12-3+1=10
  seeds[1].num_cuts = 1;                // score = 10-4=6
  REQUIRE(seeds[1].get_est_score() == 6);

  seeds[2].begin_read_value = 1 << 1ull;
  seeds[2].end_read_value = 20 << 1ull; // length 20-1+1=20
  seeds[2].num_cuts = 2;                // score = 20-8=12
  REQUIRE(seeds[2].get_est_score() == 12);

  std::vector<int> orders = get_seed_est_score_sorted_order_indices(seeds);
  REQUIRE(orders.size() == seeds.size());
  REQUIRE(orders[0] == 2); // best score is at index 2
  REQUIRE(orders[1] == 0);
  REQUIRE(orders[2] == 1); // worst score is at index 1
}

//! \cond TESTS
TEST_CASE("Tests for sr_seed.cpp.", "[gfa]")
{
  test_get_seeds_with_pair();
  test_do_you_see_me();
  test_get_seed_estimated_score();
}
//! \endcond

} // namespace weaver
