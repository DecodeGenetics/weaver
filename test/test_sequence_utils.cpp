#include <string> // std::string

#include <weaver/sequence_utils.hpp>

#include <catch2/catch.hpp>

namespace weaver::test
{
//! Tests for complement()
static void test_complement()
{
  REQUIRE(complement('A') == 'T');
  REQUIRE(complement('T') == 'A');
  REQUIRE(complement('C') == 'G');
  REQUIRE(complement('A') == 'T');
  REQUIRE(complement('N') == 'N');
  REQUIRE(complement('a') == 'T');
}

//! Tests for nt5_char_to_ull() and complement_nt5_char_to_ull
static void test_nt5_char_to_ull_and_complement_nt5_char_to_ull()
{
  // acgtn to their index
  REQUIRE(nt5_char_to_ull('A') == 0);
  REQUIRE(nt5_char_to_ull('C') == 1);
  REQUIRE(nt5_char_to_ull('G') == 2);
  REQUIRE(nt5_char_to_ull('T') == 3);
  REQUIRE(nt5_char_to_ull('U') == 3);
  REQUIRE(nt5_char_to_ull('N') == 4);

  // complement version
  REQUIRE(complement_nt5_char_to_ull('A') == 3);
  REQUIRE(complement_nt5_char_to_ull('C') == 2);
  REQUIRE(complement_nt5_char_to_ull('G') == 1);
  REQUIRE(complement_nt5_char_to_ull('T') == 0);
  REQUIRE(complement_nt5_char_to_ull('U') == 0);
  REQUIRE(complement_nt5_char_to_ull('N') == 4);

  for (char c : {'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'})
  {
    REQUIRE(complement_nt5_char_to_ull(c) == nt5_char_to_ull(complement(c)));
  }
}

//! Tests for isACGT()
static void test_isACGT()
{
  for (char c : {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'})
  {
    REQUIRE(isACGT(c));
  }

  REQUIRE(!isACGT('N'));
  REQUIRE(!isACGT('n'));
  REQUIRE(!isACGT('X'));
  REQUIRE(!isACGT('B'));
}

//! Tests for get_reverse_complement()
static void test_get_reverse_complement()
{
  REQUIRE(get_reverse_complement(std::string("A")) == "T");
  REQUIRE(get_reverse_complement(std::string("AG")) == "CT");
  REQUIRE(get_reverse_complement(std::string("AGN")) == "NCT");
  REQUIRE(get_reverse_complement(std::string("AGNGTC")) == "GACNCT");
}

//! Tests for count_prefix_mismatches()
static void test_count_prefix_mismatches()
{
  REQUIRE(count_prefix_mismatches(std::string("A"), std::string("A")) == 0);
  REQUIRE(count_prefix_mismatches(std::string("A"), std::string("T")) == 1);
  REQUIRE(count_prefix_mismatches(std::string("A"), std::string("")) == 0);
  REQUIRE(count_prefix_mismatches(std::string(""), std::string("A")) == 0);
  REQUIRE(count_prefix_mismatches(std::string(""), std::string("ACGTACGT")) == 0);
  REQUIRE(count_prefix_mismatches(std::string("ACGAACGA"), std::string("ACGTACGT")) == 2);
  REQUIRE(count_prefix_mismatches(std::string("ACGAACGAAAA"), std::string("ACGTACGT")) == 2);
  REQUIRE(count_prefix_mismatches(std::string("ACGTACGTAACGAAAA"), std::string("ACGTACGT")) == 0);
  REQUIRE(count_prefix_mismatches(std::string("ACGTACGT"), std::string("ACGTACGTAACGAAAA")) == 0);
}

static void test_count_prefix_mismatches_s2rev()
{
  REQUIRE(count_prefix_mismatches_s2rev(std::string("A"), std::string("A")) == 1);
  REQUIRE(count_prefix_mismatches_s2rev(std::string("A"), std::string("T")) == 0);
  REQUIRE(count_prefix_mismatches_s2rev(std::string("A"), std::string("")) == 0);
  REQUIRE(count_prefix_mismatches_s2rev(std::string(""), std::string("A")) == 0);
  REQUIRE(count_prefix_mismatches_s2rev(std::string(""), std::string("ACGTACGT")) == 0);
  REQUIRE(count_prefix_mismatches_s2rev(std::string("ACGAACGA"), get_reverse_complement(std::string("ACGTACGT"))) == 2);

  REQUIRE(count_prefix_mismatches_s2rev(std::string("ACGTACGTAACGAAAA"),                // sequence 1
                                        get_reverse_complement(std::string("ACGTACGT")) // sequence 2
                                        ) == 0);

  REQUIRE(count_prefix_mismatches_s2rev(std::string("ACGTACGT"),                                // sequence 1
                                        get_reverse_complement(std::string("ACGTACGTAACGAAAA")) // sequence 2
                                        ) == 0);
}

//! \cond TESTS
TEST_CASE("Tests for sequence_utils", "")
{
  test_complement();
  test_count_prefix_mismatches();
  test_count_prefix_mismatches_s2rev();
  test_isACGT();
  test_get_reverse_complement();
  test_nt5_char_to_ull_and_complement_nt5_char_to_ull();
}
//! \endcond
} // namespace weaver::test
