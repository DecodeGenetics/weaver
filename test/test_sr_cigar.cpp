#include <parallel_hashmap/phmap.h>

#include <weaver/constants.hpp>
#include <weaver/sam_record.hpp>
#include <weaver/sr_cigar.hpp>

#include <catch2/catch.hpp>

namespace weaver
{
//! Tests for the function \c left_align_deletion().
static void test_left_align_deletion()
{
  {
    // test case with 1bp left aligning
    SAMRecord record;
    // GCAAAATT
    // GCA-AATT
    std::string ref{"GCAAAATT"};
    std::string read{"GCAAATT"};

    record.seq = read;

    record.append_cigar(3, paw::CigarOperation::MATCH);
    record.append_cigar(1, paw::CigarOperation::DELETION);
    record.append_cigar(4, paw::CigarOperation::MATCH);

    REQUIRE(record.is_cigar_valid());
    REQUIRE(record.cig.size() == 3);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::DELETION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 3);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 4);

    auto ref_it = ref.cbegin() + record.cig[0].count;
    auto read_it = read.cbegin() + record.cig[0].count;
    left_align_deletion(record, 1, ref_it, read_it);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::DELETION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 2);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 5);
  }

  {
    // test case with 3bp left aligning
    SAMRecord record;
    // GCAAAATT
    // GCAAA-TT
    std::string ref{"GCAAAATT"};
    std::string read{"GCAAATT"};

    record.seq = read;

    record.append_cigar(5, paw::CigarOperation::MATCH);
    record.append_cigar(1, paw::CigarOperation::DELETION);
    record.append_cigar(2, paw::CigarOperation::MATCH);

    REQUIRE(record.is_cigar_valid());
    REQUIRE(record.cig.size() == 3);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::DELETION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 5);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 2);

    auto ref_it = ref.cbegin() + record.cig[0].count;
    auto read_it = read.cbegin() + record.cig[0].count;
    left_align_deletion(record, 1, ref_it, read_it);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::DELETION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 2);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 5);
  }

  {
    // test case with no left aligning
    SAMRecord record;
    // GCAAAATT
    // GCAAAA-T
    std::string ref{"GCAAAATT"};
    std::string read{"GCAAAAT"};

    record.seq = read;

    record.append_cigar(6, paw::CigarOperation::MATCH);
    record.append_cigar(1, paw::CigarOperation::DELETION);
    record.append_cigar(1, paw::CigarOperation::MATCH);

    REQUIRE(record.is_cigar_valid());
    REQUIRE(record.cig.size() == 3);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::DELETION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 6);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 1);

    auto ref_it = ref.cbegin() + record.cig[0].count;
    auto read_it = read.cbegin() + record.cig[0].count;
    left_align_deletion(record, 1, ref_it, read_it);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::DELETION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 6);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 1);
  }
}

//! Tests for the function \c left_align_insertion().
static void test_left_align_insertion()
{
  {
    // test case with 1bp left aligning
    SAMRecord record;
    std::string ref{"GCAAATT"};
    std::string read{"GCAAAATT"};

    record.seq = read;

    record.append_cigar(3, paw::CigarOperation::MATCH);
    record.append_cigar(1, paw::CigarOperation::INSERTION);
    record.append_cigar(4, paw::CigarOperation::MATCH);

    REQUIRE(record.is_cigar_valid());
    REQUIRE(record.cig.size() == 3);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::INSERTION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 3);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 4);

    auto ref_it = ref.cbegin() + record.cig[0].count;
    auto read_it = read.cbegin() + record.cig[0].count;
    left_align_insertion(record, 1, ref_it, read_it);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::INSERTION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 2);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 5);
  }

  {
    // test case with 3bp left aligning
    SAMRecord record;
    std::string ref{"GCAAATT"};
    std::string read{"GCAAAATT"};

    record.seq = read;

    record.append_cigar(5, paw::CigarOperation::MATCH);
    record.append_cigar(1, paw::CigarOperation::INSERTION);
    record.append_cigar(2, paw::CigarOperation::MATCH);

    REQUIRE(record.is_cigar_valid());
    REQUIRE(record.cig.size() == 3);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::INSERTION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 5);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 2);

    auto ref_it = ref.cbegin() + record.cig[0].count;
    auto read_it = read.cbegin() + record.cig[0].count;
    left_align_insertion(record, 1, ref_it, read_it);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::INSERTION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 2);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 5);
  }

  {
    // test case where left aligning is not possible
    SAMRecord record;
    std::string ref{"GCAAAAT"};
    std::string read{"GCAAAATT"};

    record.seq = read;

    record.append_cigar(6, paw::CigarOperation::MATCH);
    record.append_cigar(1, paw::CigarOperation::INSERTION);
    record.append_cigar(1, paw::CigarOperation::MATCH);

    REQUIRE(record.is_cigar_valid());
    REQUIRE(record.cig.size() == 3);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::INSERTION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 6);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 1);

    auto ref_it = ref.cbegin() + record.cig[0].count;
    auto read_it = read.cbegin() + record.cig[0].count;
    left_align_insertion(record, 1, ref_it, read_it);

    REQUIRE(record.cig[0].operation == paw::CigarOperation::MATCH);
    REQUIRE(record.cig[1].operation == paw::CigarOperation::INSERTION);
    REQUIRE(record.cig[2].operation == paw::CigarOperation::MATCH);

    REQUIRE(record.cig[0].count == 6);
    REQUIRE(record.cig[1].count == 1);
    REQUIRE(record.cig[2].count == 1);
  }
}

//! \cond TESTS
TEST_CASE("Tests for sr_cigar.cpp.", "[sam]")
{
  test_left_align_deletion();
  test_left_align_insertion();
}
//! \endcond
} // namespace weaver
