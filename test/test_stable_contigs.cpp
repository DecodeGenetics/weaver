#include <limits>
#include <string>
#include <vector>

#include <parallel_hashmap/phmap.h>

#include <weaver/stable_contigs.hpp>

#include "test.hpp"

namespace weaver::test
{
/*! \brief Tests for StableContigs::add_contig().
 *
 * \details
 * Also tests StableContigs::get_contig_length(int snid) , StableContigs::get_contig_name(int snid)
 * StableContigs::get_contig_length(std::string const & contig_name)
 */
static void test_stable_contigs_add_contig()
{
  StableContigs scontigs(/*bucket_size=*/10'000);

  REQUIRE(scontigs.get_num_contigs() == 0); // empty contigs vector
  REQUIRE(scontigs.name2sfa_idx.size() == 0);
  REQUIRE(scontigs.snid2sfa_idx.size() == 0);

  scontigs.add_contig(/*name=*/"test1", /*snid=*/0, /*min=*/0, /*max=*/10'000, /*rank=*/0);

  REQUIRE(scontigs.get_num_contigs() == 1);
  REQUIRE(scontigs.name2sfa_idx.size() == 1);
  REQUIRE(scontigs.snid2sfa_idx.size() == 1);
  REQUIRE(scontigs.get_contig_reference("test1").max == 10'000);
  REQUIRE(scontigs.contigs[0].max == 10'000);

  scontigs.add_contig(/*name=*/"test2", /*snid=*/0, /*min=*/0, /*max=*/10'001, /*rank=*/0);

  REQUIRE(scontigs.get_num_contigs() == 2);
  REQUIRE(scontigs.name2sfa_idx.size() == 2);
  REQUIRE(scontigs.snid2sfa_idx.size() == 1); // not 2 because new contig has same snid as before
  REQUIRE(scontigs.get_contig_reference("test1").max == 10'000);
  REQUIRE(scontigs.contigs[0].max == 10'000);
  REQUIRE(scontigs.get_contig_reference("test2").max == 10'001);
  REQUIRE(scontigs.contigs[1].max == 10'001);
}

//! Tests for StableContigs::get_num_buckets() and StableContigs::get_bucket_index(uint64_t order)
static void test_stable_contigs_buckets()
{
  StableContigs scontigs(/*bucket_size=*/10'000);

  REQUIRE(scontigs.bucket_size == 10'000);
  REQUIRE(scontigs.get_num_buckets() == 0);
  REQUIRE(scontigs.get_bucket_index(0) == 0);
  REQUIRE(scontigs.get_bucket_index(std::numeric_limits<uint64_t>::max()) == 0);

  scontigs.add_contig(/*name=*/"test1", /*snid=*/0, /*min=*/0, /*max=*/10'000, /*rank=*/0);

  REQUIRE(scontigs.bucket_size == 10'000);
  REQUIRE(scontigs.get_num_buckets() == 1);
  REQUIRE(scontigs.get_bucket_index(0) == 0);
  REQUIRE(scontigs.get_bucket_index(std::numeric_limits<uint64_t>::max()) == 1);

  scontigs.add_contig(/*name=*/"test2", /*snid=*/0, /*min=*/0, /*max=*/10'001, /*rank=*/0);

  REQUIRE(scontigs.bucket_size == 10'000);
  REQUIRE(scontigs.get_num_buckets() == 3);
  REQUIRE(scontigs.get_bucket_index(0) == 0);
  REQUIRE(scontigs.get_bucket_index(1ull << 32) == 1);
  REQUIRE(scontigs.get_bucket_index(std::numeric_limits<uint64_t>::max()) == 3);

  scontigs.add_contig(/*name=*/"test3", /*snid=*/0, /*min=*/0, /*max=*/20'001, /*rank=*/0);

  REQUIRE(scontigs.bucket_size == 10'000);
  REQUIRE(scontigs.get_num_buckets() == 6);
  REQUIRE(scontigs.get_bucket_index(0) == 0);
  REQUIRE(scontigs.get_bucket_index(1ull << 32) == 1);
  REQUIRE(scontigs.get_bucket_index(2ull << 32) == 3);
  REQUIRE(scontigs.get_bucket_index(2ull << 32 | 9'999ull) == 3);
  REQUIRE(scontigs.get_bucket_index(2ull << 32 | 10'000ull) == 4);
  REQUIRE(scontigs.get_bucket_index(std::numeric_limits<uint64_t>::max()) == 6);

  scontigs.add_contig(/*name=*/"test4", /*snid=*/0, /*min=*/0, /*max=*/30'000, /*rank=*/0);

  REQUIRE(scontigs.bucket_size == 10'000);
  REQUIRE(scontigs.get_num_buckets() == 9);
  REQUIRE(scontigs.get_bucket_index(0) == 0);
  REQUIRE(scontigs.get_bucket_index(1ull << 32) == 1);
  REQUIRE(scontigs.get_bucket_index(2ull << 32) == 3);
  REQUIRE(scontigs.get_bucket_index(2ull << 32 | 9'999ull) == 3);
  REQUIRE(scontigs.get_bucket_index(2ull << 32 | 10'000ull) == 4);
  REQUIRE(scontigs.get_bucket_index(3ull << 32 | 10'000ull) == 7);
  REQUIRE(scontigs.get_bucket_index(3ull << 32 | 29'999ull) == 8);
  REQUIRE(scontigs.get_bucket_index(std::numeric_limits<uint64_t>::max()) == 9);
}

//! \cond TESTS
TEST_CASE("Tests for stable_contigs.cpp", "[]")
{
  test_stable_contigs_add_contig();
  test_stable_contigs_buckets();
}
//! \endcond
} // namespace weaver::test
