#include <iostream>
#include <string>

#include <parallel_hashmap/phmap.h> // T_icu

#include <weaver/constants.hpp>
#include <weaver/filesystem.hpp>
#include <weaver/gfa.hpp>
#include <weaver/icu.hpp>
#include <weaver/index.hpp>
#include <weaver/index_io.hpp>
#include <weaver/logging.hpp>
#include <weaver/sketch_to_string.hpp>

#include <catch2/catch.hpp>

namespace weaver::test
{
//! Tests for index serialization and deserialization with a VCF
static void test_index_vcf()
{
  std::string graph_path = static_cast<std::string>(weaver_source_dir) + std::string("/test/data/test_small.gfa");
  std::string vcf_path = static_cast<std::string>(weaver_source_dir) + std::string("/test/data/test_small.vcf.gz");
  filesystem::path index_path = graph_path;
  index_path += ".wmi";
  int k{21};
  int w{10};

  if (filesystem::exists(index_path))
    filesystem::remove(index_path);

  REQUIRE(!filesystem::exists(index_path));

  GFA gfa(graph_path);
  MMI mmi1;
  T_icu icu1;
  GFA gfa_novcf(graph_path);
  MMI mmi1_novcf;
  T_icu icu1_novcf;
  make_index(gfa, mmi1, icu1, vcf_path, k, w);
  make_index(gfa_novcf, mmi1_novcf, icu1_novcf, "", k, w);

  REQUIRE(mmi1.num_keys() == 53);
  REQUIRE(mmi1_novcf.num_keys() == 47);
  REQUIRE(icu1.size() == 28);

  // all keys of novcf should be in the other
  for (auto const & key_val : mmi1_novcf.map)
  {
    std::vector<uint64_t> ret;
    std::vector<uint64_t> ret_novcf;

    mmi1.get_values(ret, key_val.first);
    mmi1_novcf.get_values(ret_novcf, key_val.first);

    for (auto const val : ret_novcf)
      REQUIRE(std::find(ret.begin(), ret.end(), val) != ret.end());
  }
}

//! Tests for index serialization and deserialization without a VCF
static void test_index_no_vcf()
{
  std::string graph_path = static_cast<std::string>(weaver_source_dir) + std::string("/test/data/test_small.gfa");
  filesystem::path index_path = graph_path;
  index_path += ".wmi";
  int k{21};
  int w{10};

  if (filesystem::exists(index_path))
    filesystem::remove(index_path);

  REQUIRE(!filesystem::exists(index_path));

  GFA gfa(graph_path);
  MMI mmi1;
  T_icu icu1;
  make_index(gfa, mmi1, icu1, "", k, w);
  REQUIRE(mmi1.num_keys() == 47);
  REQUIRE(icu1.size() == 28);

  {
    REQUIRE(!filesystem::exists(index_path));
    serialize_index(graph_path, mmi1, icu1, k, w);
    REQUIRE(filesystem::exists(index_path));
  }

  // deserialize
  {
    MMI mmi2;
    T_icu icu2;
    deserialize_index(graph_path, mmi2, icu2, k, w);
    REQUIRE(mmi1 == mmi2);
    REQUIRE(icu1 == icu2);
  }

  // deserialize or make
  {
    MMI mmi3;
    T_icu icu3;
    deserialize_or_make_index(graph_path, gfa, mmi3, icu3, "", k, w); // should deserialize
    REQUIRE(mmi1 == mmi3);
    REQUIRE(icu1 == icu3);
  }

  if (filesystem::exists(index_path))
    filesystem::remove(index_path);

  REQUIRE(!filesystem::exists(index_path));

  // deserialize or make
  {
    // Make new indices
    MMI mmi4;
    T_icu icu4;
    deserialize_or_make_index(graph_path, gfa, mmi4, icu4, "", k, w); // should make indexes
    REQUIRE(mmi1 == mmi4);
    REQUIRE(icu1 == icu4);
  }

  serialize_index(graph_path, mmi1, icu1, k, w);
  REQUIRE(filesystem::exists(index_path));

  // symlink test
  {
    MMI mmi5;
    T_icu icu5;

    filesystem::path symlink_graph_path(graph_path);
    symlink_graph_path += ".test.gfa";
    filesystem::path symlink_index_path(symlink_graph_path);
    symlink_index_path += ".wmi";

    if (filesystem::exists(symlink_graph_path))
      filesystem::remove(symlink_graph_path);

    if (filesystem::exists(symlink_index_path))
      filesystem::remove(symlink_index_path);

    filesystem::create_symlink(graph_path, symlink_graph_path);
    filesystem::create_symlink(index_path, symlink_index_path);

    deserialize_index(symlink_graph_path, mmi5, icu5, k, w);

    REQUIRE(mmi1 == mmi5);
    REQUIRE(icu1 == icu5);

    // clean up symlinks
    filesystem::remove(symlink_graph_path);
    filesystem::remove(symlink_index_path);
  }

  REQUIRE(filesystem::exists(index_path));

  // wrong k used
  {
    MMI mmi6;
    T_icu icu6;
    int k{31};
    REQUIRE(!deserialize_index(graph_path, mmi6, icu6, k, w));

    // no data is deserialized
    REQUIRE(mmi1 != mmi6);
    REQUIRE(icu1 != icu6);
    REQUIRE(mmi6.num_keys() == 0);
    REQUIRE(icu6.size() == 0);
  }
}

//! \cond TESTS
TEST_CASE("Tests for index.cpp", "[gfa][icu][index]")
{
  // order matters here, whichever index from test comes later will be kept
  test_index_vcf();
  test_index_no_vcf();
}
//! \endcond
} // namespace weaver::test
