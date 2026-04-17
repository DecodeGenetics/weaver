/*!
 * @file index_io.cpp
 * @brief Implements the function for the index I/O operations.
 */
#include "index_io.hpp"

#include <fstream>
#include <string_view>

#include "filesystem.hpp"
#include "gfa.hpp"
#include "gfa_location.hpp"
#include "hashmap.hpp"
#include "icu.hpp"
#include "index.hpp"
#include "logging.hpp"
#include "mmi.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>

namespace
{
inline std::size_t get_file_size(weaver::filesystem::path const & path)
{
  return weaver::filesystem::file_size(weaver::filesystem::canonical(path));
}

} // namespace

namespace weaver
{
bool deserialize_index(filesystem::path const & graph_path, MMI & mmi, T_icu & icu, int & k, int & w)
{
  filesystem::path index_path(graph_path);
  index_path += ".wmi";

  if (!filesystem::exists(index_path))
  {
    print_info("No index file could be found.");
    return false;
  }

  std::ifstream is(index_path, std::ios::binary);
  cereal::BinaryInputArchive archive(is);

  // get k
  {
    int new_k{0};
    archive(new_k);
    print_info("Using k=", new_k);

    if (k > 0 && new_k != k)
    {
      print_warning("Skipping deserialization. Index kmer size is=", new_k, " but not the expected=", k, ")");
      return false;
    }

    k = new_k;
  }

  // get w
  {
    int new_w{0};
    archive(new_w);
    print_info("Using w=", new_w);

    if (w > 0 && new_w != w)
    {
      print_warning("Skipping deserialization. Index kmer size is=", new_w, " but not the expected=", w, ")");
      return false;
    }

    w = new_w;
  }

  std::size_t gfa_size = 0;
  archive(gfa_size);

  auto const current_size = get_file_size(graph_path);

  if (current_size != gfa_size)
  {
    print_error("GFA file has different size than at the time of index creation.");
    print_error("  current gfa size != previous gfa size (", current_size, " != ", gfa_size, ")");
    print_error("This may mean you have provided the wrong .wmi file or that the file is corrupted.");
    print_error("Please recreate the index with 'weaver index'.");
    std::exit(1);
  }

  if (filesystem::last_write_time(graph_path) > filesystem::last_write_time(index_path))
  {
    print_warning("GFA is newer than its index; this likely means that the index is out-of-date. ",
                  "Consider recreating the index with 'weaver index'.");
  }

  print_info("Deserializing minimizer index map...");
  std::size_t map_size;
  archive(map_size);
  mmi.map.reserve(map_size);
  archive(mmi.map);

  if (mmi.map.size() != map_size)
  {
    print_warning("MMI map had unexpected size of ", mmi.map.size(), " but ", map_size, " was expected.");
    return false;
  }

  print_info("Done. mmi.map.size()=", map_size);

  print_info("Deserializing minimizer index non-unique values...");
  archive(mmi.values);
  print_info("Done. mmi.values.size()=", mmi.values.size());

  print_info("Deserializing haplotype data...");
  archive(mmi.haplotypes);
  print_info("Done. Number of haplotype contigs=", mmi.haplotypes.size());

  print_info("Deserializing ICU index...");
  archive(icu);
  print_info("ICU index ready. icu.size()=", icu.size());
  return true;
}

void deserialize_or_make_index(filesystem::path const & graph_path, //
                               GFA const & gfa,
                               MMI & mmi,
                               T_icu & icu,
                               std::string const & vcf_fn,
                               int & k,
                               int & w)
{
  filesystem::path index_path(graph_path);
  index_path += ".wmi";

  if (filesystem::exists(index_path))
  {
    int old_k = k;
    int old_w = w;
    bool success = deserialize_index(graph_path, mmi, icu, k, w);

    if (!success)
      make_index(gfa, mmi, icu, vcf_fn, old_k, old_w);
  }
  else
  {
    print_info("No prebuilt index found.");
    make_index(gfa, mmi, icu, vcf_fn, k, w);
  }
}

void serialize_index(filesystem::path const & graph_path, MMI const & mmi, T_icu const & icu, int k, int w)
{
  filesystem::path index_path(graph_path);
  index_path += ".wmi";

  if (filesystem::exists(index_path))
    print_warning("Index file exists, overwriting.");

  std::ofstream os(index_path, std::ios::binary);
  cereal::BinaryOutputArchive archive(os);
  archive(k);
  archive(w);
  archive(get_file_size(graph_path));
  archive(mmi.map.size());
  archive(mmi.map);
  archive(mmi.values);
  archive(mmi.haplotypes);
  archive(icu);
}

int print_index_stats(filesystem::path const & index_path)
{
  if (!filesystem::exists(index_path))
  {
    print_error("No index file could be found at path '", index_path, "'");
    return 1;
  }

  std::ifstream is(index_path, std::ios::binary);
  cereal::BinaryInputArchive archive(is);

  std::cout << "Index stats:\n";

  // get k
  int k{0};
  archive(k);
  std::cout << "  k = " << k << '\n';

  // get w
  int w{0};
  archive(w);
  std::cout << "  w = " << w << '\n';

  // get gfa_size
  std::size_t gfa_size{0};
  archive(gfa_size);
  std::cout << "  graph size = " << gfa_size << " bytes\n";

  std::size_t map_size;
  archive(map_size);
  std::cout << "  map size = " << map_size << " keys" << std::endl;
  return 0;
}

} // namespace weaver
