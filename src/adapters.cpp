#include <array>
#include <utility>

#include "hashmap.hpp"
#include "logging.hpp"
#include "read_sketch.hpp"
#include "sketch.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp" // sketch_value_pos, sketch_value_strand

namespace
{
uint32_t get_pos_strand(uint64_t val1, uint64_t val2)
{
  using namespace weaver;

  return static_cast<uint32_t>((sketch_value_pos(val1) + sketch_value_pos(val2)) << 1) | //
         static_cast<uint32_t>(sketch_value_strand(val1) != sketch_value_pos(val2));
}

} // namespace

namespace weaver
{
// 0 means the adapther was found or that neither read had a hit, 1 is only read1 found a hit, 2 is only read2 found a
// hit, 3 is no hit
#ifndef NDEBUG
int check_minimizers_for_adapters(std::vector<ReadSketch> & read_sketches1,
                                  std::vector<ReadSketch> & read_sketches2,
                                  std::string_view read_name,
                                  int i1,
                                  int i2)
#else
int check_minimizers_for_adapters(std::vector<ReadSketch> & read_sketches1,
                                  std::vector<ReadSketch> & read_sketches2,
                                  std::string_view /*read_name*/,
                                  int i1,
                                  int i2)
#endif // NDEBUG
{
  ReadSketch const & min1 = read_sketches2[i1];
  ReadSketch const & min2 = read_sketches1[i2];

  if (min1.first == min2.first)
    return 0;

  auto it1 = std::find(read_sketches1.cbegin(), read_sketches1.cend(), min1.first);
  auto it2 = std::find(read_sketches2.cbegin(), read_sketches2.cend(), min2.first);

  // try to find with the first minimizer in each
  while (it1 != read_sketches1.cend() && it2 != read_sketches2.cend())
  {
    auto pos_strand1 = get_pos_strand(it1->second, min1.second);
    auto pos_strand2 = get_pos_strand(it2->second, min2.second);

    bool const is_in_sync = pos_strand1 == pos_strand2;

    // found in both reads
    print_debug(_HERE_,
                " name = ",
                read_name,
                " adt_minimizers1,2 = ",
                std::distance(it1 + 1, read_sketches1.cend()),
                " of ",
                read_sketches2.size() - 1,
                " and ",
                std::distance(it2 + 1, read_sketches2.cend()),
                " of ",
                read_sketches1.size() - 1,
                " in_sync = ",
                is_in_sync,
                " pos_strand 1, 2 = ",
                pos_strand1,
                ", ",
                pos_strand2);

    if (is_in_sync)
    {
      // resize to ignore adapter minimizers
      read_sketches1.resize(1 + std::distance(read_sketches1.cbegin(), it1));
      read_sketches2.resize(1 + std::distance(read_sketches2.cbegin(), it2));
      return 0;
    }

    // try to rescue the minimizer in case of duplicates
    if (pos_strand1 < pos_strand2)
      it1 = std::find(++it1, read_sketches1.cend(), min1.first);
    else
      it2 = std::find(++it2, read_sketches2.cend(), min2.first);
  }

  if (it1 != read_sketches1.cend())
    return 1;

  if (it2 != read_sketches2.cend())
    return 2;

  return 0;
}

void remove_adapter_sketches(std::vector<ReadSketch> & read_sketches1,
                             std::vector<ReadSketch> & read_sketches2,
                             std::string_view read_name)
{
  if (read_sketches1.size() <= 1 || read_sketches2.size() <= 1)
    return;

  // int ret =
  check_minimizers_for_adapters(read_sketches1, read_sketches2, read_name, 0, 0);

  // if (ret == 1 && read_sketches1.size() > 3)
  //{
  //  print_debug(_HERE_, " read_sketches1[0] might have an error, try read_sketches1[1] at read = ", read_name);
  //  check_minimizers_for_adapters(read_sketches1, read_sketches2, read_name, 0, 2);
  //}
  // else if (ret == 2 && read_sketches2.size() > 3)
  //{
  //  print_debug(_HERE_, " read_sketches2[0] might have an error, try read_sketches2[1] at read = ", read_name);
  //  check_minimizers_for_adapters(read_sketches1, read_sketches2, read_name, 2, 0);
  //}
}

} // namespace weaver
