#include "mmi.hpp"

#include <cstdint> // uint64_t
#include <gfa.h>   // gfa_seg_t
#include <vector>  // std::vector

#include "logging.hpp"      // print_error
#include "sketch_value.hpp" // sketch_value_rid

namespace weaver
{
void MMI::get_values(std::vector<uint64_t> & returned_values, uint64_t key) const
{
  auto find_it = map.find(key);

  if (find_it == map.end())
    return;

  uint64_t const value = find_it->second;

  if ((value & (1ull << 63)) == 0ull)
  {
    // unique hit
    returned_values.push_back(value);
    return;
  }

  // 1 bit (non-unique?) | 35 bits (start index) | 28 bits (amount)
  uint64_t const start = (value << 1) >> 29;
  uint64_t const amount = (value << 36) >> 36;

#ifndef NDEBUG
  if (start >= values.size())
  {
    print_error(_HERE_, " start >= values.size(). ", start, " >= ", values.size());
    assert(start < values.size());
  }
#endif

  assert((start + amount) <= values.size());
  returned_values.insert(returned_values.end(), &values[start], &values[start + amount]);
}

int64_t MMI::num_keys() const
{
  return static_cast<int64_t>(map.size());
}

bool MMI::operator==(MMI const & o) const
{
  return this->map == o.map && this->values == o.values;
}

bool MMI::operator!=(MMI const & o) const
{
  return this->map != o.map || this->values != o.values;
}

bool is_mmi_uniq_map_valid(MMI::T_unique_set const & mmi_uniq, int const num_segments)
{
  for (auto sketch : mmi_uniq)
  {
    if (sketch.second == std::numeric_limits<uint64_t>::max())
    {
      print_warning(_HERE_, " Found an invalid mmi sketch value");
      return false;
    }

    if (sketch_value_rid(sketch.second) >= num_segments || sketch_value_rid(sketch.second) < 0)
    {
      print_warning(_HERE_, " Found an invalid mmi sketch value rid");
      return false;
    }
  }

  return true;
}

} // namespace weaver
