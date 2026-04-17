#include "sketch_value.hpp"

#include <cstdint>

#include "gfa.hpp"
#include "logging.hpp"

namespace weaver
{
bool is_sketch_value_valid(GFA const & gfa, uint64_t const sketch_value)
{
  bool is_valid{true};
  int rid{sketch_value_rid(sketch_value)};

  if (rid >= gfa.get_num_segments())
  {
    print_warning(_HERE_, " sketch contig ID too large ", rid, " >= ", gfa.get_num_segments());
    is_valid = false;
  }
  else
  {
    gfa_seg_t const & segment = gfa.get_segment(rid);
    int pos{sketch_value_pos(sketch_value)};

    if (pos >= segment.len)
    {
      print_warning(_HERE_, " sketch pos too large ", pos, " >= ", segment.len);
      is_valid = false;
    }
  }

  return is_valid;
}

} // namespace weaver
