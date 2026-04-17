#include "read_sketch.hpp"

#include <cassert>
#include <cstdint>
#include <string_view>
#include <vector>

#include <weaver/constants.hpp> // KMIN

#include "logging.hpp"
#include "mmi.hpp"
#include "sketch.hpp"
#include "sketch_cache.hpp"

namespace weaver
{
std::vector<ReadSketch> get_read_sketches(std::string_view str, int const rid, MMI const & mmi, int k, int w)
{
  assert(rid >= 0);
  std::vector<ReadSketch> read_sketches;
  SketchCache<ReadSketch> cache;

  bool constexpr is_reading_forward{true};

  // If no k is set, then use default
  if (k <= 0)
    k = K_DEFAULT;

  // If no w is set, then use default
  if (w <= 0)
    w = W_DEFAULT;

  // Sketch the query string
  sketch<ReadSketch>(read_sketches, cache, str, w, k, static_cast<uint32_t>(rid), is_reading_forward, -1);
  push_last_sketch<ReadSketch>(read_sketches, cache, -1, k);

  // Find the sketches which are not identical to their neighbouring sketch
  if (read_sketches.size() > 0)
  {
    std::vector<ReadSketch> uniq_read_sketches;
    uniq_read_sketches.push_back(read_sketches[0]);

    if (read_sketches.size() > 1)
    {
      int skipped{0};
      int constexpr MAX_SKIPS{6};

      for (int i{1}; i < static_cast<int>(read_sketches.size() - 1); ++i)
      {
        if (skipped >= MAX_SKIPS ||                                 //
            read_sketches[i - 1].first != read_sketches[i].first || //
            read_sketches[i + 1].first != read_sketches[i].first)
        {
          uniq_read_sketches.push_back(read_sketches[i]);
          skipped = 0;
        }
        else
        {
          ++skipped;
        }
      }

      uniq_read_sketches.push_back(read_sketches[read_sketches.size() - 1]);
    }

    // get mmi results
    for (ReadSketch & read_sketch : uniq_read_sketches)
      mmi.get_values(/*returned values=*/read_sketch.mmi_results, /*key=*/read_sketch.first);

    return uniq_read_sketches;
  }
  else
  {
    return read_sketches;
  }
}

} // namespace weaver
