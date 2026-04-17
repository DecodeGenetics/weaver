#include <vector>

#include "hashmap_fwd_decl.hpp" // weaver::Tset

namespace weaver
{
class ReadSketch;

//! Remove adapter sketches from read sketches
void remove_adapter_sketches(std::vector<ReadSketch> & read_sketches1,
                             std::vector<ReadSketch> & read_sketches2,
                             std::string_view read_name);

} // namespace weaver
