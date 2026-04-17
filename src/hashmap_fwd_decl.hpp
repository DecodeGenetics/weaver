#pragma once

#include <cstdint> // uint64_t

#include <parallel_hashmap/phmap_fwd_decl.h>

namespace weaver
{
// type definition for maps
template <typename T1, typename T2>
using Tmap = phmap::flat_hash_map<T1, T2>;

// type definition for sets
template <typename T>
using Tset = phmap::flat_hash_set<T>;

// type definition for sets of uint64_t
using Tset_u64 = phmap::flat_hash_set<uint64_t>;

} // namespace weaver
