#pragma once

#include <cstdint> // std::uint64_t
#include <string>  // std::string
#include <utility> // std::pair
#include <vector>  // std::vector

#include "edit.hpp"    // EditCalls
#include "gfa.hpp"     // weaver::GFA
#include "hashmap.hpp" // phmap::flat_hash_map

namespace weaver
{
/*!
 * @brief Minimizer index data structure.
 *
 * @headerfile mmi.hpp "weaver/mmi.hpp"
 *
 * @details
 * The minimizer data stores both a \a map and a vector of \a values .
 *
 * In the map, keys are the minimizer/sketch and the values are its locations in the graph.
 * Unique and non-unique locations are handled differently. Unique locations have the
 * most significant bit set as 0, and the location is in the remaining 63 bits.
 * If there are multiple locations then the most significant bit is 1, the next 35 bits are
 * the index in \a values and the last 28 bits are the number of locations (2 or more). The locations
 * are always stored sequently in \a values so these bits form a compact span into the vector.
 *
 * @see make_mmi_index() for indexing a GFA and returning an MMI index.
 */
class MMI
{
public:
  /*!
   * @name Internal types
   * @{
   */
  //! Type that stores the site information of minimizer hits.
  using T_mmi_values = std::vector<uint64_t>;

  //! Minimizer map where multiple values are allowed for each key.
  using T_multi_value_map = phmap::flat_hash_map<uint64_t, std::vector<uint64_t>>;

  //! Minimizer map where only unique values are allowed
  using T_uniq_value_map = phmap::flat_hash_map<uint64_t, uint64_t>;

  //! Type for keeping a unique set of minimizer key and values.
  using T_unique_set = phmap::flat_hash_set<std::pair<uint64_t, uint64_t>>;

  //! Type for storing the haplotype data
  using T_haplotypes = std::vector<std::vector<EditCalls>>;

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  T_uniq_value_map map;    //!< Minimizer keys and unique values or spans to non-unique values.
  T_mmi_values values;     //!< Non-unique values.
  T_haplotypes haplotypes; //!< Haplotype data.

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  /*!
   * @brief Get all values of a key.
   *
   * @details
   * The values are not returned, only appended to a referenced vector.
   *
   * @param[out] returned_values The results are appended back of this vector.
   * @param[in] key minimizer key to query.
   */
  void get_values(std::vector<uint64_t> & returned_values, uint64_t key) const;

  //! Get the number of keys/unique kmers stored in the index.
  int64_t num_keys() const;

  /*!
   * @}
   *
   * @name Operators
   * @{
   */

  //! Check if index is equal to other index \a o.
  bool operator==(MMI const & o) const;

  //! Check if index is not equal to other index \a o.
  bool operator!=(MMI const & o) const;

  /*!
   * @}
   */
};

/*!
 * @brief Checks if the mmi_uniq data structure contains valid values.
 *
 * @details
 * For debugging purposes. Call this function before inserting values into the unique and multi maps.
 *
 * @param[in] mmi_uniq The MMI unique index map.
 * @param[in] num_segments Number of segments in the indexed graph.
 *
 * @returns False iff problems are found.
 */
bool is_mmi_uniq_map_valid(MMI::T_unique_set const & mmi_uniq, int const num_segments);

} // namespace weaver
