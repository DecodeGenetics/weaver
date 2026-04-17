#pragma once

#include <cstdint>     // uint64_t
#include <string_view> // std::string_view
#include <utility>     // std::pair
#include <vector>      // std::vector

#include "mmi.hpp"

namespace weaver
{
/*!
 * @brief Class for read sketches
 *
 * @headerfile read_sketch.hpp "weaver/read_sketch.hpp"
 */
class ReadSketch
{
public:
  /*!
   * @name Public instance variables
   * @{
   */
  uint64_t first{};  //!< read key used in mmi index
  uint64_t second{}; //!< read value

  //! Results from the weaver minimizer index (MMI)
  std::vector<uint64_t> mmi_results{};

  /*!
   * @}
   *
   * @name Constructors
   * @{
   */
  ReadSketch() = default; //!< Empty constructor defaulted.

  //!
  ReadSketch(uint64_t key, uint64_t val) : first{key}, second{val}, mmi_results{}
  {
  }

  /*!
   * @}
   *
   * @name Operators
   * @{
   */
  //! Checks if two read sketch hashes are the same.
  inline bool operator==(ReadSketch const & o) const
  {
    return first == o.first; // it is enough to check only the minimizer key hash.
  }

  //! Checks if the read sketch has a given \a key .
  inline bool operator==(uint64_t const key) const
  {
    return first == key;
  }

  //! Checks if two read sketch hashes are different.
  inline bool operator!=(ReadSketch const & o) const
  {
    return first != o.first; // it is enough to check only the minimizer key hash.
  }

  //! Checks if the read sketch does not have a given \a key .
  inline bool operator!=(uint64_t const key) const
  {
    return first != key;
  }

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */
  //! Checks how many mmi results were found.
  inline int number_of_results() const
  {
    return static_cast<int>(mmi_results.size());
  }

  /*!
   * @}
   */
};

//! Find all sketches from a read sequence
std::vector<ReadSketch> get_read_sketches(std::string_view str, int const rid, MMI const & mmi, int k, int w);

} // namespace weaver
