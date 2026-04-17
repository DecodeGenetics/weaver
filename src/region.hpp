#pragma once

#include <limits>
#include <string>
#include <vector>

namespace weaver
{
class GFA;
class SRSeed;

//! Class that parses a region string
class Region
{
public:
  Region() = default;
  Region(Region const &) = default;
  Region(Region &&) = default;
  Region & operator=(Region const &) = default;
  Region & operator=(Region &&) = default;
  ~Region() = default;

  /*!
   * @brief Construction a region from a \a region string
   *
   * @details
   * The format of the region string should be one of the following:
   *
   * 1) chrN<br>
   * 2) chrN:A, is same as chrN:A-A<br>
   * 3) chrN:A-B, where B>=A<br>
   *
   * A and B are 1-based inclusive. Internally, these will be changed to 0-based and B will be non-inclusive.
   */
  explicit Region(std::string const & region);

  std::string chr{};                        //!< Chromosome name
  int begin{0};                             //!< Begin position, 0-based inclusive
  int end{std::numeric_limits<int>::max()}; //!< End position, 0-based non-inclusive

  //! Check if the region is logical, i.e. not with a missing chromosome or with an end smaller than begin
  bool check() const;

  //! Size of the region
  int size() const;
};

//! Parse the target region command line option
void parse_target_region_option(GFA const & gfa);

//! Check if seed is in the target region
bool is_seed_in_region(GFA const & gfa, SRSeed const & seed);

//! Check if all seeds are outside the target region
bool are_all_seeds_outside_region(GFA const & gfa, std::vector<SRSeed> const & seeds);

} // namespace weaver
