#pragma once
/*!
 * @file sr_chain.hpp
 * @brief Interface of the SRChain class.
 */

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "gfa.hpp"
#include "icu.hpp"

namespace weaver
{
class SRSeed;
class SRAlignment;

/*!
 * @brief Short read chain stores indices to seeds from each read and their scores.
 *
 * @headerfile sr_chain.hpp "weaver/sr_chain.hpp"
 *
 * @details
 * The SRChain class is both used to keep track of scores while selecting best seeds to align and while selecting the
 * best alignments.
 */
class SRChain
{
public:
  //! Minimum score constant.
  int static constexpr MIN_SCORE{std::numeric_limits<int>::lowest() / 2};

  int score{MIN_SCORE};                                 //!< score of the chain
  int score1{MIN_SCORE};                                //!< score for read 1
  int score2{MIN_SCORE};                                //!< score for read 2
  int score_secondary1{MIN_SCORE};                      //!< next best score of the chain
  int score_secondary2{MIN_SCORE};                      //!< next best score of the chain
  std::vector<std::pair<int, int>> hits{};              //!< List of primary hits.
  std::vector<std::pair<int, int>> seed_secondary1_i{}; //!< List of secondary hits for read 1
  std::vector<std::pair<int, int>> seed_secondary2_i{}; //!< List of secondary hits for read 2

  SRChain() = default; //!< Explicit defaulted empty constructor

  //! Check if this score is best/secondary score
  void check_score(GFA const & gfa,
                   T_icu const & icu,
                   int new_score, //
                   int i1,
                   int i2,
                   int new_score1,
                   int new_score2,
                   std::vector<SRSeed> const & seeds1,
                   std::vector<SRSeed> const & seeds2);

  //! Represent the short read chain as a string, for debugging purposes.
  std::string to_string() const;

  //! Mapping quality for read 1 based on chains
  int get_mapping_quality1() const;

  //! Mapping quality for read 2 based on chains
  int get_mapping_quality2() const;

  //! Get the indexes for the primary and secondary hit. Missing value is -1.
  void get_primary_and_secondary_index(std::pair<int, int> & primary_index,
                                       std::pair<int, int> & secondary_index,
                                       std::vector<SRSeed> const & s1,
                                       std::vector<SRSeed> const & s2);

  //! Get the indexes for the primary hit. Missing value is -1.
  std::pair<int, int> get_primary_index(std::vector<SRSeed> const & s1, std::vector<SRSeed> const & s2) const;
};

SRChain get_chains(GFA const & gfa,
                   T_icu const & icu,
                   std::vector<SRSeed> const & seeds1,
                   std::vector<SRAlignment> const & alignments1,
                   std::vector<SRSeed> & seeds2,
                   std::vector<SRAlignment> const & alignments2);

} // namespace weaver
