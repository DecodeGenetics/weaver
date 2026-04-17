#pragma once

#include <limits>
#include <string>
#include <vector>

#include "gfa.hpp"
#include "icu.hpp"
#include "sr_seed.hpp"

namespace weaver
{
/*!
 * @brief Structure for storing indexes to a pair of SRSeed and their score.
 *
 * @headerfile sr_seed_pair.hpp "weaver/sr_seed_pair.hpp"
 *
 * @details
 * Seed pairs with higher scores are thought to be more promising mappings to the graph. The object does not contain
 * the seeds of each read per se, only indices to them in external vectors.
 */
class SRSeedPair
{
public:
  /*!
   * @name Constants
   * @{
   */
  //! Minimum score constant.
  int static constexpr MIN_SCORE{std::numeric_limits<int>::min()};

  //! Maximum number of seed pair returned from get_best_seed_pairs() .
  int static constexpr MAX_CHAINS{18};

  /*!
   * @brief The maximum number of seed pairs to check, per read in the pair.
   *
   * @details
   * Note however, that if the next seed also has the same score, that one is also checked.
   * This is necessary to remove any systematic selection of equal mapping locations.
   */
  int static constexpr MAX_CHAINS_CHECKED{1024};

  /*!
   * @brief Hard limit on the number of seed to check, per read in the pair.
   *
   * @details
   * This limit should not be reached in practice
   */
  int static constexpr MAX_CHAINS_CHECKED_HARD_LIMIT{4096};

  //! Ignore seed pairs that have a score < best_seed_pair_score - MAX_CHAIN_SCORE_DIFF
  int static constexpr MAX_CHAIN_SCORE_DIFF{70};

  /*!
   * @}
   *
   * @name Public instance variables.
   * @{
   */
  int seed_i1{-1};      //!< seed index of read 1.
  int seed_i2{-1};      //!< seed index of read 2.
  int score{MIN_SCORE}; //!< score of the chain.

  /*!
   * @}
   *
   * @name Constructors and destructor
   * @{
   */

  /*!
   * @brief Constructs an new seed pair given their indices and score.
   *
   * @param[in] i1 is seed index of read 1.
   * @param[in] i2 is seed index of read 2.
   * @param[in] sc the score of the chain.
   */
  SRSeedPair(int const i1, int const i2, int const sc) : seed_i1(i1), seed_i2(i2), score(sc)
  {
  }

  SRSeedPair() = default;  //!< Explicit defaulted empty constructor.
  ~SRSeedPair() = default; //!< Explicit defaulted destructor.

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  //! Represent the short read chain as a string, for debugging purposes.
  std::string to_string() const;

  /*!
   * @}
   */
};

/*!
 * @brief Get the best seed pairs and their score.
 *
 * @param[in] gfa Graph which the seeds were created from.
 * @param[in] icu I see you index, used for estimating the minimum distance between two seeds.
 * @param[in] seeds1 Seeds from read 1.
 * @param[in] seeds2 Seeds from read 2.
 *
 * @returns The best seed pairs and their score.
 *
 * @see filter_sr_seeds() is normally called afterwards to filter the reads based on the seed pairings.
 */
std::vector<SRSeedPair> get_best_seed_pairs(GFA const & gfa,
                                            T_icu const & icu,
                                            std::vector<SRSeed> const & seeds1,
                                            std::vector<SRSeed> const & seeds2);

/*!
 * @brief Use the chains to remove seeds which are not part of any chain.
 *
 * @param[in, out] seeds1 Seeds from read 1.
 * @param[in, out] seeds2 Seeds from read 2.
 * @param[in] chains Chains to check if the seeds are in.
 *
 * @see get_best_seed_pairs()
 */
void filter_sr_seeds(std::vector<SRSeed> & seeds1,
                     std::vector<SRSeed> & seeds2,
                     std::vector<SRSeedPair> const & chains);

/*!
 * @brief Finds the best seed pairs and removes seeds that are only in pair with very poor score.
 *
 * @param[in] gfa Graph which the seeds were created from.
 * @param[in] icu I see you index, used for estimating the minimum distance between two seeds.
 * @param[in,out] seeds1 Seeds from read 1.
 * @param[in,out] seeds2 Seeds from read 2.
 *
 * @see get_best_seed_pairs()
 * @see filter_sr_seeds()
 */
void get_sr_chains_and_filter(GFA const & gfa,
                              T_icu const & icu,
                              std::vector<SRSeed> & seeds1,
                              std::vector<SRSeed> & seeds2);

} // namespace weaver
