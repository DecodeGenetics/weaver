#pragma once

#include <limits> // std::numeric_limits

#include "icu.hpp"          // T_icu
#include "sr_alignment.hpp" // SRAlignment
#include "sr_seed.hpp"      // SRSeed

namespace weaver
{
class GFA;

/*!
 * @brief Selects the best alignments and calculates the MAPQ of the alignment.
 *
 * @headerfile sr_alignment_sel.hpp "weaver/sr_alignment_sel.hpp"
 *
 * @details
 * The object stores various metrics for the purposes of selecting the best alignments. This includes the indexes
 * to best/secondary best hits of the input SRAlignment and SRSeed vectors for each read of the pair, s1+s2 and
 * s1_2nd+s2_2nd and multiple different scores.
 */
class SRAlignmentSel
{
public:
  /*!
   * @name Constants
   * @{
   */

  //! Value for missing score
  int constexpr static MISSING_SCORE{std::numeric_limits<int>::lowest() / 2};

  //! Value for a missing log_sum_exp
  double constexpr static MISSING_LOG_SUM_EXP{std::numeric_limits<double>::lowest()};

  //! Value for index that has not been set
  int constexpr static MISSING_INDEX{-1};

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  // Indexes to primary and secondary hits.
  int s1{MISSING_INDEX};     //!< Read 1 index of best hit
  int s2{MISSING_INDEX};     //!< Read 2 index of best hit
  int s1_2nd{MISSING_INDEX}; //!< Read 1 index to secondary best hit. TODO use this as a tag
  int s2_2nd{MISSING_INDEX}; //!< Read 2 index to secondary best bit. TODO use this as a tag

  // Score metrics
  int best_pair_score{MISSING_SCORE};             //!< Best pair score among any two selected alignments.
  int second_best_pair_score{MISSING_SCORE};      //!< Secondary best pair score (one might be shared with best).
  int other_read1_best_read_score{MISSING_SCORE}; //!< Best read score when read 1 is not shared with best.
  int other_read2_best_read_score{MISSING_SCORE}; //!< Best read score when read 2 is not shared with best.
  double log_sum_exp{MISSING_LOG_SUM_EXP};        //!< Score sum of exponents. Used in calculating mapping quality.

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  //! Get the score difference between the best and second best matches.
  int get_score_diff() const;

  /*!
   * @brief Get mapping quality for read 1.
   *
   * @param[in] best_read1_score Read score of read 1 the primary pair.
   * @param[in] read1_length Length of read 1.
   *
   * @returns The mapping quality for read 1 in range [0, 60]
   */
  int get_mapq1(int best_read1_score, int const read1_length) const;

  /*!
   * @brief Get mapping quality for read 2.
   *
   * @param[in] best_read2_score Read score of read 2 the primary pair.
   * @param[in] read2_length Length of read 2.
   *
   * @returns The mapping quality for read 2 in range [0, 60]
   */
  int get_mapq2(int best_read2_score, int const read2_length) const;

  /*!
   * @}
   */
};

/*!
 * @brief Select the primary and secondary alignments.
 *
 * @details
 * The counts of SAM records and seeds must match for each read, i.e. sam_records1.size() == seeds1.size() and same for
 * read 2. Also, either the aligments/seeds for read 1 or read 2 should be non-empty.
 *
 * @param[in] gfa           Input graph.
 * @param[in] icu           The "I see you" index containing shortest distances in the graph between vertices.
 * @param[in] seeds1        All seed objects for read 1.
 * @param[in] sam_records1  All SAM records for read 1.
 * @param[in] seeds2        All seed objects for read 2.
 * @param[in] sam_records2  All SAM records for read 2.
 *
 * @returns Object containing indexes to the primary and secondary alignments and related scores.
 *
 * @see SRAlignmentSel for more details on the returned object.
 */
SRAlignmentSel select_alignment(GFA const & gfa,
                                T_icu const & icu,
                                std::vector<SRSeed> const & seeds1,
                                std::vector<SAMRecord> const & sam_records1,
                                std::vector<SRSeed> const & seeds2,
                                std::vector<SAMRecord> const & sam_records2);

} // namespace weaver
