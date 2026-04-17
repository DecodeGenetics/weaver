#pragma once
/*!
 * @file sr_alignment.hpp
 *
 * @brief Defines the SRAlignment structure and several functions related to short read alignments.
 */

#include <gfa.h>       // gfa_arc_t
#include <limits>      // std::numeric_limits
#include <memory>      // std::unique_ptr
#include <string_view> // std::string_view
#include <tuple>       // std::tuple
#include <vector>      // std::vector

#include <paw/align/alignment_results.hpp> // paw::AlignmentResults

#include "icu.hpp"

namespace weaver
{
class GFA;
class SRSeed;
class SAMRecord;

/*!
 * @brief Stores short read begin and end extension alignments of a seed.
 *
 * @headerfile sr_alignment.hpp "weaver/sr_alignment.hpp"
 *
 * @see extend_and_align() generates this object from a seed.
 */
class SRAlignment
{
public:
  /*!
   * @name Constants
   * @{
   */
  int constexpr static MISSING_SCORE{std::numeric_limits<int>::lowest() / 2};

  /*!
   * @}
   * @name Public instance variables
   * @{
   */
  std::unique_ptr<paw::AlignmentResults> begin_alignment_ext{}; //!< Begin alignment extension results
  std::unique_ptr<paw::AlignmentResults> end_alignment_ext{};   //!< End alignment extension results
  std::vector<gfa_arc_t const *> begin_arcs;                    //!< Arcs in begin alignment extension
  std::vector<gfa_arc_t const *> end_arcs;                      //!< Arcs in end alignment extension

  /*!
   * @brief Alignment score
   *
   * @details
   * The total score of the seed and extensions. Set as the smallest integer if it is missing or unknown.
   */
  int score{MISSING_SCORE};

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  //! Clear all instance variables.
  void clear();

  //! Clear begin alignment only.
  void clear_begin_alignment_extension();

  //! Clear end alignment only.
  void clear_end_alignment_extension();

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  /*! @brief Get the total score of both extensions.
   *
   * @details
   * Does not include the score of the seed, that is stored in \a SRSeed . The instance variable \a score includes
   * the score of the seed.
   *
   * @returns Total score of both extensions
   */
  int get_ext_score() const;

  /*!
   * @}
   */
};

/*!
 * @brief Extend and align a short read seed.
 *
 * @param[in] seed Seed to extend and align.
 * @param[in] query_seq View of the query sequence.
 * @param[in] query_seq_rev View of the reverse complement of the query sequence.
 *
 * @returns SRAlignment object containing the alignemnt results and arcs crossed in the graph.
 */
SRAlignment extend_and_align(SRSeed & seed, std::string_view query_seq, std::string_view query_seq_rev);

/*!
 * @brief Creates a string containing alignment results.
 *
 * @details
 * For debugging purposes.
 *
 * @param[in] ar The alignment results to make a string of.
 *
 * @returns A string with the alignment results.
 */
std::string alignment_results_to_string(paw::AlignmentResults const & ar);

/*!
 * @brief Orders seeds and their extensions such that larger scores have priority.
 *
 * @param[in] s Seed to check score for.
 * @param[in] sa Extension alignment to check score for.
 * @param[in] o Other seed to check.
 * @param[in] oa Other Extension alignment to check.
 *
 * @returns True iff s+sa has a greater score than o+oa.
 */
bool seed_alignment_order_gt(SRSeed const & s, SRAlignment const & sa, SRSeed const & o, SRAlignment const & oa);

/*!
 * @brief Get an estimation of the seed alignment score.
 *
 * @details
 * The score includes both the seed and extension score. The extension score are exact, the seed score is estimated
 * with SRSeed::get_est_score() .
 *
 * @param[in] seed Seed to get the score of
 * @param[in] sa Seed extension alignment.
 *
 * @returns The estimated score.
 */
int get_est_seed_alignment_score(SRSeed const & seed, SRAlignment const & sa);

void remove_seeds_with_no_score(std::vector<SRSeed> & seeds, std::vector<SRAlignment> & aln);

void remove_duplicates(GFA const & gfa, T_icu const & icu, std::vector<SRSeed> & seeds, std::vector<SRAlignment> & aln);

void remove_duplicates(GFA const & gfa,
                       T_icu const & icu,
                       std::vector<SRSeed> & seeds,
                       std::vector<SRAlignment> & aln,
                       std::vector<SAMRecord> & recs);

} // namespace weaver
