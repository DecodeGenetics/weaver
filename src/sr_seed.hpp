#pragma once
/*!
 * @file sr_seed.hpp
 * @brief Defines the SRSeed class.
 */

#include <cstdint>     // uint64_t
#include <limits>      // std::numeric_limits
#include <string>      // std::string
#include <string_view> // std::string_view
#include <vector>      // std::vector

#include <paw/align/alignment_results.hpp> // paw::AlignmentResults

#include <weaver/constants.hpp>

#include "gfa.hpp" // GFA
#include "icu.hpp" // T_icu
// #include "mmi.hpp" // T_mmi

namespace weaver
{
class ReadSketch;

/*!
 * @brief Links the query sequence with a hit on the reference graph.
 *
 * @headerfile sr_seed.hpp "weaver/sr_seed.hpp"
 *
 * @details
 * The sketch values of {begin,end}_{ref,read}_value have inclusive interval indices:<br>
 *   [begin_ref_value, end_ref_value]   <- reference graph indices.<br>
 *   [begin_read_value, end_read_value] <- read indices.<br>
 *
 *  The motivation behing using inclusive intervals instead of more common half-open intervals, is that the seed can
 *  change strand efficiently. Half-open intervals would require a much more expensive check of next index on the graph.
 *  Half-open intervals can have length 0 but seed will never be of length 0 so it is not an issue.
 *
 *  The seed has a field named \a num_cuts which indicates how many sections of the read are inexact matches.
 *  In each cut there might several mismatches/indels.
 */
class SRSeed
{
public:
  /*!
   * @name Sketch values
   * @{
   */

  uint64_t begin_ref_value{0};  //!< begin sketch value of "ref_value" so we know where the seed starts.
  uint64_t end_ref_value{0};    //!< end sketch value for reference sequence.
  uint64_t begin_read_value{0}; //!< begin sketch value of "read_value", i.e. where the seed starts on the read.
  uint64_t end_read_value{0};   //!< end sketch value for read sequence.

private:
  /*!
   * @}
   *
   * @name Private member variables
   * @{
   */

public:
  /*!
   * @}
   *
   * @name Public member variables
   * @{
   */
  //! Mapping score of the seed. Set as the smallest integer when it has not been set.
  int score{std::numeric_limits<int>::min()};

  int length{0};                                         //!< Length of the seed.
  int begin_unaccounted_read_bases{0};                   //!< Number of unaccounted read bases while extending at begin.
  int end_unaccounted_read_bases{0};                     //!< Number of unaccounted read bases while extending at end.
  int num_cuts{0};                                       //!< Number of cuts in the seed.
  bool is_extending_cut{false};                          //!< True iff the seed has been cut.
  std::vector<gfa_arc_t const *> arcs{};                 //!< which arcs the seed crossed.
  std::unique_ptr<std::vector<paw::Cigar>> cig{nullptr}; //!< CIGAR alignment between seed and graph.

  /*!
   * @}
   *
   * @name Constructors and destructor
   * @{
   */

  // explicit c/dtors
  //! Constructor for a new seed from a reference and read minimizer values.
  SRSeed(uint64_t ref_val, uint64_t read_val) noexcept;

  //! Copy constructor for a new seed.
  SRSeed(SRSeed const & o) noexcept;

  //! Copy assignment for seeds.
  SRSeed & operator=(SRSeed const & o) noexcept;

  // implicit c/dtors
  SRSeed() noexcept = default;                        //!< Default empty constructor.
  SRSeed(SRSeed && o) noexcept = default;             //!< Move constructor for a new seed.
  SRSeed & operator=(SRSeed && o) noexcept = default; //!< Move assignment for seeds.
  ~SRSeed() noexcept = default;                       //!< Default deconstructor for seeds.

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  /*!
   * @brief Align the query sequence to the graph and calculate a score/cigar (if not all M)
   *
   * @param[in] query_seq View of the query sequence to align.
   */
  void align_to_graph(std::string_view query_seq);

  /*!
   * @brief Clear all essential instance variables.
   *
   * @details
   * Essential variables are the variable storing the read and reference values, the score, and the length of the seed.
   */
  void clear_essentials();

  //! Extend the seed end via setting a new end ref and end read values.
  void extend_end(int l, uint64_t new_end_ref_value, uint64_t new_end_read_value);

  //! Shrink the reference of the seed begin
  int shrink_reference_begin(int by);

  //! Shrink the reference on the seed end. Returns how many bases it could not shrink.
  int shrink_reference_end(int by);

  //! Shrink the beginning of the seed.
  void shrink_read_begin(int by);

  //! Shrink the ending of the seed.
  void shrink_read_end(int by);

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  //! True iff there is an indel in the cigar string. False if seed has no cigar.
  bool has_indel_in_cigar() const;

  //! True iff there both a deletion and an insertion in the cigar string. False if seed has no cigar.
  bool has_deletion_and_insertion_in_cigar() const;

  //! Return whether the read begin is reversed compared to the reference
  bool is_begin_reversed() const;

  //! Return whether the read end is reversed compared to the reference
  bool is_end_reversed() const;

  //! Check if seed is empty (no kmers)
  bool is_empty() const;

  //! Get the length of the seed on the read, i.e. determined from the begin and end read values only.
  int get_read_length() const;

  /*!
   * @brief Estimate a score for the seed.
   *
   * @details
   * The estimation is based on the number of cuts in the seed, rather than calculating an exact score by
   * checking the actual bases on the sequence and on the read. Seeds with higher scores are typically more
   * likely to be mapped correctly. This score does not take into account any extensions of the seed.
   *
   * @returns The estimated score for the seed.
   *
   * @see get_score() for the exact score if it has been calculated.
   * @see test_get_seed_estimated_score() tests functionality.
   */
  int get_est_score() const;

  /*!
   * @brief Get a score for the seed.
   *
   * @details
   * Seeds with higher scores are typically more likely to be mapped correctly. This score does not take
   * into account any extensions.
   *
   * @see get_est_score() for getting an estimated score when this exact score has not been determined.
   */
  int get_score() const;

  /*!
   * @brief Get reference sequence from graph.
   *
   * @details
   * The graph sequence that the seed covers will be extracted and returned. The seed should already have
   * a pointer to the GFA graph to use.
   */
  std::string get_ref_sequence() const;

  /*!
   * @brief Get read sequence from graph.
   *
   * @details
   * The read sequence that the seed covers will be extracted.
   */
  std::string get_read_sequence(std::string_view query_seq) const;

  /*!
   * @brief Check whether this seed sees the seed \a to.
   *
   * @details
   * "Seeing" means that the seed \a to is within certain distance threshold in the shortest walks.
   *
   * @param[in] gfa graph.
   * @param[in] icu the I see you index
   * @param[in] to the seed to check if this seeds sees it.
   * @param[in] max_dist Maximum distance to check for. Should not be larger than the size the ICU index used.
   *
   * @returns true iff this seeds sees \a to.
   *
   * @see is_within_distance(GFA const & gfa, T_icu const & icu, uint64_t from, uint64_t to, int distance).
   */
  bool do_you_see_me(GFA const & gfa, T_icu const & icu, SRSeed const & to, int const max_dist = 1500) const;

  /*!
   * @}
   *
   * @name Methods for debugging
   * @{
   */

  //! Create a string containing the contents of the seed.
  std::string to_string() const;

  //! Checks if the current state of the seed is valid.
  /*! By valid, I mean that the positions are real between begin and end we have all the arcs in the correct order. This
   *  function is useful for debugging.
   *
   *  @returns true iff the seed is in a valid state.
   */
  bool is_valid() const;

  /*!
   * @}
   */
};

/*!
 * @brief Defines ordering between two seeds \a s and \a o, such that the most promising seeds come first.
 *
 * @details
 * The function is used when many final seeds have been generated and they need to ordered to check the
 * most promsing hits first. The ordering is calculated from the estimated seed scores.
 *
 * @param[in] s First seed in the comparison.
 * @param[in] o Second seed in the comparison.
 *
 * @returns true iff \a s is greater/more promising than \a o. False is returned if they are equal.
 *
 * @see SRSeed::get_est_score() for documentation on the estimated score.
 * @see test_get_seed_estimated_score() tests functionality.
 */
bool seed_est_score_order_gt(SRSeed const & s, SRSeed const & o);

/*!
 * @brief Help function to get the sort order of the best seeds, greatest score first.
 *
 * @param[in] seeds to get sorted order of
 *
 * @returns Vector containing the sorted order. The returned vector has the same size as the input vector.
 *
 * @see test_get_seed_estimated_score() tests functionality.
 */
std::vector<int> get_seed_est_score_sorted_order_indices(std::vector<SRSeed> const & seeds);

/*!
 * @brief Get seeds from \a read sketches1 and \a read_sketches2
 *
 * @param[out] seeds1 Seeds made from \a read_sketches1 .
 * @param[out] seeds2 Seeds made from \a read_sketches2 .
 * @param[in] read_sketches1 Sketches made from read 1 of the pair.
 * @param[in] read_sketches2 Sketches made from read 2 of the pair.
 * @param[in] gfa GFA graph.
 * @param[in] icu The "I see you" index.
 * @param[in] is_debug Set as true to get extra debugging messages.
 */
void get_seeds_with_pair(std::vector<SRSeed> & seeds1,
                         std::vector<SRSeed> & seeds2,
                         std::vector<ReadSketch> const & read_sketches1,
                         std::vector<ReadSketch> const & read_sketches2,
                         GFA const & gfa,
                         T_icu const & icu,
                         bool const is_debug = false);

} // namespace weaver
