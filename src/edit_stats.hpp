#pragma once

#include <cstdint>
#include <limits>
#include <vector>

namespace weaver
{
class Edit;
class EditCalls;

//! Observed statistics such as read counts and haplotype likelihoods for variants stored in EditCalls objects.
class EditStats
{
public:
  /*!
   * @name Constructors and destructors
   * @{
   */

  /*!
   * Construct a new EditStats object given a variant.
   *
   * @brief
   * Stats for the reference allele are added at construction.
   */
  explicit EditStats(EditCalls const & edit_calls);

  EditStats() = default;                  //!< Default empty constructor.
  EditStats(EditStats const &) = default; //!< Default copy constructor.
  EditStats(EditStats &&) = default;      //!< Default move constructor.
  ~EditStats() = default;                 //! Default deconstructor.

  /*!
   * @}
   *
   * @name Operators
   * @{
   */

  EditStats & operator=(EditStats const &) = default; //!< Default copy assignment
  EditStats & operator=(EditStats &&) = default;      //!< Default move assignment

  //! Check if two edit stats have the same stats.
  bool operator==(EditStats const &) const;

  //! Check if two edit stats have different stats.
  bool operator!=(EditStats const & o) const;

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  //! Number of haplotypes with a non-missing all. Includes the reference haplotype which is set at initialization.
  uint16_t non_missing_haplotypes{0};

  //! Records how often the edit is seen in read data.
  uint16_t edit_read_count{0};

  //! Records how often a no edit is seen in read data.
  uint16_t no_edit_read_count{0};

  //! Number of haplotypes with each edit.
  uint16_t count{0};

  double edit_prob{0.0};

  //! If genotype is homozygous, get the homozygous call. Otherwise return Edit::MISSING_CALL
  uint8_t get_hom_call() const;

  double get_hom_no_edit_eps() const;
  double get_hom_edit_eps() const;

  // double alpha_edit_prob_log{std::numeric_limits<double>::lowest()};
  // double beta_edit_prob_log{std::numeric_limits<double>::lowest()};
  double alpha_beta_log_sum{std::numeric_limits<double>::lowest()};

  std::vector<double> alpha_beta;

  //! Posterior probabilities of each likelihoods of based on HMM
  // std::vector<double> post_prob_log{};

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  //! Get the PHRED probability of homozygous edit
  double get_hom_edit_phred() const;

  //! Get the PHRED probability of heterozygousity
  double get_het_phred() const;

  //! Get the PHRED probability of homozygous no edit
  double get_hom_no_edit_phred() const;

  //! Get the frequency of an edit.
  double get_edit_frequency() const;

  //! Get the frequency of no edit.
  double get_no_edit_frequency() const;

  //! Get the likelihood that an edit is present
  double get_likelihood_of_edit() const;

  //! Get the likelihood that an edit is NOT present
  double get_likelihood_of_no_edit() const;

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  //! Merge other edit stats into this object.
  void merge_with(EditStats const & other_edit_stats);

  /*!
   * @}
   */
};

} // namespace weaver
