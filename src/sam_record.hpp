#pragma once

#include <cstdint> // uint16_t
#include <limits>  // std::numeric_limits
#include <string>  // std::string
#include <vector>  // std::vector

#include <paw/align/cigar.hpp>

#include <weaver/constants.hpp>

namespace weaver
{
class FastqData;

//! Class that has all fields for writing a SAM record
class SAMRecord
{
public:
  /*!
   * @name Constants
   * @{
   */

  //! Missing tag values will be set as this constant.
  int static constexpr MISSING_TAG{std::numeric_limits<int>::lowest() / 2};

  //! Missing position value.
  int static constexpr MISSING_POS{-1};

  /*!
   * @}
   *
   * @name Required SAM fields
   * @{
   */

  std::string_view qname{};      //!< query template name.
  std::string_view extra_tags{}; //!< Extra tags, for example a read group (RG) tag.
  std::string_view seq{};        //!< query sequence.
  std::string_view qual{};       //!< ASCII of Phred-scaled base qualities.

  int sfa_idx{-1};      // stable contig index.
  int snid{-1};         // stable contig name index.
  int pos{MISSING_POS}; //!< 1-based leftmost mapping position.
  int mapq{0};          //!< Mapping quality.
  int tlen{0};          //!< observed template length. 0 if missing

  std::vector<paw::Cigar> cig{}; //!< CIGAR string
  uint16_t flags{0};             //!< bitwise flags.

  /*!
   * @}
   *
   * @name Optional SAM tags
   * @{
   */

  int num_edits{MISSING_TAG};              //!< NM tag. The number of edits between read and reference.
  int alignment_score{MISSING_TAG};        //!< AS tag. The alignment score.
  int weaver_score{MISSING_TAG};           //!< WS tag. Weaver's alignments score, used for calculating mapping quality.
  int second_alignment_score{MISSING_TAG}; //!< XS tag. The secondary alignment score, may be larger than AS.

  /*!
   * @}
   *
   * @name Constructors and destructor
   * @{
   */

  SAMRecord() = default;
  SAMRecord(SAMRecord const & other_sam_record) = default;
  explicit SAMRecord(FastqData const & fastq_data);
  ~SAMRecord() = default;

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  /*!
   * @brief Appends the cigar string.
   *
   * @details
   * If the same operation is seen twice in a row, they are merged instead of adding new entries.
   *
   * @see append_cigar(paw::Cigar const & cigar)
   */
  void append_cigar(int count, paw::CigarOperation op);

  /*!
   * @brief Appends the cigar string.
   *
   * @details
   * If the same operation is seen twice in a row, they are merged instead of adding new entries.
   *
   * @see append_cigar(uint32_t count, paw::CigarOperation op)
   */
  void append_cigar(paw::Cigar const & cigar);

  /*!
   * @brief Replace the beginning of a CIGAR based on the amount of reference bases with soft clip.
   *
   * @details
   * Replacing bases with a soft clip will change the position of the record.
   *
   * @returns the amount of query bases clipped.
   */
  int soft_clip_cigar_begin_based_on_ref(int num_ref_bases);

  /*!
   * @brief Replace the end of a CIGAR based on the amount of reference bases with soft clip.
   *
   * @returns the amount of query bases clipped.
   */
  int soft_clip_cigar_end_based_on_ref(int num_ref_bases);

  /*!
   * @brief Sets flags and clears fields that indicate an unmapped read alignment.
   *
   * @details
   * Useful if after aligning there are some clear problems detected. A reference to the other sam record is needed
   * such that flags are correctly set in both records.
   */
  void make_unmapped(SAMRecord & other_sam_record);

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  //! Generate a std::string represention of the cigar string.
  std::string get_cigar() const;

  //! Get how many query bases are in the cigar string.
  int get_cigar_query_length() const;

  //! Get how many reference bases are in the cigar string.
  int get_cigar_reference_length() const;

  /*!
   * @brief Return the SAM record order value.
   *
   * @details
   * The SAM record order value specifies how the record should be sorted. I.e. records with lower order should appear
   * first in the final SAM output. If two records have the same order, then it is undefined which should appear first.
   *
   * The order value is made of two 32-bit values: sfa_idx (32 bits) | position (32 bits)
   */
  uint64_t get_sam_order() const;

  //! Return how far the query reaches on the reference.
  int get_reference_reach() const;

  //! Returns true iff read was clipped on both ends with more \a min_clip_length bases
  bool is_double_clipped(int const min_clip_length = 10) const;

  //! Checks the flags if this read is unmapped and its mate is mapped
  bool is_unmapped_with_mapped_mate() const;

  /*!
   * @}
   *
   * @name Debugging methods
   * @{
   */

  /*!
   * @brief Checks the record cigar string for any unregularities.
   *
   * @details
   * Will return true if the CIGAR string is valid or empty.
   */
  bool is_cigar_valid() const;

  /*!
   * @brief Check if the SAM record pair looks valid
   *
   * @returns false iff problems are found in the SAM record.
   * @see is_valid() is also called for running other validity checks. Call is_valid() instead, if the record is missing
   *                 paired-end information.
   */
  bool is_pair_valid(SAMRecord const & other) const;

  /*!
   * @brief Check if the SAM record looks valid.
   *
   * @returns false iff problems are found in the SAM record.
   * @see is_pair_valid() should be called instead if the record has paired-end information.
   */
  bool is_valid() const;

  /*!
   * @}
   */
};

} // namespace weaver
