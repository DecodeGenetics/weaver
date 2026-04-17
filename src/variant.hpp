#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "io.hpp"

namespace weaver
{
/*!
 * @brief Small variant site data.
 *
 * @headerfile small_variant.hpp "weaver/small_variant.hpp"
 *
 * @details
 * Stores alleles and calls of a small variant.
 */
class Variant
{
public:
  /*!
   * @name Constants
   * @{
   */

  //! The maximum allele size for an allele to be considered small.
  int static constexpr MAX_SIZE{70};
  uint16_t static constexpr MISSING_CALL{std::numeric_limits<uint16_t>::max()};

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  //! Position of the variant (1-based indexing). -1 if missing.
  int pos{-1};

  //! Allele sequence at variant site.
  std::vector<std::vector<char>> seqs{};

  /*!
   * Indicates which allele which haplotype has calls.
   *
   * @brief
   * Its length is the number of haplotypes. Missing call is denoted with \a this->MISSING_CALL
   */
  std::vector<uint16_t> calls;

  /*!
   * @}
   *
   * @name Read-only methods
   */

  //! Determines the ordering of small variants for ordered containers. Call data is not considered for ordering.
  bool operator<(Variant const & o) const
  {
    return pos < o.pos || (pos == o.pos && seqs < o.seqs);
  }

  //! Determines equality of small variants. Call data is not considered for checking equality.
  bool operator==(Variant const & o) const
  {
    return pos == o.pos && seqs == o.seqs;
  }

  //! Determines inequality of small variants. Call data is not considered for checking inequality.
  bool operator!=(Variant const & o) const
  {
    return pos != o.pos || seqs != o.seqs;
  }

  //! Make a string describing the variant. For debugging.
  inline std::string to_string() const
  {
    std::string str;
    str += std::to_string(pos) + ' ';

    if (seqs.size() > 0)
    {
      str += std::string(seqs[0].begin(), seqs[0].end());

      for (int i{1}; i < static_cast<int>(seqs.size()); ++i)
        str += ',' + std::string(seqs[i].begin(), seqs[i].end());

      if (calls.size() > 0)
      {
        str += ' ' + std::to_string(calls[0]);

        for (int i{1}; i < static_cast<int>(calls.size()); ++i)
          str += ',' + std::to_string(calls[i]);
      }
    }

    return str;
  }

  //! This method lets cereal know which data members to serialize.
  template <class Archive>
  void serialize(Archive & archive)
  {
    archive(pos, seqs, calls);
  }

  /*!
   * @}
   */
};

/*!
 * @brief Read VCF data and extract all variant site data.
 *
 * @param[in] in_vcf Smart pointer to a hts file containing VCF data.
 * @param[in] in_tbx Smart pointer to tabix data.
 * @param[in] in_it Smart pointer to tabix query information.
 *
 * @returns All variants of the region.
 */
std::vector<Variant> get_variants_in_a_region(hts_file_ptr const & in_vcf,
                                              tbx_t_ptr const & in_tbx,
                                              hts_itr_t_ptr const & in_it);

} // namespace weaver
