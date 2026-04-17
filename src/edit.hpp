#pragma once

#include <cassert>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "logging.hpp"
#include "sequence_utils.hpp" // isACGTN

namespace weaver
{
class Edit
{
public:
  /*!
   * @name Constants
   * @{
   */

  //! Value to use when edit call is missing.
  uint8_t static constexpr MISSING_CALL{std::numeric_limits<uint8_t>::max()};

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  int pos{-1};             //!< Edit chromasomal position, 0-based indexing. -1 if missing.
  char type{'\0'};         //!< Types are: 0=NA, X=SNP, I=insertion, D=deletion.
  std::vector<char> seq{}; //!< Edit sequence. Expected to be A,C,G,T or N.

  /*!
   * @}
   *
   * @name Constructors
   * @{
   */

  //! Empty constructor is defaulted
  Edit() = default;

  /*!
   * @brief Construct a new edit instance.
   *
   * @param[in] _pos  0-based position of the edit.
   * @param[in] _type type of the edit: X=SNP, I=insertion, D=deletion.
   * @param[in] _seq  Sequence of the edit. Expected to be A,C,G,T or N.
   */

  Edit(int _pos, char _type, std::vector<char> const & _seq) : pos(_pos), type(_type), seq(_seq)
  {
#ifndef NDEBUG
    for (char c : seq)
      assert(isACGTN(c));
#endif // NDEBUG
  }

  //! @brief Construct a new edit instance of a SNP.
  Edit(int _pos, char _type, char _snp) : pos(_pos), type(_type), seq(1, _snp)
  {
    assert(_type == 'X');
    assert(isACGTN(_snp));
  }

  /*!
   * @}
   *
   * @name Read-only methods
   */

  //! Checks if the edit type is a SNP
  inline bool is_snp() const
  {
    return type == 'X';
  }

  //! Checks if the edit type is an insertion
  inline bool is_insertion() const
  {
    return type == 'I';
  }

  //! Checks if the edit type is a deletion
  inline bool is_deletion() const
  {
    return type == 'D';
  }

  //! Checks if the edit type is a missing call
  inline bool is_missing() const
  {
    return type == '\0';
  }

  inline bool operator==(Edit const & b) const
  {
    return pos == b.pos && type == b.type && seq == b.seq;
  }

  inline bool operator!=(Edit const & b) const
  {
    return !(*this == b);
  }

  inline bool operator<(Edit const & b) const
  {
    // sorted order of edit is insertion, deletion, snp
    auto const order_a = (type == 'D') + 2 * (type == 'X');
    auto const order_b = (b.type == 'D') + 2 * (b.type == 'X');

    return pos < b.pos || (pos == b.pos && order_a < order_b) || (pos == b.pos && order_a == order_b && seq < b.seq);
  }

  inline int get_order() const
  {
    return (type == 'D') + 2 * (type == 'X');
  }

  std::string to_string() const;

  /*!
   * @}
   */
};

class EditHash
{
public:
  std::size_t operator()(Edit const & e) const;
};

class EditCalls
{
public:
  /*!
   * @name Constants
   * @{
   */
  uint8_t static constexpr MISSING_CALL{std::numeric_limits<uint8_t>::max()};

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  int pos{-1};                //!< Edit chromasomal position, 0-based indexing. -1 if missing.
  char type{'\0'};            //!< Types are: 0=NA, X=SNP, I=insertion, D=deletion.
  std::vector<char> seq{};    //!< Edit sequence.
  std::vector<uint8_t> calls; //!< Edit haplotype calls

  /*!
   * @}
   *
   * @name Constructors
   * @{
   */

  EditCalls() = default;

  EditCalls(int _pos, char _type, std::vector<char> const & _seq, std::vector<uint8_t> && _calls) :
    pos(_pos), type(_type), seq(_seq), calls(std::forward<std::vector<uint8_t>>(_calls))
  {
#ifndef NDEBUG
    for (char c : seq)
      assert(isACGTN(c));
#endif // NDEBUG
  }

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  inline bool operator==(EditCalls const & b) const
  {
    return pos == b.pos && type == b.type && seq == b.seq;
  }

  inline bool operator==(Edit const & b) const
  {
    return pos == b.pos && type == b.type && seq == b.seq;
  }

  inline bool operator!=(EditCalls const & b) const
  {
    return !(*this == b);
  }

  inline bool operator!=(Edit const & b) const
  {
    return !(*this == b);
  }

  inline bool operator<(EditCalls const & b) const
  {
    // sorted order of edit is insertion, deletion, snp
    auto const order_a = (type == 'D') + 2 * (type == 'X');
    auto const order_b = (b.type == 'D') + 2 * (b.type == 'X');

    return pos < b.pos || (pos == b.pos && order_a < order_b) || (pos == b.pos && order_a == order_b && seq < b.seq);
  }

  std::string to_string() const;

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  //! This method lets cereal know which data members to serialize.
  template <class Archive>
  void serialize(Archive & archive)
  {
    for (char c : seq)
    {
      if (!isACGTN(c))
        print_error(_HERE_, " bad base ", c);
    }

    archive(pos, type, seq, calls);

    for (char c : seq)
    {
      if (!isACGTN(c))
        print_error(_HERE_, " bad base ", c);
    }
  }

  /*!
   * @}
   */
};

//! Add a pair of haplotype matches.
void add_hap_matches(std::vector<double> & hap_matches,
                     std::vector<double> const & snp_hap_matches,
                     std::vector<double> const & nonsnp_hap_matches);

/*!
 * @brief Get the number of haplotype matches to a collection of edit calls.
 *
 * @details
 * In case of missing allele data in haplotype, they are also counted as matches.
 *
 * @returns The number of reference matches.
 */
int get_hap_matches(std::vector<double> & hap_matches,
                    std::vector<EditCalls> const & ec,
                    std::vector<int> const & edits);

/*!
 * @brief Get the number of matches
 */
int get_ref_matches(std::vector<int> const & edits);

/*!
 * @brief Method for sorting edit indexes with std::sort
 *
 * @details
 * Special care needs to be taking to handling also the no edits/negative values.
 */
inline bool is_less_edit_index(int a, int b)
{
  return (a < b && (-a - 1) < b) || (a < (-b - 1) && -a - 1 < (-b - 1));
}

/*!
 * @brief Alternative method for sorting edit indexes with std::sort
 *
 * @see is_less_edit_index
 */
inline bool is_less_edit_index_alt(int a, int b)
{
  return (a < 0 ? -a - 1 : a) < (b < 0 ? -b - 1 : b);
}

} // namespace weaver
