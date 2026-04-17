#pragma once

/*!
 * @file sr_cigar.hpp
 * @brief Defines methods for handling the cigar string in sam records.
 */

#include <paw/align/alignment_results.hpp>

#include "sam_record.hpp"

namespace weaver
{
/*!
 * @brief Left align a deletion cigar operation at \a main_record.cig[c].

 * @param[in,out] main_record The SAM record containg the cigar string to modify.
 * @param[in] c index of the deletion cigar operation.
 * @param[in,out] ref_it Iterator to the reference sequence.
 * @param[in,out] read_it Iterator to the read sequence.
 *
 * @see left_align_insertion()
 */
template <typename Tit>
void left_align_deletion(SAMRecord & main_record, int c, Tit & ref_it, Tit & read_it);

/*!
 * @brief Left align a insertion cigar operation at \a main_record.cig[c].

 * @param[in,out] main_record The SAM record containg the cigar string to modify.
 * @param[in] c index of the deletion cigar operation.
 * @param[in,out] ref_it Iterator to the reference sequence.
 * @param[in,out] read_it Iterator to the read sequence.
 *
 * @see left_align_deletion()
 */
template <typename Tit>
void left_align_insertion(SAMRecord & main_record, int c, Tit & ref_it, Tit & read_it);

/*!
 * @brief Process the cigar string in a trivial manner
 *
 * @details
 * Trivial means either all are matches or everything is soft clipped.
 *
 * @returns True iff the cigar could be processed trivially.
 */
bool process_cigar_trivial(SAMRecord & main,
                           paw::AlignmentResults const & ar,
                           int const num_read_bases,
                           bool is_pos_affected);

//! Process cigar for the beginning of a forward oriented read.
void process_cigar_forward_begin(SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_begin_bases);

//! Process cigar for the end of a forward oriented read.
void process_cigar_forward_end(SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_end_bases);

//! Process cigar for the beginning of a reverse oriented read.
void process_cigar_reverse_begin(SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_end_bases);

//! Process cigar for the end of a reverse oriented read.
void process_cigar_reverse_end(SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_end_bases);

//! Returns the query length of a CIGAR string.
int get_cigar_query_length(std::vector<paw::Cigar> const & cigar_string);

//! Returns the query length of a CIGAR string.
int get_cigar_reference_length(std::vector<paw::Cigar> const & cigar_string);

inline std::string cigar_element2string(paw::Cigar const & c)
{
  std::string str;
  str += std::to_string(c.count);
  str += cigar2char(c.operation);
  return str;
}

//! Append a cigar operation \a op with count \a count to the back of cigar string \a cig
inline void append_cigar_back(std::vector<paw::Cigar> & cig, uint32_t count, paw::CigarOperation op)
{
  auto const n_cigar{cig.size()};

  if (n_cigar == 0 || cig[n_cigar - 1].operation != op)
    cig.emplace_back(count, op);
  else
    cig[n_cigar - 1].count += count;
}

//! Append a cigar operation \a op with count \a count to the front of cigar string \a cig
inline void append_cigar_front(std::vector<paw::Cigar> & cig, uint32_t count, paw::CigarOperation op)
{
  if (cig.size() == 0)
    cig.emplace_back(count, op);
  else if (cig[0].operation != op)
    cig.emplace(cig.begin(), count, op);
  else
    cig[0].count += count;
}

template <typename Tit>
inline std::string inv_cigar2string(Tit begin, Tit end)
{
  std::string str;

  while (begin != end)
  {
    paw::Cigar const & c = *begin;
    str += std::to_string(c.count);
    str += inv_cigar2char(c.operation);
    ++begin;
  }

  return str;
}

template <typename Tit>
inline std::string cigar2string(Tit begin, Tit end)
{
  std::string str;

  while (begin != end)
  {
    paw::Cigar const & c = *begin;
    str += std::to_string(c.count);
    str += cigar2char(c.operation);
    ++begin;
  }

  return str;
}

} // namespace weaver
