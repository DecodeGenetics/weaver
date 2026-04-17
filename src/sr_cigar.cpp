/*!
 * @file sr_cigar.cpp
 * @brief Implements methods for handling the cigar striing in sam records.
 */

#include <paw/align/alignment_results.hpp> // paw::AlignmentResults
#include <paw/align/cigar.hpp>             // paw::Cigar

#include "logging.hpp"        // print_info
#include "sam_record.hpp"     // SAMRecord
#include "sequence_utils.hpp" // isACGT
#include "sr_alignment.hpp"   //

namespace
{
bool constexpr POS_IS_AFFECTED{true};
bool constexpr POS_IS_NOT_AFFECTED{false};

} // namespace

namespace weaver
{
template <typename Tit>
void left_align_deletion(SAMRecord & main_record, int c, Tit & ref_it, Tit & read_it)
{
  // we must have 'M' operation before and after
  if (c == 0 ||                                                         // not first operation
      c == (static_cast<int>(main_record.cig.size()) - 1) ||            // not last operation
      main_record.cig[c - 1].operation != paw::CigarOperation::MATCH || // operation before must be 'M'
      main_record.cig[c + 1].operation != paw::CigarOperation::MATCH)   // operation after must be 'M'
  {
    return;
  }

  paw::Cigar & prev_cigar = main_record.cig[c - 1];
  int const max_shift = static_cast<int>(prev_cigar.count) - 1;

  assert(max_shift >= 0);

  if (max_shift == 0)
    return;

  paw::Cigar const & curr_cigar = main_record.cig[c];
  assert(curr_cigar.operation == paw::CigarOperation::DELETION);
  paw::Cigar & next_cigar = main_record.cig[c + 1];

  // go one step back
  Tit ref_before_it = std::next(ref_it, -1);
  Tit read_before_it = std::next(read_it, -1);
  Tit ref_after_it = std::next(ref_before_it, curr_cigar.count);

  if (*ref_before_it != *ref_after_it && *ref_before_it == *read_before_it)
  {
    // We can't move to the left
    print_debug(_HERE_,
                " ref_before, read_before, ref_after, read_after=",
                *ref_before_it,
                ", ",
                *read_before_it,
                ", ",
                *ref_after_it,
                ", ",
                *(std::next(read_before_it, 1)));

    print_debug(_HERE_, " del left align NOT possible for qname=", main_record.qname);
    return;
  }

  // adjust the counts of 'M' before and after the insertion
  --prev_cigar.count;
  ++next_cigar.count;

  // update ref_it and read_it to match the left alignment
  std::advance(ref_it, -1);
  std::advance(read_it, -1);

  for (int step{1}; step <= max_shift; ++step)
  {
    // go one step back
    std::advance(ref_before_it, -1);
    std::advance(read_before_it, -1);
    std::advance(ref_after_it, -1);

    if (*ref_before_it != *ref_after_it && *ref_before_it == *read_before_it)
    {
      // print_info(_HERE_, " del left align NOT possible for qname=", main_record.qname);
      break;
    }
    else
    {
      if (step == max_shift)
      {
        print_debug(_HERE_, " left aligned would have resulted in zero count cigar for qname=", main_record.qname);
        print_debug(_HERE_, " prev cigar:", std::to_string(prev_cigar.count), cigar2char(prev_cigar.operation));
        print_debug(_HERE_, " next cigar:", std::to_string(next_cigar.count), cigar2char(next_cigar.operation));
        break;
      }
      // print_info(_HERE_, " del left align possible for qname=", main_record.qname);

      // adjust the counts of 'M' before and after the insertion
      --prev_cigar.count;
      ++next_cigar.count;

      // update ref_it and read_it to match the left alignment
      std::advance(ref_it, -1);
      std::advance(read_it, -1);
    }
  }
}

template <typename Tit>
void left_align_insertion(SAMRecord & main_record, int c, Tit & ref_it, Tit & read_it)
{
  // we must have 'M' operation before and after
  if (c == 0 ||                                                         // not first operation
      c == (static_cast<int>(main_record.cig.size()) - 1l) ||           // not last operation
      main_record.cig[c - 1].operation != paw::CigarOperation::MATCH || // operation before must be 'M'
      main_record.cig[c + 1].operation != paw::CigarOperation::MATCH)   // operation after must be 'M'
  {
    return;
  }

  paw::Cigar & prev_cigar = main_record.cig[c - 1];

  if (prev_cigar.count == 1)
    return;

  paw::Cigar const & curr_cigar = main_record.cig[c];
  assert(curr_cigar.operation == paw::CigarOperation::INSERTION);
  paw::Cigar & next_cigar = main_record.cig[c + 1];

  // go one step back
  Tit ref_before_it = std::next(ref_it, -1);
  Tit read_before_it = std::next(read_it, -1);
  Tit read_after_it = std::next(read_before_it, curr_cigar.count);

  if (*ref_before_it != *read_after_it && *ref_before_it == *read_before_it)
    return;

  // adjust the counts of 'M' before and after the insertion
  --prev_cigar.count;
  ++next_cigar.count;

  // update ref_it and read_it to match the left alignment
  std::advance(ref_it, -1);
  std::advance(read_it, -1);

  while (prev_cigar.count > 1)
  {
    // go one step back
    std::advance(ref_before_it, -1);
    std::advance(read_before_it, -1);
    std::advance(read_after_it, -1);

    if (*ref_before_it != *read_after_it && *ref_before_it == *read_before_it)
      break;

    // adjust the counts of 'M' before and after the insertion
    --prev_cigar.count;
    ++next_cigar.count;

    // update ref_it and read_it to match the left alignment
    std::advance(ref_it, -1);
    std::advance(read_it, -1);
  }
}

bool process_cigar_trivial(weaver::SAMRecord & main,
                           paw::AlignmentResults const & ar,
                           int const num_read_bases,
                           bool is_pos_affected)
{
  if (ar.cigar_string_ptr == nullptr)
  {
    // cigar_string_ptr==nullptr iff alignment is only match or mismatch ('M' operation)
    if (is_pos_affected)
      main.pos -= num_read_bases;

    main.append_cigar(num_read_bases, paw::CigarOperation::MATCH);
    return true;
  }

  std::vector<paw::Cigar> const & cigar_string = *(ar.cigar_string_ptr);

  if (cigar_string.size() == 0 || ar.clip_begin > 0 || ar.clip_end == 0)
  {
    // Soft clip the remaining read bases
    main.append_cigar(num_read_bases, paw::CigarOperation::SOFT_CLIP);
    return true; // Indicates that the cigar has been processed.
  }

  return false;
}

void process_cigar_forward_begin(weaver::SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_bases)
{
  if (process_cigar_trivial(main, ar, num_read_bases, POS_IS_AFFECTED))
  {
    print_debug(_HERE_, " Processed forward begin cigar trivially.");
    assert(main.get_cigar_reference_length() == num_read_bases || main.get_cigar_reference_length() == 0);
    assert(main.get_cigar_query_length() == num_read_bases);
    return;
  }

  assert(ar.cigar_string_ptr != nullptr);
  std::vector<paw::Cigar> const & cigar_string = *(ar.cigar_string_ptr);

#ifndef NDEBUG
  print_debug(_HERE_, " alignment: ", alignment_results_to_string(ar));
  int query_count{0};
#endif // NDEBUG

  for (auto cig_it = cigar_string.begin(); cig_it != cigar_string.end(); ++cig_it)
  {
    int cigar_count = cig_it->count;
    paw::CigarOperation inv_op = inv_cigar(cig_it->operation);
    main.append_cigar(cigar_count, inv_op);

#ifndef NDEBUG
    if (paw::advances_query(inv_op))
      query_count += cigar_count;
#endif // NDEBUG

    if (paw::advances_ref(inv_op))
      main.pos -= cigar_count;
  }

#ifndef NDEBUG
  print_debug(_HERE_,
              " (fwd begin) qname=",
              main.qname,
              " cigar=",
              main.get_cigar(),
              " num_read_bases=",
              num_read_bases,
              " query_count=",
              query_count);

  assert(query_count == num_read_bases);
#endif // NDEBUG
}

void process_cigar_forward_end(weaver::SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_bases)
{
  if (process_cigar_trivial(main, ar, num_read_bases, POS_IS_NOT_AFFECTED))
  {
    print_debug(_HERE_, " (fw end) qname=", main.qname, " processed forward end cigar trivially.");
    return;
  }

  assert(ar.cigar_string_ptr != nullptr);
  std::vector<paw::Cigar> & cigar_string = *(ar.cigar_string_ptr);

#ifndef NDEBUG
  print_debug(_HERE_, " alignment: ", alignment_results_to_string(ar));
  int query_count{0};
#endif // NDEBUG

  for (auto cig_it = cigar_string.rbegin(); cig_it != cigar_string.rend(); ++cig_it)
  {
    int cigar_count = cig_it->count;
    paw::CigarOperation inv_op = inv_cigar(cig_it->operation);

#ifndef NDEBUG
    if (paw::advances_query(inv_op))
      query_count += cigar_count;
#endif // NDEBUG

    main.append_cigar(cigar_count, inv_op);
  }

#ifndef NDEBUG
  print_debug(_HERE_,
              " (fwd end) qname=",
              main.qname,
              " cigar=",
              main.get_cigar(),
              " num_read_bases=",
              num_read_bases,
              " query_count=",
              query_count);

  assert(query_count == num_read_bases);
#endif // NDEBUG
}

void process_cigar_reverse_end(weaver::SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_bases)
{
  if (process_cigar_trivial(main, ar, num_read_bases, POS_IS_AFFECTED))
  {
    print_debug(_HERE_, " Processed reverse end cigar trivially.");
    return;
  }

  assert(ar.cigar_string_ptr != nullptr);
  std::vector<paw::Cigar> const & cigar_string = *(ar.cigar_string_ptr);

#ifndef NDEBUG
  print_debug(_HERE_, " alignment: ", alignment_results_to_string(ar));
  int query_count{0};
#endif // NDEBUG

  for (auto cig_it = cigar_string.begin(); cig_it != cigar_string.end(); ++cig_it)
  {
    int cigar_count = cig_it->count;
    paw::CigarOperation inv_op = inv_cigar(cig_it->operation);
    main.append_cigar(cigar_count, inv_op);

#ifndef NDEBUG
    if (paw::advances_query(inv_op))
      query_count += cigar_count;
#endif // NDEBUG

    if (paw::advances_ref(inv_op))
      main.pos -= cigar_count;
  }

#ifndef NDEBUG
  print_debug(_HERE_,
              " (rev end) qname=",
              main.qname,
              " cigar=",
              main.get_cigar(),
              " num_read_bases=",
              num_read_bases,
              " query_count=",
              query_count);

  assert(query_count == num_read_bases);
#endif // NDEBUG
}

void process_cigar_reverse_begin(weaver::SAMRecord & main, paw::AlignmentResults const & ar, int const num_read_bases)
{
  if (process_cigar_trivial(main, ar, num_read_bases, POS_IS_NOT_AFFECTED))
  {
    print_debug(_HERE_, " Processed reverse begin cigar trivially.");
    return;
  }

  assert(ar.cigar_string_ptr != nullptr);
  std::vector<paw::Cigar> const & cigar_string = *(ar.cigar_string_ptr);

#ifndef NDEBUG
  print_debug(_HERE_, " alignment: ", alignment_results_to_string(ar));
  int query_count{0};
#endif // NDEBUG

  for (auto cig_it = cigar_string.rbegin(); cig_it != cigar_string.rend(); ++cig_it)
  {
    int cigar_count = cig_it->count;
    paw::CigarOperation inv_op = inv_cigar(cig_it->operation);

#ifndef NDEBUG
    if (paw::advances_query(inv_op))
      query_count += cigar_count;
#endif // NDEBUG

    main.append_cigar(cigar_count, inv_op);
  }

#ifndef NDEBUG
  print_debug(_HERE_,
              " (rev begin) qname=",
              main.qname,
              " cigar=",
              main.get_cigar(),
              " num_read_bases=",
              num_read_bases,
              " query_count=",
              query_count);

  assert(query_count == num_read_bases);
#endif // NDEBUG
}

int get_cigar_query_length(std::vector<paw::Cigar> const & cigar_string)
{
  int query_len{0};

  for (paw::Cigar const & cigar : cigar_string)
  {
    if (paw::advances_query(cigar.operation))
      query_len += cigar.count;
  }

  return query_len;
}

int get_cigar_reference_length(std::vector<paw::Cigar> const & cigar_string)
{
  int ref_len{0};

  for (paw::Cigar const & cigar : cigar_string)
  {
    if (paw::advances_ref(cigar.operation))
      ref_len += cigar.count;
  }

  return ref_len;
}

// Explicit instantiations
//
//! Explicit instantiation of left_align_deletion() for strings.
template void left_align_deletion(SAMRecord & main_record,
                                  int c,
                                  std::string::const_iterator & ref_it,
                                  std::string::const_iterator & read_it);

//! Explicit instantiation of left_align_insertion() for strings.
template void left_align_insertion(SAMRecord & main_record,
                                   int c,
                                   std::string::const_iterator & ref_it,
                                   std::string::const_iterator & read_it);

} // namespace weaver
