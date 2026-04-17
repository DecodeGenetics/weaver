/*!
 * @file sequence_utils.hpp
 * @brief Implements various functions for operating on sequences.
 */

#include "sequence_utils.hpp"

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/pairwise_alignment.hpp>

#include "logging.hpp"
#include "options.hpp"
#include "sr_cigar.hpp"

namespace weaver
{
template <typename Tstring>
std::string get_reverse_complement(Tstring const & fwd)
{
  std::string rev(fwd.size(), '\0');
  std::transform(fwd.rbegin(), fwd.rend(), rev.begin(),
                 complement); // reverse complement
  return rev;
}

template <typename Tstring1, typename Tstring2>
int count_prefix_mismatches(Tstring1 const & s1, Tstring2 const & s2)
{
  int const min_size = std::min(s1.size(), s2.size());
  int count{0};

  for (int i{0}; i < min_size; ++i)
    count += (s1[i] != s2[i] && isACGT(s1[i]) && isACGT(s2[i]));

  return count;
}

template <typename Tstring1, typename Tstring2>
int count_prefix_mismatches_s2rev(Tstring1 const & s1, Tstring2 const & s2)
{
  int const min_size = std::min(s1.size(), s2.size());
  int count{0};

  for (int i{0}; i < min_size; ++i)
  {
    int j = s2.size() - 1 - i;
    count += (s1[i] != complement(s2[j]) && isACGT(s1[i]) && isACGT(s2[j]));
  }

  return count;
}

template <typename Tstring>
std::vector<std::string_view> split_string(Tstring const & str, char const delimiter, int const max_occurance)
{
  std::vector<std::string_view> output;
  std::string_view strv(str);
  auto first = strv.cbegin();

  while (first != strv.cend())
  {
    auto const second = std::find(first, strv.cend(), delimiter);

    if (first != second)
    {
      std::size_t const pos = std::distance(strv.cbegin(), first);
      output.emplace_back(strv.substr(pos, second - first));
    }

    if (second == strv.cend() || static_cast<int>(output.size()) == max_occurance)
      break;

    first = std::next(second);
  }

  return output;
}

void find_mem_with_max_length(int & mem_read_start, //
                              int & mem_ref_start,
                              int & mem_length,
                              int read_start,
                              int ref_start,
                              int const max_length,
                              std::string const & read,
                              std::string const & ref)
{
  if (max_length <= mem_length)
    return; // small optimiziation becaues in this case it is impossible that this segment can be the mem

  int cur_length{0};
  int const read_end = read_start + max_length;

  while (read_start < read_end)
  {
    assert(read_start < static_cast<int>(read.size()));
    assert(ref_start < static_cast<int>(ref.size()));

    if (read[read_start] == ref[ref_start])
    {
      ++cur_length; // match
    }
    else
    {
      // mismatch
      if (cur_length > mem_length)
      {
        mem_read_start = read_start - cur_length;
        mem_ref_start = ref_start - cur_length;
        mem_length = cur_length;
      }

      cur_length = 0;
    }

    ++read_start;
    ++ref_start;
  }

  if (cur_length > mem_length)
  {
    mem_read_start = read_start - cur_length;
    mem_ref_start = ref_start - cur_length;
    mem_length = cur_length;
  }
}

void find_mem(int & mem_read_start, //
              int & mem_ref_start,
              int & mem_length,
              std::string const & read,
              std::string const & ref)
{
  assert(mem_read_start <= 0);
  assert(mem_length <= 0);

  if (read.size() == ref.size())
  {
    find_mem_with_max_length(mem_read_start,
                             mem_ref_start,
                             mem_length,
                             /*read_start=*/0,
                             /*ref_start=*/0,
                             /*max_length=*/read.size(),
                             read,
                             ref);

    // print_info(_HERE_, " mem read start=", mem_read_start, " ref start=", mem_ref_start, " length=", mem_length);
  }
  else
  {
    //
    paw::AlignmentOptions<uint8_t> aln_opts;
    Options const & copts = *(Options::const_instance());

    aln_opts.set_match(copts.match)
      .set_mismatch(copts.mismatch)
      .set_gap_open(copts.gap_open)
      .set_gap_extend(copts.gap_extend);

    // aln_opts.get_aligned_strings = true;
    aln_opts.get_cigar_string = true;

    paw::pairwise_alignment(ref, read, aln_opts);

    assert(aln_opts.ar);
    assert(aln_opts.ar->cigar_string_ptr);

    std::vector<paw::Cigar> const & cigar_string = *(aln_opts.ar->cigar_string_ptr);
    int read_start{0};
    int ref_start{0};

    for (auto it = cigar_string.rbegin(); it != cigar_string.rend(); ++it)
    {
      paw::CigarOperation const & op = it->operation;
      int const max_length = static_cast<int>(it->count);
      // print_debug(_HERE_, " c=", max_length, paw::cigar2char(op));

      if (op == paw::CigarOperation::MATCH)
      {
        find_mem_with_max_length(mem_read_start,
                                 mem_ref_start,
                                 mem_length,
                                 read_start,
                                 ref_start,
                                 max_length,
                                 read,
                                 ref);

        read_start += max_length;
        ref_start += max_length;
      }
      else
      {
        if (paw::advances_query(op))
          read_start += max_length;

        if (paw::advances_ref(op))
          ref_start += max_length;
      }
    }

    /*
    print_info(_HERE_, " mem read start=", mem_read_start, " ref start=", mem_ref_start, " length=", mem_length);

    print_info(_HERE_,
               " read=",
               read,
               " ref=",
               ref,
               " cigar=",
               inv_cigar2string(cigar_string.rbegin(), cigar_string.rend()));
    */
    //
  }
}

int get_prefix_size(std::vector<char> const & s1, std::vector<char> const & s2)
{
  int const min_size = std::min(s1.size(), s2.size());

  for (int i{0}; i < min_size; ++i)
  {
    if (s1[0] != s2[0])
      return i;
  }

  return min_size;
}

int get_suffix_size(std::vector<char> const & s1, std::vector<char> const & s2)
{
  int const n1 = static_cast<int>(s1.size());
  int const n2 = static_cast<int>(s2.size());
  int const min_size = std::min(n1, n2);

  for (int s{1}; s <= min_size; ++s)
  {
    if (s1[n1 - s] != s2[n2 - s])
      return s - 1;
  }

  return min_size;
}

// explicit instantiation
//
template std::string get_reverse_complement<std::string>(std::string const & fwd);
template std::string get_reverse_complement<std::string_view>(std::string_view const & fwd);

template int count_prefix_mismatches(std::string const & s1, std::string const & s2);
template int count_prefix_mismatches(std::string_view const & s1, std::string const & s2);
template int count_prefix_mismatches(std::string const & s1, std::string_view const & s2);
template int count_prefix_mismatches(std::string_view const & s1, std::string_view const & s2);

template int count_prefix_mismatches_s2rev(std::string const & s1, std::string const & s2);
template int count_prefix_mismatches_s2rev(std::string_view const & s1, std::string const & s2);
template int count_prefix_mismatches_s2rev(std::string const & s1, std::string_view const & s2);
template int count_prefix_mismatches_s2rev(std::string_view const & s1, std::string_view const & s2);

template std::vector<std::string_view> split_string(std::string const & str,
                                                    char const delimiter,
                                                    int const max_occurance);

template std::vector<std::string_view> split_string(std::string_view const & str,
                                                    char const delimiter,
                                                    int const max_occurnace);

} // namespace weaver
