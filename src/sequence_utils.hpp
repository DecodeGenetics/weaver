#pragma once
/*!
 * @file sequence_utils.hpp
 * @brief Defines various function for operating on sequences.
 */

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace weaver
{
//! Complements a base
inline char complement(char c)
{
  switch (c)
  {
  case 'a':
  case 'A':
    return 'T';

  case 'c':
  case 'C':
    return 'G';

  case 'g':
  case 'G':
    return 'C';

  case 't':
  case 'T':
    return 'A';

  default:
    return c;
  }
}

//! Maps a->0,c->1,g->2,t|u->3,others->4
inline uint64_t nt5_char_to_ull(char nt)
{
  switch (nt)
  {
  case 'a':
  case 'A':
    return 0;

  case 'c':
  case 'C':
    return 1;

  case 'g':
  case 'G':
    return 2;

  case 't':
  case 'T':
  case 'u':
  case 'U':
    return 3;

  default:
    return 4;
  }
}

//! Complement version of nt5_char_to_ull(char nt)
inline uint64_t complement_nt5_char_to_ull(char nt)
{
  switch (nt)
  {
  case 't':
  case 'T':
  case 'u':
  case 'U':
    return 0;

  case 'g':
  case 'G':
    return 1;

  case 'c':
  case 'C':
    return 2;

  case 'a':
  case 'A':
    return 3;

  default:
    return 4;
  }
}

//! Returns true iff character \a c is either A or C or G or T (case insensitive).
inline bool isACGT(char c)
{
  switch (c)
  {
  case 'a':
  case 'A':
  case 'c':
  case 'C':
  case 'g':
  case 'G':
  case 't':
  case 'T':
    return true;

  default:
    return false;
  }

  return false;
}

//! Returns true iff character \a c is either A or C or G or T or N (case insensitive).
inline bool isACGTN(char c)
{
  switch (c)
  {
  case 'a':
  case 'A':
  case 'c':
  case 'C':
  case 'g':
  case 'G':
  case 't':
  case 'T':
  case 'n':
  case 'N':
    return true;

  default:
    return false;
  }

  return false;
}

//! Complements a DNA base. i.e. A<>T, C<>G
char complement(char c);

//! return the reverse complement of a string
template <typename Tstring>
std::string get_reverse_complement(Tstring const & fwd);

//! Count the number of bases match in the maximal prefix
template <typename Tstring1, typename Tstring2>
int count_prefix_mismatches(Tstring1 const & s1, Tstring2 const & s2);

//! Same as count_prefix_mismatches but reverses \a s2.
template <typename Tstring1, typename Tstring2>
int count_prefix_mismatches_s2rev(Tstring1 const & s1, Tstring2 const & s2);

template <typename Tstring>
std::vector<std::string_view> split_string(Tstring const & s1, char const delimiter, int const max_occurance = -1);

//! Find the MEM between the read and reference
void find_mem(int & mem_read_start, //
              int & mem_ref_start,
              int & mem_length,
              std::string const & read,
              std::string const & ref);

/*!
 * @brief Get the size of a completely matching prefix sequences
 *
 * @details
 * The prefix size can never be larger than either s1 or s2.
 */
int get_prefix_size(std::vector<char> const & s1, std::vector<char> const & s2);

/*!
 * @brief Get the size of a completely matching suffix sequences
 *
 * @see get_prefix_size()
 */
int get_suffix_size(std::vector<char> const & s1, std::vector<char> const & s2);

} // namespace weaver
