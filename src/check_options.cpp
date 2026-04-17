#include "check_options.hpp"

/*!
 * @file check_options.cpp
 * @brief Implements functions for checking user defined options.
 *
 * @see options.hpp contains the object to strore the options.
 */

#include "logging.hpp"

namespace weaver
{
void is_k_option_valid(bool & is_valid, int const k)
{
  if (k > 0 && (((k % 2) == 0) || k < 17))
  {
    print_error("Invalid k=", k, ", only odd k>=17 supported.");
    is_valid = false;
  }
}

void is_w_option_valid(bool & is_valid, int const w)
{
  if (w > 0)
  {
    if (w == 1)
    {
      print_error("w==1 is not supported");
      is_valid = false;
    }
    else if (w > static_cast<int>(MAX_W))
    {
      print_error("w too large. Maximum w currently is ",
                  MAX_W,
                  ". Change 'MAX_W' in in.constants.hpp and recompile weaver to be able to use w=",
                  w);
      is_valid = false;
    }
  }
}

void is_read_group_header_line_valid(bool & is_valid, std::string const & read_group_header_line)
{
  if (read_group_header_line.size() > 0)
  {
    if (read_group_header_line.size() <= 8)
    {
      print_error(_HERE_, " Bad read_group_header_line=", read_group_header_line, " is too small.");
      print_error(_HERE_, " Should on the format: @RG\\tID:foo\tSM:bar");
      is_valid = false;
    }

    std::string const check_str = "@RG\\tID:";
    int cmp_val = read_group_header_line.compare(0, 8, check_str, 0, 8);

    if (cmp_val != 0)
    {
      print_error(_HERE_, " Bad read_group_header_line='", read_group_header_line, "' != '", check_str, "'");
      print_error(_HERE_, " returned comparison value=", cmp_val);
      is_valid = false;
    }

    if (read_group_header_line.find("\\tSM:") == std::string::npos)
    {
      print_error(_HERE_, " Bad read_group_header_line='", read_group_header_line, " that has no SM: tag");
      is_valid = false;
    }
  }
}

void are_fastq_options_valid(bool & is_valid, std::string const & fastq1, std::string const & fastq2)
{
  if (fastq1 == fastq2)
  {
    print_error(" Both FASTQs are the same path. Skip -2,--fq2 option if the file is interleaved.");
    is_valid = false;
  }
}

} // namespace weaver
