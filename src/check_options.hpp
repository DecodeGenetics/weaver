#pragma once

/*!
 * @file check_options.hpp
 * @brief Defines functions for checking user defined options.
 *
 * @see options.hpp contains the object to strore the options.
 */
#include <string> // std::string

namespace weaver
{
/*! @brief Checks if the options for \a k is valid.
 *
 * @param[out] is_valid Set to false if \a k is not valid.
 * @param[in]  k kmer size to check.
 */
void is_k_option_valid(bool & is_valid, int const k);

/*! @brief Checks if the options for \a w is valid.
 *
 * @param[out] is_valid Set to false if \a w is not valid.
 * @param[in]  w window size to check.
 */
void is_w_option_valid(bool & is_valid, int const w);

/*! @brief Checks if the the \a read_group_header_line is valid.
 *
 * @param[out] is_valid Set to false if \a read_group_header_line is not valid.
 * @param[in] read_group_header_line Read group header line to check.
 */
void is_read_group_header_line_valid(bool & is_valid, std::string const & read_group_header_line);

/*! @brief Checks if the options for \a fastq1 and \a fastq2 are valid.
 *
 * @param[out] is_valid set to false if \a fastq1 and \a fastq2 are not valid.
 * @param[in]  fastq1 FASTQ[.gz]1 file path to check.
 * @param[in]  fastq2 FASTQ[.gz]2 file path to check.
 */
void are_fastq_options_valid(bool & is_valid, std::string const & fastq1, std::string const & fastq2);

} // namespace weaver
