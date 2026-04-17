#pragma once
/*!
 * @file index.hpp
 * @brief Declares the function for the weaver's index subcommand.
 *
 * @see index_io.hpp for I/O operations using the index.
 */

#include <string>

#include "gfa.hpp"
#include "icu.hpp"
#include "mmi.hpp"

namespace weaver
{
/*! @brief Make an index that contains both minimizer and ICU index.
 *
 * @param[in]  gfa GFA graph.
 * @param[out] mmi The minimizer index.
 * @param[out] icu The "I see you" index.
 * @param[in]  vcf_fn Filename path to a VCF file. Empty if no VCF.
 * @param[in]  k minimizer index kmer size to use in the index.
 * @param[in]  w minimizer index window size.
 */
void make_index(GFA const & gfa, MMI & mmi, T_icu & icu, std::string const & vcf_fn, int k, int w);

} // namespace weaver
