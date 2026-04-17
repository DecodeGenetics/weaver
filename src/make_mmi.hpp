#pragma once

#include <string>

#include "gfa.hpp"
#include "mmi.hpp"

namespace weaver
{
/*!
 * @brief Creates a minimizer index (MMI) for the GFA / rGFA graph.
 *
 * @details
 * The MMI keys are the minimizers (also refered to as "sketches") and their associated value, which is a list of
 * locations in the graph.
 *
 * @param[in] gfa graph.
 * @param[in] vcf_fn VCF filename containing small variants in a phased VCF.
 * @param[in] k kmer size. k<=0 uses default (KMIN in constants.hpp).
 * @param[in] w minimizer window size.
 *
 * @returns The minimizer index.
 *
 * @see MMI the returned index type.
 */
MMI make_mmi_index(GFA const & gfa, std::string const & vcf_fn = "", int k = -1, int w = -1);

} // namespace weaver
