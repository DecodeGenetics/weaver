#pragma once
/*!
 * @file pesr.hpp
 * @brief Defines the functions for the \c pesr subcommand.
 */

#include <string>

namespace weaver
{
void pesr(std::string const & graph_fn, //
          std::string const & vcf_fn,
          std::string const & fastq1,
          std::string const & fastq2,
          int const k,
          int const w);

} // namespace weaver
