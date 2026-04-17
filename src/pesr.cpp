/*!
 * @file pesr.cpp
 * @brief Contains functions for the \c pesr (paired-end short reads) subcommand.
 */
#include "pesr.hpp"

#include <iostream> // std::cout
#include <string>   // std::string

#include "edit_stats.hpp" // EditStats
#include "gfa.hpp"
#include "gfa_location.hpp"
#include "hashmap.hpp"
#include "index_io.hpp"   // deserialize_index()
#include "logging.hpp"    // print_info
#include "options.hpp"    // Options
#include "region.hpp"     // Region
#include "sam_writer.hpp" // SAMWriter
#include "sr_align.hpp"   // align_fastq

namespace weaver
{
void pesr(std::string const & graph_fn, //
          std::string const & vcf_fn,
          std::string const & fastq1,
          std::string const & fastq2,
          int k,
          int w)
{
  if (fastq1.empty())
  {
    print_error(_HERE_, " missing input FASTQ.");
    return;
  }

  if (fastq1 == fastq2)
    print_info("Mapping ", fastq1, " to ", graph_fn);
  else
    print_info("Mapping ", fastq1, " and ", fastq2, " to ", graph_fn);

  log_singleton->flush_stream();
  GFA gfa(graph_fn);
  print_info("Done reading gfa.");
  log_singleton->flush_stream();

  parse_target_region_option(gfa);

  MMI mmi;   // minimizer index
  T_icu icu; // "I see you" index

  print_info("Reading index '", graph_fn, ".wmi'");
  bool const is_success = deserialize_index(graph_fn, mmi, icu, k, w);

  if (is_success)
  {
    print_info("Done reading index with k = ", k, " w = ", w);
  }
  else
  {
    print_info("Did not find an index to read, creating one with k=", k, ", w=", w);
    log_singleton->flush_stream();
    deserialize_or_make_index(graph_fn, gfa, mmi, icu, vcf_fn, k, w);
  }

  SAMWriter sam_writer(std::cout);
  sam_writer.write_header(gfa);

#ifndef NDEBUG
  for (auto it = icu.begin(); it != icu.end(); ++it)
  {
    print_debug(_HERE_, " ", it->first >> 32, "|", (it->first << 32) >> 32, " ", it->second);
  }
#endif // NDEBUG

  // start aligning the reads in the fastqs
  print_info("Align starts");
  sr_align(gfa, mmi, icu, sam_writer, fastq1, fastq2, k, w);
  print_info("Align ends");
}

} // namespace weaver
