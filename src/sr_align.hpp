#pragma once

#include <string>

#include "fastq_data.hpp"
#include "gfa.hpp"
#include "hashmap_fwd_decl.hpp"
#include "icu.hpp"
#include "mmi.hpp"
#include "sam_writer.hpp"

namespace weaver
{
class HaplotypeStats;

/*!
 * @brief Align paired-end short reads in files \a fastq1 and \a fastq2 to a \a gfa graph.
 *
 * @param[in] gfa GFA input.
 * @param[in] mmi minimizer index.
 * @param[in] icu The "I see you" index.
 * @param[in] sam_writer Write alignment SAM records to this writer
 * @param[in] fastq1 First FASTQ filename. May be a gzipped file.
 * @param[in] fastq2 Second FASTQ filename. May be a gzipped file.
 * @param[in] k kmer size.
 * @param[in] w window size.
 *
 * @see align_pair
 */
void sr_align(GFA const & gfa,            // graph
              MMI const & mmi,            // index
              T_icu const & icu,          // The I see you index
              SAMWriter & sam_writer,     // Writer for alignments in SAM format
              std::string const & fastq1, // fastq1 filename
              std::string const & fastq2, // fastq2 filename
              int const k,                // kmer size
              int const w);               // window size

/*!
 * @brief Align a pair of reads in \a record1 and \a record2 to \a gfa graph.
 *
 * @param[in] is_unique_mapping Set iff the read is only mapped if it is mapped uniquely to adjust haplotype stats.
 * @param[in] gfa GFA input.
 * @param[in] mmi minimizer index.
 * @param[in] icu The "I see you" index.
 * @param[in,out] record1 Record containing the first read, output is written to the record as well.
 * @param[in,out] record2 Record containing the second read, output is written to the record as well.
 * @param[in] k kmer size.
 * @param[in] w window size.
 * @param[in] s1 Sequence of the first read.
 * @param[in] s2 Sequence of the second read.
 * @param[in] is_debug Set to add additionally debugging messages.
 *
 * @returns True iff read should be written.
 *
 * @see sr_align() for aligning all read pairs from FASTQ data.
 */
bool sr_align_pair(bool const is_unique_mapping,
                   GFA const & gfa,
                   MMI const & mmi,
                   T_icu const & icu,
                   HaplotypeStats & hap_stats,
                   SAMRecord & record1,
                   SAMRecord & record2,
                   int const k,
                   int const w,
                   std::string const & s1,
                   std::string const & s2,
                   bool const is_debug);
} // namespace weaver
