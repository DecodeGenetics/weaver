#pragma once
/*!
 * @file sr_sam.hpp
 * @brief Defines methods for handling sam records with short reads.
 */
#include <string>
#include <vector>

#include "gfa.hpp"
#include "icu.hpp"

namespace weaver
{
// Forward declarations
class SRAlignment;
class SRSeed;
class SAMRecord;
class SAMWriter;
class HaplotypeStats;

//! Return true iff the seed and alignment extension does not cross any arcs of another stable sequence.
bool is_alignment_on_same_stable_sequence(GFA const & gfa, SRSeed const & seed, SRAlignment const & alignment);

//! Return true iff the seed does not cross any arcs which are not part of the same stable sequence.
bool is_seed_on_same_stable_sequence(GFA const & gfa, SRSeed const & seed);

/*!
 * @brief Go through the SAM record cigar and left align indels.
 *
 * @param[in] main_record The record to check the cigar of.
 * @param[in] full_read The full query read in the record in its aligned orientation.
 * @param[in] full_ref The full reference which the record aligned to.
 */
void left_align_record(SAMRecord & main_record, std::string const & full_read, std::string const & full_ref);

/*!
 * @brief Generate a SAM record from an alignment.
 *
 * @param[in] gfa Graph.
 * @param[in] hap_stats Haplotype statistics.
 * @param[in,out] sam_record The generated SAM record.
 * @param[in] alignment The alignment to create a SAM record from.
 * @param[in] seq_fwd The read sequence in forward orientation.
 * @param[in] seq_rev The read sequence in reverse orientation.
 *
 * @returns The edits between the alignment and the reference.
 */
void process_sam_record(GFA const & gfa, //
                        HaplotypeStats const & hap_stats,
                        SAMRecord & sam_record,
                        SRSeed const & seed,
                        SRAlignment const & alignment,
                        std::string const & seq_fwd,
                        std::string const & seq_rev);

//! Generate a pair of SAM records from alignments.
void process_sam_record_pair(GFA const & gfa, //
                             T_icu const & icu,
                             SAMRecord & read1,
                             SAMRecord & read2,
                             SRSeed const & seed1,
                             SRSeed const & seed2);

/*!
 * @brief Prepare a SAMRecord to be written.
 *
 * @param[in] gfa Graph.
 * @param[in] icu The "I see you" index.
 * @param[in] seed1 First record seed.
 * @param[in,out] record1 First record.
 * @param[in] seed2 Second record seed.
 * @param[in,out] record2 Second record.
 */
void prepare_sam_write(GFA const & gfa,
                       T_icu const & icu,
                       SRSeed const & seed1,
                       SAMRecord & record1,
                       SRSeed const & seed2,
                       SAMRecord & record2);

/*!
 * @brief Write a pair of SAMRecord to \a sam_writer .
 *
 * @param[out] sam_writer Write data into this writer.
 * @param[in] record1 First record to write.
 * @param[in] record2 Second record to write.
 */
void sam_write(SAMWriter & sam_writer, //
               SAMRecord & record1,
               SAMRecord & record2);

/*!
 * @brief Append the SAM lines to \a bucket_of_sam_lines
 *
 * @details
 * The bucket id will be determined and the lines are appending to that bucket.
 *
 * @param[out] bucket_of_sam_lines The two records will be written into one of these buckets.
 * @param[in] record1 First record to write.
 * @param[in] record2 Second record to write.
 */
void append_sam_lines(std::vector<std::vector<std::pair<uint64_t, std::string>>> & bucket_of_sam_lines,
                      SAMRecord & record1,
                      SAMRecord & record2);

} // namespace weaver
