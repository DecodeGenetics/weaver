#include "sr_align.hpp"

#include <algorithm> // find_if
#include <cassert>   // assert
#include <iterator>  // std::advance
#include <kseq.h>    // gzread
#include <numeric>
#include <sstream> // std::ostringstream
#include <string>  // std::string
#include <vector>
#include <zlib.h> // gzFile

#include <paw/station.hpp>

#include "adapters.hpp" // add_adapter_sketches
#include "alignment_utils.hpp"
#include "fastq_data.hpp"
#include "filesystem.hpp"
#include "gfa.hpp"             // GFA
#include "haplotype_stats.hpp" // HaplotypeStats
#include "hashmap.hpp"         // weaver::Tset
#include "io.hpp"
#include "logging.hpp" // print_info
#include "options.hpp"
#include "read_sketch.hpp"
#include "region.hpp"
#include "sam_record.hpp"
#include "sam_writer.hpp"
#include "sam_zstd_reader.hpp"
#include "sam_zstd_writer.hpp"
#include "sequence_utils.hpp"
#include "sketch.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp"
#include "sr_alignment.hpp"
#include "sr_alignment_sel.hpp"
#include "sr_chain.hpp"
#include "sr_sam.hpp"
#include "sr_seed.hpp"
#include "sr_seed_pair.hpp"
#include "stable_contigs.hpp"
#include "system.hpp"

namespace
{
static void destroy_kseq_t(kseq_t * s)
{
  if (s != nullptr)
    kseq_destroy(s);
}

} // namespace

namespace weaver
{
bool sr_align_pair(bool const,      //  is_unique_mapping
                   GFA const & gfa, //
                   MMI const & mmi,
                   T_icu const & icu,
                   HaplotypeStats & hap_stats,
                   SAMRecord & sam_record1,
                   SAMRecord & sam_record2,
                   int const k,
                   int const w,
                   std::string const & seq1,
                   std::string const & seq2,
                   bool const is_debug)
{
  int const num_segments = gfa.get_num_segments();

  // placeholders
  SRSeed empty_seed;
  SRAlignment empty_alignment;
  SAMRecord empty_record1(sam_record1);
  SAMRecord empty_record2(sam_record2);

  std::vector<SRSeed> seeds1;
  std::vector<SRSeed> seeds2;

  {
    // get the read sketches
    std::vector<ReadSketch> read_sketches1 = get_read_sketches(seq1, num_segments, mmi, k, w);
    std::vector<ReadSketch> read_sketches2 = get_read_sketches(seq2, num_segments, mmi, k, w);

    if (is_debug)
    {
      print_info("sketches before adapter removal: ", read_sketches1.size(), " ", read_sketches2.size());
    }

    // check for matches between reads and ignores adapters
    remove_adapter_sketches(read_sketches1, read_sketches2, sam_record1.qname);

    if (is_debug)
    {
      print_info("sketches after adapter removal: ", read_sketches1.size(), " ", read_sketches2.size());
    }

    get_seeds_with_pair(seeds1, seeds2, read_sketches1, read_sketches2, gfa, icu, is_debug);
  }

  if (is_debug)
  {
    print_info(_HERE_, " num seeds read1, read2=", seeds1.size(), ", ", seeds2.size());

    auto order1 = get_seed_est_score_sorted_order_indices(seeds1);
    auto order2 = get_seed_est_score_sorted_order_indices(seeds2);

    for (int o1{0}; o1 < std::min(128, static_cast<int>(order1.size())); ++o1)
      print_info("[all] seed1: ", order1[o1], " ", seeds1[order1[o1]].to_string());

    for (int o2{0}; o2 < std::min(128, static_cast<int>(order2.size())); ++o2)
      print_info("[all] seed2: ", order2[o2], " ", seeds2[order2[o2]].to_string());
  }

  get_sr_chains_and_filter(gfa, icu, seeds1, seeds2);

  if (is_debug)
    print_info(_HERE_, " num shrinked seeds read1, read2=", seeds1.size(), ", ", seeds2.size());

  Options const & copts = *(Options::const_instance());

  // Stop if outside target region
  if (copts.region_sfa_idx >= 0)
  {
    if (are_all_seeds_outside_region(gfa, seeds1) && are_all_seeds_outside_region(gfa, seeds2))
    {
      print_debug(_HERE_, " skipped writing a read");
      return false; // Not in target region, skip writing read
    }
  }
  else if (seeds1.empty() && seeds2.empty())
  {
    prepare_sam_write(gfa, icu, empty_seed, sam_record1, empty_seed, sam_record2);
    return true; // No seeds found, no need to try a remap
  }

  // sequence 1 and 2 reversed
  std::string const & seq1_rev = get_reverse_complement(seq1);
  std::string const & seq2_rev = get_reverse_complement(seq2);

  std::vector<SRAlignment> alignments1(seeds1.size());

  for (size_t s{0}; s < seeds1.size(); ++s)
    alignments1[s] = extend_and_align(seeds1[s], seq1, seq1_rev);

  std::vector<SRAlignment> alignments2(seeds2.size());

  for (size_t s{0}; s < seeds2.size(); ++s)
    alignments2[s] = extend_and_align(seeds2[s], seq2, seq2_rev);

  // remove seeds with no score
  remove_seeds_with_no_score(seeds1, alignments1);
  remove_seeds_with_no_score(seeds2, alignments2);

  // Must be after extend and align
  remove_duplicates(gfa, icu, seeds1, alignments1);
  remove_duplicates(gfa, icu, seeds2, alignments2);

  if (seeds1.empty() && seeds2.empty())
  {
    prepare_sam_write(gfa, icu, empty_seed, sam_record1, empty_seed, sam_record2);
    return true; // No seeds found, no need to try a remap
  }

  if (is_debug)
  {
    print_info(_HERE_, " num seeds after duplicate rm read1, read2=", seeds1.size(), ", ", seeds2.size());

    for (int o1{0}; o1 < static_cast<int>(seeds1.size()); ++o1)
      print_info("seed1: ", o1, " ", seeds1[o1].to_string());

    print_info("");

    for (int o2{0}; o2 < static_cast<int>(seeds2.size()); ++o2)
      print_info("seed2: ", o2, " ", seeds2[o2].to_string());
  }

  std::vector<SAMRecord> sam_records1(seeds1.size(), empty_record1);
  std::vector<SAMRecord> sam_records2(seeds2.size(), empty_record2);

  // Generate a CIGAR string
  for (int s1{0}; s1 < static_cast<int>(seeds1.size()); ++s1)
  {
    SRSeed const & seed1 = seeds1[s1];
    SRAlignment & alignment1 = alignments1[s1];
    SAMRecord & record1 = sam_records1[s1];
    process_sam_record(gfa, hap_stats, record1, seed1, alignment1, seq1, seq1_rev);
    assert(record1.weaver_score != std::numeric_limits<int>::lowest());
    alignment1.score = record1.weaver_score;
  }

  for (int s2{0}; s2 < static_cast<int>(seeds2.size()); ++s2)
  {
    SRSeed const & seed2 = seeds2[s2];
    SRAlignment & alignment2 = alignments2[s2];
    SAMRecord & record2 = sam_records2[s2];
    process_sam_record(gfa, hap_stats, record2, seed2, alignment2, seq2, seq2_rev);
    assert(record2.weaver_score != std::numeric_limits<int>::lowest());
    alignment2.score = record2.weaver_score;
  }

  SRAlignmentSel sas = select_alignment(gfa, icu, seeds1, sam_records1, seeds2, sam_records2);

  SRSeed & seed1 = get_element_reference(sas.s1, seeds1, empty_seed);
  SRSeed & seed2 = get_element_reference(sas.s2, seeds2, empty_seed);

  // Stop if outside region
  if (copts.region_sfa_idx >= 0)
  {
    bool is_in_region = is_seed_in_region(gfa, seed1) || is_seed_in_region(gfa, seed2);

    if (!is_in_region)
    {
      // print_info(_HERE_, " skipping read because it's outside the target region.");
      return false;
    }
  }

  if ((seed1.is_empty() || is_seed_on_same_stable_sequence(gfa, seed1)) &&
      (seed2.is_empty() || is_seed_on_same_stable_sequence(gfa, seed2)))
  {
    SRAlignment & alignment1 = get_element_reference(sas.s1, alignments1, empty_alignment);
    SRAlignment & alignment2 = get_element_reference(sas.s2, alignments2, empty_alignment);

    // Make the output records reference the primary one
    sam_record1 = get_element_reference(sas.s1, sam_records1, empty_record1);
    sam_record2 = get_element_reference(sas.s2, sam_records2, empty_record2);

    // read 1 qc
    sam_record1.mapq = sas.get_mapq1(alignment1.score, /*read1_length=*/static_cast<int>(seq1.size()));
    // sam_record1.weaver_score = alignment1.score;

    // read 2 qc
    sam_record2.mapq = sas.get_mapq2(alignment2.score, /*read2_length=*/static_cast<int>(seq2.size()));
    // sam_record2.weaver_score = alignment2.score;

    assert(SAMRecord::MISSING_TAG == SRAlignmentSel::MISSING_SCORE);

    if (sas.other_read1_best_read_score != SRAlignmentSel::MISSING_SCORE)
      sam_record1.second_alignment_score = sas.other_read1_best_read_score;

    if (sas.other_read2_best_read_score != SRAlignmentSel::MISSING_SCORE)
      sam_record2.second_alignment_score = sas.other_read2_best_read_score;

    prepare_sam_write(gfa, icu, seed1, sam_record1, seed2, sam_record2);

    // filtering bad alignments
    if (sam_record1.is_double_clipped() && sam_record2.is_double_clipped())
    {
      print_debug("Making a record pair unmapped due to double clipping.");
      sam_record1.make_unmapped(/*other_sam_record=*/sam_record2);
      sam_record2.make_unmapped(/*other_sam_record=*/sam_record1);
    }

    if (sam_record1.alignment_score < copts.mininum_score_to_output)
      sam_record1.make_unmapped(sam_record2);

    if (sam_record2.alignment_score < copts.mininum_score_to_output)
      sam_record2.make_unmapped(sam_record1);

    // problem checking
    if (not sam_record1.is_cigar_valid())
    {
      print_warning(_HERE_, " read 1 with invalid cigar, qname=", sam_record1.qname);
      sam_record1.make_unmapped(/*other_sam_record=*/sam_record2);
    }

    if (not sam_record2.is_cigar_valid())
    {
      print_warning(_HERE_, " read 2 with invalid cigar, qname=", sam_record2.qname);
      sam_record2.make_unmapped(/*other_sam_record=*/sam_record1);
    }

    if (sam_record1.is_unmapped_with_mapped_mate())
    {
      sam_record1.sfa_idx = sam_record2.sfa_idx;
      sam_record1.pos = sam_record2.pos;
    }
    else if (sam_record2.is_unmapped_with_mapped_mate())
    {
      sam_record2.sfa_idx = sam_record1.sfa_idx;
      sam_record2.pos = sam_record1.pos;
    }

    assert(sam_record1.is_pair_valid(sam_record2));
    assert(sam_record2.is_pair_valid(sam_record1));
  }

  return true;
}

void parallel_printer_st(SAMWriter & sam_writer, std::vector<filesystem::path> const & sam_zstd_paths)
{
  int const num_buckets = stable_contigs.get_num_buckets();

  {
    ParallelSAMZstdReader preader(sam_zstd_paths.size());
    print_debug(_HERE_, " writing directly");

    for (int b{0}; b < num_buckets; ++b)
    {
      for (int i{0}; i < static_cast<int>(sam_zstd_paths.size()); ++i)
        preader.open(sam_zstd_paths[i] / std::to_string(b), i);

      std::string sam_line;

      while (preader.read_next_line(sam_line))
        sam_writer.write_line(sam_line);
    }
  }

  // Last bucket contains only unmapped reads, no need to use the parallel reader
  // print_info(_HERE_, " writing directly last bucket st. b=", num_buckets);
  SAMZstdReader reader(/*file_index=*/0); // opens a new zstd compressed sam file

  for (filesystem::path const & sam_zstd_path : sam_zstd_paths)
  {
    filesystem::path const p = sam_zstd_path / std::to_string(num_buckets);

    if (filesystem::exists(p))
    {
      reader.open(p);
      SAMOrderLine const * ol = reader.read_line();

      while (ol != nullptr)
      {
        sam_writer.write_line(ol->sam_line);
        ol = reader.read_line();
      }
    }
  }

  print_debug(_HERE_, " DONE writing directly. last bucket b=", num_buckets);
}

void parallel_align(int const thread_id,
                    GFA const * gfa_ptr,
                    MMI const * mmi_ptr,
                    T_icu const * icu_ptr,
                    std::vector<std::string> const * sam_lines_fn_ptr,
                    std::vector<uint32_t> * sam_lines_counter_ptr,
                    std::vector<HaplotypeStats> * p_hap_stats_ptr,
                    std::vector<FastqData> * reads_ptr,
                    int const k,
                    int const w,
                    bool const is_unique_mapping)
{
  assert(gfa_ptr);
  assert(mmi_ptr);
  assert(icu_ptr);
  assert(sam_lines_fn_ptr);
  assert(thread_id < static_cast<int>(sam_lines_fn_ptr->size()));
  assert(sam_lines_counter_ptr);
  assert(thread_id < static_cast<int>(sam_lines_counter_ptr->size()));
  assert(p_hap_stats_ptr);
  assert(thread_id < static_cast<int>(p_hap_stats_ptr->size()));
  assert(reads_ptr);

  // in this struct we store a full SAM line and its order value for sorting
  using Tsam_line = std::pair<uint64_t, std::string>;

  // Get references from the pointers
  GFA const & gfa = *gfa_ptr;
  MMI const & mmi = *mmi_ptr;
  T_icu const & icu = *icu_ptr;
  std::string const & sam_lines_fn = (*sam_lines_fn_ptr)[thread_id];
  uint32_t & sam_lines_counter = (*sam_lines_counter_ptr)[thread_id];
  HaplotypeStats & hap_stats = (*p_hap_stats_ptr)[thread_id];

  int const num_buckets = stable_contigs.get_num_buckets();
  std::vector<std::vector<Tsam_line>> bucket_of_sam_lines;
  bucket_of_sam_lines.resize(num_buckets + 1); // the "+ 1" is for the bucket with unmapped reads
  int const num_reads = reads_ptr->size();
  std::vector<FastqData> && reads = std::move(*reads_ptr);
  bool is_read_name_warning_printed{false};

  // Unique mapping
  for (int i{1}; i < num_reads; i += 2)
  {
    Options const & copts = *(Options::const_instance());
    FastqData & read_data1 = reads[i - 1];
    FastqData & read_data2 = reads[i];

    read_data1.remove_slash_from_name();
    read_data2.remove_slash_from_name();

    if (copts.rta3_quals)
    {
      read_data1.rta3_bin_quals();
      read_data2.rta3_bin_quals();
    }

    assert(read_data1.name == read_data2.name);

    if (!is_read_name_warning_printed && read_data1.name != read_data2.name)
    {
      print_warning("Found a pair of reads with different read names. Please check your FASTQ input files.");
      print_warning("Read 1 name = ", read_data1.name);
      print_warning("Read 2 name = ", read_data2.name);
      is_read_name_warning_printed = true;
    }

    SAMRecord sam_record1(read_data1);
    SAMRecord sam_record2(read_data2);

    bool const is_debug = copts.debug_read_name.size() > 0;

    if (is_debug)
    {
      if (read_data1.name != copts.debug_read_name)
        continue;

      print_info(_HERE_, " debug read found=", copts.debug_read_name);
    }

    bool const is_writing = sr_align_pair(is_unique_mapping,
                                          gfa,
                                          mmi,
                                          icu,
                                          hap_stats,
                                          sam_record1,
                                          sam_record2,
                                          k,
                                          w,
                                          read_data1.seq,
                                          read_data2.seq,
                                          is_debug);

    if (is_writing)
      append_sam_lines(bucket_of_sam_lines, sam_record1, sam_record2);
  }

  delete reads_ptr; // Read data has been used, delete here to free memory

  filesystem::path path = sam_lines_fn + std::to_string(sam_lines_counter);
  ++sam_lines_counter;

  if (filesystem::exists(path))
    print_warning("Path=', path, ' exists, overwriting.");

  filesystem::create_directory(path);
  auto first_bucket_it = std::find_if_not(bucket_of_sam_lines.begin(),
                                          bucket_of_sam_lines.end(),
                                          [](std::vector<Tsam_line> const & sl) { return sl.empty(); });

  if (first_bucket_it == bucket_of_sam_lines.end())
    return;

  int64_t const first_bucket_index = std::distance(bucket_of_sam_lines.begin(), first_bucket_it);
  SAMZstdWriter sam_zstd_writer(path / std::to_string(first_bucket_index));

  // b is the bucket index
  for (int64_t b{first_bucket_index}; b < static_cast<int64_t>(bucket_of_sam_lines.size()); ++b)
  {
    std::vector<Tsam_line> & sam_lines = bucket_of_sam_lines[b];

    // Go to next bucket if this one is empty
    if (sam_lines.empty())
      continue;

    if (b > first_bucket_index)
    {
      // Reopen the file with a new path
      sam_zstd_writer.close();
      sam_zstd_writer.open(path / std::to_string(b));
    }

    // Sort the sam lines before printing, unless this is the bucket with unmapped reads
    bool const is_bucket_with_unmapped_reads = (b + 1ll) == static_cast<int64_t>(bucket_of_sam_lines.size());

    if (not is_bucket_with_unmapped_reads)
    {
      std::sort(sam_lines.begin(),
                sam_lines.end(),
                [](Tsam_line const & line1, Tsam_line const & line2) { return line1.first < line2.first; });
    }

    // Print to file
    for (Tsam_line const & sam_line : sam_lines)
      sam_zstd_writer.write_line(sam_line);
  }

  print_info(_HERE_, " DONE writing to sam_lines_fn=", path);
}

void sr_align(GFA const & gfa,            // graph
              MMI const & mmi,            // index
              T_icu const & icu,          // The I see you index
              SAMWriter & sam_writer,     // Writer for alignments in SAM format
              std::string const & fastq1, // fastq1 filename
              std::string const & fastq2, // fastq2 filename
              int const k,                // kmer size
              int const w)                // window size
{
  using kseq_t_ptr = std::unique_ptr<kseq_t, void (*)(kseq_t *)>; //! type definition of a smart kseq_t pointer

  gzFile fp1 = fastq1 == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fastq1.c_str(), "r");

  if (fp1 == nullptr)
  {
    print_error("Could not open file FASTQ file: \"", fastq1, '"');
    std::exit(1);
  }

  gzFile fp2 = nullptr;
  bool const is_interleaved = fastq1 == fastq2;

  if (not is_interleaved)
  {
    fp2 = gzopen(fastq2.c_str(), "r");

    if (fp2 == nullptr)
    {
      print_error("Could not open file FASTQ file 2: \"", fastq2, '"');
      std::exit(1);
    }
  }

  Options const & copts = *(Options::const_instance());
  int const threads = copts.threads <= 1 ? 1 : copts.threads;
  kseq_t_ptr seq1 = kseq_t_ptr(kseq_init(fp1), ::destroy_kseq_t);
  kseq_t_ptr seq2 = kseq_t_ptr(is_interleaved ? nullptr : kseq_init(fp2), ::destroy_kseq_t);

  std::string tmp_dir = create_temp_dir();
  std::vector<uint32_t> sam_lines_counter(threads);
  std::vector<std::string> sam_lines_fn(threads, tmp_dir);

  for (int t{0}; t < threads; ++t)
    sam_lines_fn[t] += "/thread" + std::to_string(t) + ".chunk";

  std::vector<HaplotypeStats> p_hap_stats(threads); // Create a hap stats object for each thread
  bool constexpr IS_USING_HMM{false};

  for (auto & hap_stats : p_hap_stats)
    hap_stats.set_haps_ptr(&mmi.haplotypes);

  {
    paw::Station read_station(threads, /*max_queue_size=*/3);
    int num_read_buffers{0};

    while (num_read_buffers < copts.max_read_buffers)
    {
      std::vector<FastqData> * seqs{nullptr};

      if (is_interleaved)
        seqs = seq_store_get(seq1.get());
      else
        seqs = seq_store_get(seq1.get(), seq2.get());

      assert(seqs != nullptr);

      if (seqs->empty())
      {
        delete seqs; // free the last buffer
        break;       // stop here, we have seen the last buffer
      }

      ++num_read_buffers;

      if (num_read_buffers == copts.max_read_buffers || static_cast<int>(seqs->size()) < copts.fastq_data_buffer_size)
      {
        // should be the last buffer, make the boss work on it
        read_station.add_to_thread_with_thread_id(static_cast<std::size_t>(threads - 1),
                                                  parallel_align,
                                                  &gfa,
                                                  &mmi,
                                                  &icu,
                                                  &sam_lines_fn,
                                                  &sam_lines_counter,
                                                  &p_hap_stats,
                                                  seqs,
                                                  k,
                                                  w,
                                                  /*is_unique_mapping=*/IS_USING_HMM);
      }
      else
      {
        // the common case, we have a full buffer
        read_station.add_work_with_thread_id(parallel_align, //
                                             &gfa,
                                             &mmi,
                                             &icu,
                                             &sam_lines_fn,
                                             &sam_lines_counter,
                                             &p_hap_stats,
                                             seqs,
                                             k,
                                             w,
                                             /*is_unique_mapping=*/IS_USING_HMM);
      }
    }

    std::string thread_info = read_station.join();
    print_info(_HERE_, " finished mapping, thread info=", thread_info);
  }

  // Get a path to every directory containing zstd compressed SAM files
  std::vector<filesystem::path> const sam_zstd_paths = get_sam_zstd_paths(sam_lines_fn, sam_lines_counter);

  // Write SAM data
  if (sam_zstd_paths.size() > 0)
  {
    parallel_printer_st(sam_writer, sam_zstd_paths);
    sam_writer.flush(); // flush the remaining data
    print_info(_HERE_, " All SAM data written.");
  }

  // Remove temporary files if there is no "no-cleanup" flag set
  if (copts.no_cleanup)
  {
    print_info("Temporary files left in ", tmp_dir);
  }
  else
  {
    print_info("Cleaning up files in ", tmp_dir);
    filesystem::remove_all(tmp_dir);
  }

  if (fp1 != nullptr)
    gzclose(fp1);

  if (fp2 != nullptr)
    gzclose(fp2);
}

} // namespace weaver
