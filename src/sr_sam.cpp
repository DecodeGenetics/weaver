/*!
 * @file sr_sam.cpp
 * @brief Implements methods for handling sam records with short reads.
 */

#include "sr_sam.hpp"

#include <cstdint>     // uint32_t
#include <iterator>    // std::next, std::advance
#include <string>      // std::string
#include <string_view> // std::string_view

#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/pairwise_alignment.hpp>

#include "alignment_utils.hpp"  // add_log
#include "edit.hpp"             // Variant
#include "edit_stats.hpp"       // VariantStats
#include "gfa.hpp"              // GFA
#include "gfa_location.hpp"     // GFALocation
#include "haplotype_stats.hpp"  // HaplotypeStats
#include "logging.hpp"          // print_debug
#include "options.hpp"          // Options
#include "sam_record.hpp"       // SAMRecord
#include "sam_writer.hpp"       // SAMWriter
#include "sequence_utils.hpp"   // get_reverse_complement, isACGT
#include "sketch_to_string.hpp" // sketch_value_to_string
#include "sketch_value.hpp"     // sketch_value_pos
#include "sr_alignment.hpp"     // SRAlignment
#include "sr_cigar.hpp"         // process_cigar
#include "sr_seed.hpp"          // SRSeed
#include "stable_contigs.hpp"   // stable_contigs

namespace
{
uint64_t get_edit_order(int pos, char type)
{
  if (pos < 0 || type == 'N')
    return std::numeric_limits<uint64_t>::max();
  else
    return (static_cast<uint64_t>(pos) << 2ull) | static_cast<uint64_t>((type == 'D') + 2 * (type == 'X'));
}

int edit_lower_bound(std::vector<weaver::EditCalls> const & edit_calls, int pos)
{
  int first{0};
  int count{static_cast<int>(edit_calls.size())};

  while (count > 0)
  {
    int step{count / 2};
    int i{first + step};

    if (edit_calls[i].pos < pos)
    {
      first = i + 1;
      count -= step + 1;
    }
    else
    {
      count = step;
    }
  }

  return first;
}

int edit_upper_bound(std::vector<weaver::EditCalls> const & edit_calls, int v_hi, int pos_end)
{
  // linear search for upper bound since it is almost never than a couple of checks away
  while (v_hi < static_cast<int>(edit_calls.size()) && edit_calls[v_hi].pos < pos_end)
    ++v_hi;

  return v_hi;
}

void advance_from_seed(int by, weaver::GFALocation & location, std::vector<gfa_arc_t const *> & arcs, std::string & seq)
{
  int unaccounted = location.advance_when_same_contig(1, arcs);

  if (unaccounted == 0)
  {
    unaccounted = location.advance_when_same_contig_and_get_sequence(by, arcs, seq);

    if (unaccounted > 0)
    {
      location.get_base(seq);
      seq += std::string(unaccounted - 1, 'X');
    }
  }
  else
  {
    seq += std::string(by, 'X');
  }
}

} // namespace

namespace weaver
{
bool is_seed_on_same_stable_sequence(GFA const & gfa, SRSeed const & seed)
{
  return gfa.are_arcs_with_same_rank_as_segments(seed.arcs); // seed arcs
}

bool is_alignment_on_same_stable_sequence(GFA const & gfa, SRSeed const & seed, SRAlignment const & alignment)
{
  return is_seed_on_same_stable_sequence(gfa, seed) &&                    // seed arcs
         gfa.are_arcs_with_same_rank_as_segments(alignment.begin_arcs) && // begin extension arcs
         gfa.are_arcs_with_same_rank_as_segments(alignment.end_arcs);     // end extension arcs
}

void process_sam_record_seq_reverse_flag(SAMRecord & main_record, SAMRecord & other_record, SRSeed const & seed)
{
  if (seed.is_empty())
    return;

  if (sketch_value_strand(seed.end_ref_value) == 1u)
  {
    // reversed
    main_record.flags |= SAMFlags::IS_SEQ_REVERSED;
    other_record.flags |= SAMFlags::IS_MATE_SEQ_REVERSED;
  }
}

void left_align_record(SAMRecord & main_record, std::string const & full_read, std::string const & full_ref)
{
  int const num_cigar_operatations{static_cast<int>(main_record.cig.size())};
  std::string::const_iterator ref_it = full_ref.cbegin();
  std::string::const_iterator read_it = full_read.cbegin();

  for (int c{0}; c < num_cigar_operatations; ++c)
  {
    paw::Cigar const & cigar = main_record.cig[c];
    assert(cigar.count > 0);

    switch (cigar.operation)
    {
    case paw::CigarOperation::MATCH:
    case paw::CigarOperation::EQUAL:
    case paw::CigarOperation::DIFFERENT:
    {
      std::advance(ref_it, cigar.count);
      std::advance(read_it, cigar.count);
      break;
    }

    case paw::CigarOperation::INSERTION:
    {
      left_align_insertion(main_record, c, ref_it, read_it);
      std::advance(read_it, cigar.count);
      break;
    }

    case paw::CigarOperation::DELETION:
    {
      left_align_deletion(main_record, c, ref_it, read_it);
      std::advance(ref_it, cigar.count);
      break;
    }

    case paw::CigarOperation::SOFT_CLIP:
    {
      std::advance(read_it, cigar.count);
      break;
    }

    case paw::CigarOperation::HARD_CLIP:
    {
      break;
    }

    default:
    {
      print_warning(_HERE_, " unexpected cigar operation: ", static_cast<int>(cigar.operation));
      break;
    }
    } // switch end

    assert(std::distance(ref_it, full_ref.cend()) >= 0);
    assert(std::distance(read_it, full_read.cend()) >= 0);
  } // for end

  assert(ref_it == full_ref.cend());
  assert(read_it == full_read.cend());
}

void calculate_alignment_score(HaplotypeStats const & hap_stats,
                               SAMRecord & sam_record,
                               std::string const & ref_seed,
                               std::string && ref_before,
                               std::string const & ref_after,
                               bool const is_strand_forward,
                               std::string const & seq_fwd,
                               std::string const & seq_rev)
{
  std::vector<int> edits;        // edits in the alignment, negative values v indicate no edit at for edit -v-1
  std::vector<int> nonsnp_edits; // nonsnp edits are handled separately

  assert(nonsnp_edits.size() == 0);
  assert(edits.size() == 0);

  print_debug(_HERE_, " pos=", sam_record.pos);
  print_debug(_HERE_, " ref_before=", ref_before);
  print_debug(_HERE_, " ref_seed=", ref_seed);
  print_debug(_HERE_, " ref_after=", ref_after);
  std::string full_ref = std::move(ref_before);

  if (is_strand_forward)
    full_ref += ref_seed;
  else
    full_ref += get_reverse_complement(ref_seed);

  full_ref += ref_after;

  std::string const & full_read = is_strand_forward ? seq_fwd : seq_rev;

  print_debug(_HERE_, " read_name=", sam_record.qname);
  print_debug(_HERE_, "  full_ref=", full_ref);
  print_debug(_HERE_, " full_read=", full_read);
  print_debug(_HERE_, "     cigar=", cigar2string(sam_record.cig.begin(), sam_record.cig.end()));
  assert(static_cast<int>(full_read.size()) == sam_record.get_cigar_query_length());
  assert(static_cast<int>(full_ref.size()) == sam_record.get_cigar_reference_length());
  assert(sam_record.is_cigar_valid());

  left_align_record(sam_record, full_read, full_ref);

  assert(sam_record.snid >= 0);
  std::vector<EditCalls> const * ec_ptr{nullptr};
  int e_lo{0}; // edit lower bound
  int e_hi{0}; // edit higher bound
  assert(hap_stats.haps_ptr != nullptr);
  MMI::T_haplotypes const & haps = *(hap_stats.haps_ptr);

  if (sam_record.snid < static_cast<int>(haps.size()))
  {
    assert(haps.size() > 0);
    assert(haps.size() == hap_stats.snid_haps_stats.size());
    ec_ptr = &haps[sam_record.snid];
    e_lo = edit_lower_bound(*ec_ptr, sam_record.pos);
    e_hi = edit_upper_bound(*ec_ptr, e_lo, sam_record.pos + static_cast<int>(full_ref.size()));

    // find upper bound with linear search
    print_debug(_HERE_, " sam_record.pos=", sam_record.pos, " e_lo, e_hi=", e_lo, ", ", e_hi);
  }

  // scan read again to calculate the score
  std::string::const_iterator ref_it = full_ref.cbegin();
  std::string::const_iterator read_it = full_read.cbegin();
  int ref_pos = sam_record.pos;
  int const num_cigar_operatations = sam_record.cig.size();
  sam_record.alignment_score = 0;
  sam_record.num_edits = 0;
  Options const & copts = *(Options::const_instance());
  // int n_edits_with_no_as_effect{0};
  // int ref_edit_matches{0};
  // int ref_edit_mismatches{0};

  auto get_next_edit_order_at_e_lo = [&ec_ptr, &e_lo, e_hi](uint64_t & next_edit_order)
  {
    if (e_lo < e_hi)
    {
      assert(ec_ptr != nullptr);
      EditCalls const & ne = (*ec_ptr)[e_lo]; // next edit
      assert(ne.pos >= 0);
      assert(ne.type != 'N');
      next_edit_order = ::get_edit_order(ne.pos, ne.type);
    }
    else
    {
      next_edit_order = std::numeric_limits<uint64_t>::max();
    }
  };

  uint64_t next_edit_order{std::numeric_limits<uint64_t>::max()};
  get_next_edit_order_at_e_lo(next_edit_order); // Gets the first edit

  auto update_next_edit = [&next_edit_order, &e_lo, &get_next_edit_order_at_e_lo](uint64_t ref_pos, char type)
  {
    uint64_t const current_pos_order = ::get_edit_order(ref_pos, type);

    while (next_edit_order < current_pos_order)
    {
      ++e_lo;
      get_next_edit_order_at_e_lo(next_edit_order);
    }
  };

  for (int c{0}; c < num_cigar_operatations; ++c)
  {
    paw::Cigar const & cigar = sam_record.cig[c];
    int const c_count = cigar.count;
    assert(c_count > 0);

    switch (cigar.operation)
    {
    case paw::CigarOperation::INSERTION:
    {
      int const score_penalty = copts.gap_open + copts.gap_extend * (c_count - 1);
      // print_info(_HERE_, " insertion at in cigar c = ", c, " ", c_count, " score=", sam_record.alignment_score);
      update_next_edit(ref_pos, 'I');

      while (next_edit_order == ::get_edit_order(ref_pos, 'I'))
      {
        assert(ec_ptr);
        assert(e_lo < static_cast<int>(ec_ptr->size()));
        assert(e_lo < e_hi);

        EditCalls const & edit_calls = (*ec_ptr)[e_lo];
        assert(edit_calls.pos == ref_pos);
        assert(edit_calls.type == 'I');

        if (static_cast<int>(edit_calls.seq.size()) == c_count &&
            std::equal(edit_calls.seq.begin(), edit_calls.seq.end(), read_it))
        {
          print_debug(_HERE_, " exact matched INS with ", edit_calls.to_string(), " index=", e_lo);
          nonsnp_edits.push_back(e_lo);
        }

        ++e_lo;
        get_next_edit_order_at_e_lo(next_edit_order); // advance to the next edit
      }

      sam_record.num_edits += c_count;
      assert(copts.gap_open >= 0);
      assert(copts.gap_extend >= 0);
      assert(c_count >= 1);
      sam_record.alignment_score -= score_penalty;
      std::advance(read_it, c_count);
      break;
    }

    case paw::CigarOperation::DELETION:
    {
      int const score_penalty = copts.gap_open + copts.gap_extend * (c_count - 1);
      // print_info(_HERE_, " deletion at in cigar c=", c;
      update_next_edit(ref_pos, 'D');

      while (next_edit_order == ::get_edit_order(ref_pos, 'D'))
      {
        assert(ec_ptr);
        assert(e_lo < static_cast<int>(ec_ptr->size()));
        assert(e_lo < e_hi);

        EditCalls const & edit_calls = (*ec_ptr)[e_lo];
        assert(edit_calls.pos == ref_pos);
        assert(edit_calls.type == 'D');

        if (static_cast<int>(edit_calls.seq.size()) == c_count)
        {
          print_debug(_HERE_, " exact matched DEL with ", edit_calls.to_string(), " index=", e_lo);
          nonsnp_edits.push_back(e_lo);
        }

        ++e_lo;
        get_next_edit_order_at_e_lo(next_edit_order); // advance to the next edit
      }

      sam_record.num_edits += c_count;
      sam_record.alignment_score -= score_penalty;
      std::advance(ref_it, c_count);
      ref_pos += c_count;
      break;
    }

    case paw::CigarOperation::MATCH:
    case paw::CigarOperation::EQUAL:
    case paw::CigarOperation::DIFFERENT:
    {
      for (int i{0}; i < c_count; ++i)
      {
        assert(ref_it != full_ref.cend());
        assert(read_it != full_read.cend());

        update_next_edit(ref_pos, 'X');

        if (*ref_it == *read_it)
        {
          // print_info(_HERE_, " match at i=", i, " in cigar c=", c, " score=", sam_record.alignment_score);
          while (next_edit_order == ::get_edit_order(ref_pos, 'X'))
          {
            assert(ec_ptr);
            assert(e_lo < static_cast<int>(ec_ptr->size()));
            assert(e_lo < e_hi);
            edits.push_back(-e_lo - 1);

#ifndef NDEBUG
            EditCalls const & edit_calls = (*ec_ptr)[e_lo];
            assert(edit_calls.pos == ref_pos);
            assert(edit_calls.type == 'X');
            assert(edit_calls.seq.size() == 1);
            assert(isACGT(edit_calls.seq[0]));

            print_debug(_HERE_, " found a read with NO edit ", edit_calls.to_string());
#endif // NDEBUG

            ++e_lo;
            get_next_edit_order_at_e_lo(next_edit_order); // advance to the next edit
          }

          sam_record.alignment_score += copts.match;
        }
        else if (!isACGT(*ref_it) || !isACGT(*read_it))
        {
          if (*ref_it != *read_it)
            ++sam_record.num_edits; // according to specs of the AS tag

          // Ambigous base gets a match
          sam_record.alignment_score += copts.match;
        }
        else
        {
          // print_info(_HERE_, " mismatch at i=", i, " in cigar c=", c, " score=", sam_record.alignment_score);
          while (next_edit_order == ::get_edit_order(ref_pos, 'X'))
          {
            // is_at_edit = true;
            assert(ec_ptr);
            assert(e_lo < static_cast<int>(ec_ptr->size()));
            assert(e_lo < e_hi);

            EditCalls const & edit_calls = (*ec_ptr)[e_lo];
            assert(edit_calls.pos == ref_pos);
            assert(edit_calls.type == 'X');
            assert(edit_calls.seq.size() == 1);
            assert(isACGT(edit_calls.seq[0]));

            // Check if the read matches the edit
            if (*read_it == edit_calls.seq[0])
            {
              edits.push_back(e_lo);

              // DEBUG
              print_debug(_HERE_, " found a read with edit ", (*ec_ptr)[e_lo].to_string());
              // DEBUG ends
            }

            ++e_lo;
            get_next_edit_order_at_e_lo(next_edit_order); // advance to the next edit
          }

          sam_record.alignment_score -= copts.mismatch;
          ++sam_record.num_edits;
        }

        ++ref_it;
        ++read_it;
        ++ref_pos;
      }

      break;
    }

    case paw::CigarOperation::SOFT_CLIP:
    {
      sam_record.alignment_score -= copts.clip;
      std::advance(read_it, c_count);

      break;
    }

    case paw::CigarOperation::HARD_CLIP:
    {
      break;
    }

    default:
    {
      print_warning(_HERE_, " unexpected cigar operation: ", static_cast<int>(cigar.operation));
      break;
    }
    }
  }

  assert(ref_it == full_ref.cend());
  assert(read_it == full_read.cend());

  if (edits.size() == 0 && nonsnp_edits.size() == 0)
  {
    sam_record.weaver_score = sam_record.alignment_score;
  }
  else
  {
    assert(ec_ptr);
    assert(ec_ptr->size() > 0);

    double const c = 1.05; // controls the weights

    // Gather the haplotype matches for SNP and non-SNP separately
    std::vector<double> snp_hap_matches;
    double const snp_ref_matches = get_hap_matches(snp_hap_matches, *ec_ptr, edits);

    std::vector<double> nonsnp_hap_matches;

#ifndef NDEBUG
    {
      double const nonsnp_ref_matches = get_hap_matches(nonsnp_hap_matches, *ec_ptr, nonsnp_edits);
      assert(nonsnp_ref_matches == 0);
    }
#else
    get_hap_matches(nonsnp_hap_matches, *ec_ptr, nonsnp_edits);
#endif // NDEBUG

    std::vector<double> hap_matches;
    add_hap_matches(hap_matches, snp_hap_matches, nonsnp_hap_matches);

    // because nonsnp_ref_matches == 0 trivially, then ref_matches == snp_ref_matches
    double const w_log_sum = add_log(c * snp_ref_matches, log_sum_base(/*vals=*/hap_matches, /*base=*/c));
    double w_matches = snp_ref_matches * std::exp(c * snp_ref_matches - w_log_sum);
    double w_nonsnp_matches{0.0};
    assert(hap_matches.size() > 0);

    for (int h{0}; h < static_cast<int>(hap_matches.size()); ++h)
    {
      auto const nonsnp_hap_match = (h < static_cast<int>(nonsnp_hap_matches.size())) ? nonsnp_hap_matches[h] : 0;
      auto const snp_hap_match = (h < static_cast<int>(snp_hap_matches.size())) ? snp_hap_matches[h] : 0;
      auto const hap_match = hap_matches[h];

      w_matches += snp_hap_match * std::exp(c * hap_match - w_log_sum);
      w_nonsnp_matches += nonsnp_hap_match * std::exp(c * hap_match - w_log_sum);
    }

    // Number of SNP and nonSNP edits
    double const n_edits{static_cast<double>(edits.size())};
    double const n_nonsnp_edits{static_cast<double>(nonsnp_edits.size())};

    double const w_score = static_cast<double>(copts.match) * w_matches - //
                           static_cast<double>(copts.mismatch) * (n_edits - w_matches) -
                           static_cast<double>(copts.gap_open) * (n_nonsnp_edits - w_nonsnp_matches);

    double const snp_ref_mismatches = n_edits - snp_ref_matches;
    assert(snp_ref_mismatches >= 0.0);

    // alignment score without any of the edits from the graph
    int const as_without_edits = sam_record.alignment_score +                     //
                                 std::round(snp_ref_mismatches * copts.mismatch - //
                                            snp_ref_matches * copts.match +       //
                                            /*nonsnp_ref_mismatches=*/n_nonsnp_edits * copts.gap_open);

    sam_record.weaver_score = as_without_edits + std::round(w_score);

    // // DEBUG
    // print_info(_HERE_,
    //            " alignment_score=",
    //            sam_record.alignment_score,
    //            " weaver_score=",
    //            sam_record.weaver_score,
    //            " log_sum=",
    //            matches_log_sum,
    //            " ref_match,ref_mismatch=",
    //            static_cast<int>(ref_matches),
    //            ", ",
    //            std::round(edits.size() - ref_matches),
    //            " weighted_matches=",
    //            weighted_matches,
    //            " weighted_score=",
    //            weighted_score);
  }
}

void process_sam_record(GFA const & gfa,
                        HaplotypeStats const & hap_stats,
                        SAMRecord & sam_record,
                        SRSeed const & seed,
                        SRAlignment const & alignment,
                        std::string const & seq_fwd,
                        std::string const & seq_rev)
{
  if (seed.is_empty())
  {
    print_warning(_HERE_, " ignoring record with an empty seed.");
    return;
  }

  if (!is_seed_on_same_stable_sequence(gfa, seed))
  {
    print_warning(_HERE_, " ignoring record with a seed that spans multiple stable sequences.");
    return;
  }

  GFALocation begin_ref(seed.begin_ref_value);
  GFALocation end_ref(seed.end_ref_value);
  std::vector<gfa_arc_t const *> begin_arcs;
  std::vector<gfa_arc_t const *> end_arcs;

  gfa_seg_t const & begin_segment = gfa.get_segment(begin_ref.rid);
  gfa_seg_t const & end_segment = gfa.get_segment(end_ref.rid);

  int const begin_sfa_idx = stable_contigs.get_sfa_idx(begin_segment.snid, begin_segment.soff);
  int const end_sfa_idx = stable_contigs.get_sfa_idx(end_segment.snid, end_segment.soff);

  if (begin_sfa_idx != end_sfa_idx)
    return;

  /*
  if (begin_segment.snid != 0 || end_segment.snid != 0)
  {
    print_info(_HERE_,
               " DEBUG begin name=",
               begin_segment.name,
               " rid=",
               begin_ref.rid,
               " snid=",
               begin_segment.snid,
               " soff=",
               begin_segment.soff,
               " rank=",
               begin_segment.rank,
               " sfa_idx=",
               begin_sfa_idx);

    print_info(_HERE_,
               " DEBUG end name=",
               end_segment.name,
               " rid=",
               end_ref.rid,
               " snid=",
               end_segment.snid,
               " soff=",
               end_segment.soff,
               " rank=",
               end_segment.rank,
               " sfa_idx=",
               end_sfa_idx);
  }
  */

  Contig const & begin_contig = stable_contigs.contigs[begin_sfa_idx];
  Contig const & end_contig = stable_contigs.contigs[end_sfa_idx];

  // gfa_sseq_t const & begin_stable_sequence = gfa.get_stable_sequence(begin_segment);
  // gfa_sseq_t const & end_stable_sequence = gfa.get_stable_sequence(end_segment);

  // print_warning(_HERE_, " begin stable rank=", begin_stable_sequence.rank, " and min=", begin_stable_sequence.min);
  // print_warning(_HERE_, " begin_segment pos=", begin_ref.pos, " begin_segment.soff=", begin_segment.soff);
  int const begin_stable_pos = begin_ref.pos + begin_segment.soff - begin_contig.min;
  int const end_stable_pos = end_ref.pos + end_segment.soff - end_contig.min;

  // int const begin_stable_pos = begin_ref.pos + begin_segment.soff - begin_stable_sequence.min;
  // int const end_stable_pos = end_ref.pos + end_segment.soff - end_stable_sequence.min;

  assert(sam_record.seq.size() == seq_fwd.size());
  assert(sam_record.seq.size() == seq_rev.size());
  int const num_read_begin_bases = sketch_value_pos(seed.begin_read_value);
  int const num_read_end_bases = static_cast<int>(seq_fwd.size()) - 1 - sketch_value_pos(seed.end_read_value);

  print_debug(_HERE_,
              " begin/end_rid=",
              begin_ref.rid,
              " ",
              end_ref.rid,
              " begin/end_stable_pos=",
              begin_stable_pos,
              " ",
              end_stable_pos,
              " num_read_begin/end_bases=",
              num_read_begin_bases,
              " ",
              num_read_end_bases,
              " seed=",
              seed.to_string());

  assert(begin_ref.strand != end_ref.strand);
  assert(sam_record.cig.empty());
  assert(seed.is_valid());

  std::string ref_seed = seed.get_ref_sequence();
  std::string ref_before;
  std::string ref_after;

  // Use strand on "end_ref" because begin_ref is already reversed
  if (end_ref.is_strand_forward())
  {
    print_debug(_HERE_, " read is in forward direction");
    assert(begin_segment.soff >= 0);
    sam_record.pos = begin_stable_pos;
    sam_record.snid = begin_segment.snid;
    sam_record.sfa_idx = begin_sfa_idx;

    if (alignment.begin_alignment_ext != nullptr)
    {
      paw::AlignmentResults const & ar_begin = *alignment.begin_alignment_ext;
      print_debug(_HERE_, " ar_begin:", alignment_results_to_string(ar_begin));
      process_cigar_forward_begin(sam_record, ar_begin, num_read_begin_bases);
    }
    else
    {
      print_debug(_HERE_, " no begin alignment, softclipping ", num_read_begin_bases);
      sam_record.append_cigar(num_read_begin_bases, paw::CigarOperation::SOFT_CLIP);
    }

    // advance begin portion
    {
      int const advanced_ref_begin = sam_record.get_cigar_reference_length();
      std::string ref_before_rev;
      advance_from_seed(advanced_ref_begin, begin_ref, begin_arcs, ref_before_rev);
      assert(static_cast<int>(ref_before_rev.size()) == advanced_ref_begin);
      ref_before = get_reverse_complement(ref_before_rev);
    }

    assert(sam_record.get_cigar_query_length() == num_read_begin_bases);

    if (seed.cig == nullptr)
    {
      sam_record.append_cigar(seed.get_read_length(), paw::CigarOperation::MATCH); // Add cigar from seed
    }
    else
    {
#ifndef NDEBUG
      int query_count{0};
#endif // NDEBUG

      for (auto cig_it = seed.cig->rbegin(); cig_it != seed.cig->rend(); ++cig_it)
      {
        int const cigar_count = cig_it->count;
        paw::CigarOperation inv_op = inv_cigar(cig_it->operation);

#ifndef NDEBUG
        if (paw::advances_query(inv_op))
          query_count += cigar_count;
#endif // NDEBUG

        sam_record.append_cigar(cigar_count, inv_op);
      }

      assert(query_count == seed.get_read_length());
    }

    assert(sam_record.get_cigar_query_length() == num_read_begin_bases + seed.get_read_length());

    if (alignment.end_alignment_ext != nullptr)
    {
      paw::AlignmentResults const & ar_end = *alignment.end_alignment_ext;
      print_debug(_HERE_, " printing ar_end:", alignment_results_to_string(ar_end));
      process_cigar_forward_end(sam_record, ar_end, num_read_end_bases);
    }
    else
    {
      print_debug(_HERE_, " no end alignment, softclipping ", num_read_end_bases);
      sam_record.append_cigar(num_read_end_bases, paw::CigarOperation::SOFT_CLIP);
    }

    assert(sam_record.get_cigar_query_length() == num_read_begin_bases + seed.get_read_length() + num_read_end_bases);

    // advance end portion
    {
      int const advanced_ref_end = sam_record.get_cigar_reference_length() - ref_before.size() - seed.length;
      advance_from_seed(advanced_ref_end, end_ref, end_arcs, ref_after);
      assert(static_cast<int>(ref_after.size()) == advanced_ref_end);
    }
  }
  else
  {
    print_debug(_HERE_, " read is in reverse direction");
    assert(end_segment.soff >= 0);
    sam_record.pos = end_stable_pos;
    sam_record.snid = end_segment.snid;
    sam_record.sfa_idx = end_sfa_idx;

    if (alignment.end_alignment_ext)
    {
      paw::AlignmentResults const & ar_end = *alignment.end_alignment_ext;
      print_debug(_HERE_, " printing ar_end:", alignment_results_to_string(ar_end));
      process_cigar_reverse_end(sam_record, ar_end, num_read_end_bases);
    }
    else
    {
      print_debug(_HERE_, " no end alignment, softclipping ", num_read_end_bases);
      sam_record.append_cigar(num_read_end_bases, paw::CigarOperation::SOFT_CLIP);
    }

    // advance end part
    {
      int const advanced_ref_end = sam_record.get_cigar_reference_length();
      std::string ref_before_rev;
      advance_from_seed(advanced_ref_end, end_ref, end_arcs, ref_before_rev);
      assert(static_cast<int>(ref_before_rev.size()) == advanced_ref_end);
      ref_before = get_reverse_complement(ref_before_rev);
    }

    assert(sam_record.get_cigar_query_length() == num_read_end_bases);

    if (seed.cig == nullptr)
    {
      // This means there is perfect match between the query and the reference
      sam_record.append_cigar(seed.length, paw::CigarOperation::MATCH); // Add cigar from seed
    }
    else
    {
#ifndef NDEBUG
      int query_count{0};
#endif // NDEBUG

      for (auto cig_it = seed.cig->begin(); cig_it != seed.cig->end(); ++cig_it)
      {
        int const cigar_count = cig_it->count;
        paw::CigarOperation inv_op = inv_cigar(cig_it->operation);

#ifndef NDEBUG
        if (paw::advances_query(inv_op))
          query_count += cigar_count;
#endif // NDEBUG

        sam_record.append_cigar(cigar_count, inv_op);
      }

      assert(query_count == seed.get_read_length());
    }

    assert(sam_record.get_cigar_query_length() == num_read_end_bases + seed.get_read_length());

    if (alignment.begin_alignment_ext)
    {
      paw::AlignmentResults const & ar_begin = *alignment.begin_alignment_ext;
      print_debug(_HERE_, " printing ar_begin:", alignment_results_to_string(ar_begin));
      process_cigar_reverse_begin(sam_record, ar_begin, num_read_begin_bases);
    }
    else
    {
      print_debug(_HERE_, " no begin alignment, softclipping ", num_read_begin_bases);
      sam_record.append_cigar(num_read_begin_bases, paw::CigarOperation::SOFT_CLIP);
    }

    assert(sam_record.get_cigar_query_length() == num_read_end_bases + seed.get_read_length() + num_read_begin_bases);

    // advance begin part
    {
      int const advanced_ref_begin = sam_record.get_cigar_reference_length() - ref_before.size() - seed.length;
      assert(ref_after.empty()); // ref_after should not be set
      advance_from_seed(advanced_ref_begin, begin_ref, begin_arcs, ref_after);
      assert(static_cast<int>(ref_after.size()) == advanced_ref_begin);
    }
  }

  // Check if the begin and end positions are outside of the stable fasta contig sequence
  {
    assert(sam_record.sfa_idx >= 0);
    assert(sam_record.sfa_idx < static_cast<int>(stable_contigs.contigs.size()));
    int const contig_length = stable_contigs.contigs[sam_record.sfa_idx].get_length();
    int const ref_reach = sam_record.get_reference_reach();
    assert(contig_length > 0);

    // end position
    if (ref_reach > contig_length)
    {
      int const clip_amount = std::min(static_cast<int>(ref_after.size()), ref_reach - contig_length);

      sam_record.soft_clip_cigar_end_based_on_ref(clip_amount);
      // ref_after.erase(ref_after.size() - clip_amount, clip_amount);
      ref_after.resize(ref_after.size() - clip_amount);
    }

    // begin position
    if (sam_record.pos < 0)
    {
      int const clip_amount = std::min(static_cast<int>(ref_before.size()), -sam_record.pos);
      assert(clip_amount > 0);

      sam_record.soft_clip_cigar_begin_based_on_ref(clip_amount);
      assert(sam_record.pos == 0); // pos is adjusted in SAMRecord::soft_clip_cigar_begin_based_on_ref()
      ref_before.erase(0, clip_amount);
    }
  }

  assert(sam_record.pos != SAMRecord::MISSING_POS);
  assert(sam_record.is_cigar_valid());

  // Check alignment score and if indels are left aligned
  calculate_alignment_score(hap_stats,
                            sam_record,
                            ref_seed,
                            std::move(ref_before),
                            ref_after,
                            end_ref.is_strand_forward(),
                            seq_fwd,
                            seq_rev);
}

void adapter_removal(SAMRecord & main_record, SAMRecord const & other_record, int const max_removal)
{
  if (max_removal <= 0)
    return;

  bool const seq1_reversed = (main_record.flags & SAMFlags::IS_SEQ_REVERSED) != 0u;
  bool const seq2_reversed = (main_record.flags & SAMFlags::IS_MATE_SEQ_REVERSED) != 0u;

  // check for sequence adapters
  if ((not seq1_reversed) && seq2_reversed)
  {
    // forward => check if end expands over the other record
    int const this_ref_reach = main_record.get_reference_reach();
    int const other_ref_reach = other_record.get_reference_reach();
    int const clip_amount = this_ref_reach - other_ref_reach;

    if (clip_amount > 0 && clip_amount <= max_removal)
      main_record.soft_clip_cigar_end_based_on_ref(clip_amount);
  }
  else if (seq1_reversed && (not seq2_reversed))
  {
    // reversed
    int const clip_amount = other_record.pos - main_record.pos; // negative indicates nothing should clipped

    if (clip_amount > 0 && clip_amount <= max_removal)
      main_record.soft_clip_cigar_begin_based_on_ref(clip_amount);
  }
}

void process_sam_record_pair(GFA const & gfa, //
                             T_icu const & icu,
                             SAMRecord & sam_record1,
                             SAMRecord & sam_record2,
                             SRSeed const & seed1,
                             SRSeed const & seed2)
{
  // check adapters
  Options const & copts = *(Options::const_instance());

  // Check if record expands over the other record
  if (not copts.is_no_adapter_removal && sam_record1.sfa_idx == sam_record2.sfa_idx)
  {
    int const max_removal =
      std::min(get_cigar_reference_length(sam_record1.cig), get_cigar_reference_length(sam_record2.cig)) - 10;
    adapter_removal(sam_record1, sam_record2, max_removal);
    adapter_removal(sam_record2, sam_record1, max_removal);
  }

  // stop here if either read has been marked as unmapped
  if ((sam_record1.flags & SAMFlags::IS_UNMAPPED) != 0 || (sam_record2.flags & SAMFlags::IS_UNMAPPED) != 0)
    return;

  int const sdist = estimate_shortest_distance(gfa, icu, seed1.begin_ref_value ^ 1, seed2.begin_ref_value);

#ifndef NDEBUG
  assert(sdist == estimate_shortest_distance(gfa, icu, seed2.begin_ref_value ^ 1, seed1.begin_ref_value));

  print_debug(_HERE_,
              " sdist=",
              sdist,
              " reversed=",
              estimate_shortest_distance(gfa, icu, seed2.begin_ref_value ^ 1, seed1.begin_ref_value));
#endif // NDEBUG

  if (sdist < 0)
  {
    print_debug(_HERE_,
                " ",
                sketch_value_to_string(seed1.begin_ref_value ^ 1),
                " -> ",
                sketch_value_to_string(seed2.begin_ref_value));
    return;
  }

  if (sdist <= 1500)
  {
    sam_record1.flags |= SAMFlags::IS_PROPER_PAIR;
    sam_record2.flags |= SAMFlags::IS_PROPER_PAIR;
  }

  if (sam_record1.sfa_idx == sam_record2.sfa_idx)
  {
    if (sam_record1.pos <= sam_record2.pos)
    {
      int dist = sam_record2.pos + sam_record2.get_cigar_reference_length() - sam_record1.pos;
      sam_record1.tlen = dist;
      sam_record2.tlen = -dist;
    }
    else
    {
      int dist = sam_record1.pos + sam_record1.get_cigar_reference_length() - sam_record2.pos;
      sam_record1.tlen = -dist;
      sam_record2.tlen = dist;
    }
  }
  else
  {
    // TODO get a more accurate tlen for different sfa_idx
    if (sam_record1.sfa_idx < sam_record2.sfa_idx)
    {
      sam_record1.tlen = sdist;
      sam_record2.tlen = -sdist;
    }
    else
    {
      sam_record1.tlen = -sdist;
      sam_record2.tlen = sdist;
    }
  }
}

void prepare_sam_write(GFA const & gfa,
                       T_icu const & icu,
                       SRSeed const & seed1,
                       SAMRecord & record1,
                       SRSeed const & seed2,
                       SAMRecord & record2)
{
  // set flags
  record1.flags = SAMFlags::IS_PAIRED | SAMFlags::IS_FIRST_IN_PAIR;
  record2.flags = SAMFlags::IS_PAIRED | SAMFlags::IS_SECOND_IN_PAIR;

  bool const read1_unmapped = seed1.get_read_length() == 0;
  bool const read2_unmapped = seed2.get_read_length() == 0;

  if (read1_unmapped)
  {
    print_debug(_HERE_, " No seed for read1=", record1.qname);
    record1.make_unmapped(record2);
  }

  if (read2_unmapped)
  {
    print_debug(_HERE_, " No seed for read2=", record2.qname);
    record2.make_unmapped(record1);
  }

  // Make sure the seed does not jump between stable sequences.
  // The two seeds can still be on different stable sequence.
  if (is_seed_on_same_stable_sequence(gfa, seed1) && //
      is_seed_on_same_stable_sequence(gfa, seed2))
  {
    // set is seq reverse flags
    process_sam_record_seq_reverse_flag(record1 /*record of seed*/, record2 /*other record*/, seed1);
    process_sam_record_seq_reverse_flag(record2 /*record of seed*/, record1 /*other record*/, seed2);

    // set paired read flags
    process_sam_record_pair(gfa, icu, record1, record2, seed1, seed2);
  }
  else
  {
    print_warning(_HERE_, " Marked ", record1.qname, " as unpaired as it was on different stable contigs.");
    record1.make_unmapped(record2);
    record2.make_unmapped(record1);
  }

#ifndef NDEBUG
  if constexpr (true) // set to true to check output reads
  {
    // These should never happen but they do, TODO investigate
    if (!record1.is_pair_valid(record2))
    {
      print_debug("Making a record unampped due to a sanity check.");
      record1.make_unmapped(record2);
    }

    if (!record2.is_pair_valid(record1))
    {
      print_debug("Making a record unampped due to a sanity check.");
      record2.make_unmapped(record1);
    }

    // check the sam records
    assert(record1.is_pair_valid(record2));
    assert(record2.is_pair_valid(record1));
  }
#endif // NDEBUG
}

void sam_write(SAMWriter & sam_writer, //
               SAMRecord & record1,
               SAMRecord & record2)
{
  // get cigar strings for both reads
  std::string cigar1 = record1.cig.empty() ? std::string(1, '*') : cigar2string(record1.cig.begin(), record1.cig.end());
  std::string cigar2 = record2.cig.empty() ? std::string(1, '*') : cigar2string(record2.cig.begin(), record2.cig.end());

  // record1.mpos = record2.pos;
  // record1.mmapq = record2.mapq;
  //
  // record2.mpos = record1.pos;
  // record2.mmapq = record1.mapq;

  sam_writer.write_sr_record(record1, record2, cigar1, cigar2);
  sam_writer.write_sr_record(record2, record1, cigar2, cigar1);
}

void append_sam_lines(std::vector<std::vector<std::pair<uint64_t, std::string>>> & bucket_of_sam_lines, //
                      SAMRecord & record1,
                      SAMRecord & record2)
{
  // Get cigar strings for both reads
  std::string cigar1 = record1.cig.empty() ? std::string(1, '*') : cigar2string(record1.cig.begin(), record1.cig.end());
  std::string cigar2 = record2.cig.empty() ? std::string(1, '*') : cigar2string(record2.cig.begin(), record2.cig.end());

  // record1.mpos = record2.pos;
  // record1.mmapq = record2.mapq;
  //
  // record2.mpos = record1.pos;
  // record2.mmapq = record1.mapq;

  // Read 1
  {
    uint64_t const ppos1 = record1.get_sam_order();
    int const bucket_index1 = stable_contigs.get_bucket_index(ppos1);
    assert(bucket_index1 < static_cast<int>(bucket_of_sam_lines.size()));

    bucket_of_sam_lines[bucket_index1].emplace_back(
      ppos1,
      get_sam_string(record1, record2, /*cigar=*/cigar1, /*other_cigar=*/cigar2));
  }

  // Read 2
  {
    uint64_t const ppos2 = record2.get_sam_order();
    int const bucket_index2 = stable_contigs.get_bucket_index(ppos2);
    assert(bucket_index2 < static_cast<int>(bucket_of_sam_lines.size()));

    bucket_of_sam_lines[bucket_index2].emplace_back(
      ppos2,
      get_sam_string(record2, record1, /*cigar=*/cigar2, /*other_cigar=*/cigar1));
  }
}

} // namespace weaver
