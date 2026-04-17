#include "sr_alignment.hpp"

#include <cassert>
#include <cstdint>
#include <sstream>

#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/pairwise_alignment.hpp>
#include <paw/align/pairwise_ext_alignment.hpp>

#include "gfa.hpp"
#include "gfa_location.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "sam_record.hpp"
#include "segment.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp"
#include "sr_cigar.hpp" // cigar2string
#include "sr_sam.hpp"
#include "sr_seed.hpp"

namespace
{
using Tuint = uint16_t;
int constexpr PAD_UNTIL_MISMATCH{12};

uint64_t extend_no_indel(weaver::GFALocation gfa_location,
                         std::vector<gfa_arc_t const *> & arcs,
                         std::string_view query,
                         int & mismatches,
                         int const threshold)
{
  using namespace weaver;

  if (query.size() <= 0ull)
    return 0ull;

  assert(gfa_location.is_valid());
  int unaccounted = gfa_location.advance_when_same_contig(1, arcs);

  if (unaccounted > 0)
    return 0ull;

  assert(gfa_location.is_valid());
  std::string seq;
  unaccounted = gfa_location.advance_when_same_contig_and_get_sequence(query.size() - 1, arcs, seq);

  assert(gfa_location.is_valid());
  gfa_location.get_base(seq);
  assert(unaccounted + seq.size() == query.size());
  print_debug(_HERE_, " query begin/end=", query, " graph begin/end=", seq);

  for (uint64_t b{0}; b < seq.size(); ++b)
  {
    assert(b < query.size());

    if (seq[b] == 'X' || query[b] == 'X' || seq[b] != query[b])
    {
      ++mismatches;

      if (mismatches > threshold)
        return b;
    }
  }

  return seq.size();
}

template <typename Tuint>
void alignment(paw::AlignmentOptions<Tuint> & aln_opts,
               std::string_view after,
               std::string_view ref,
               std::unique_ptr<paw::AlignmentResults> & results_ptr)
{
  weaver::Options const & copts = *(weaver::Options::const_instance());

  aln_opts.set_match(copts.match)
    .set_mismatch(copts.mismatch)
    .set_gap_open(copts.gap_open)
    .set_gap_extend(copts.gap_extend)
    .set_clip(copts.clip);

#ifndef NDEBUG
  aln_opts.get_aligned_strings = true;
#endif
  aln_opts.get_cigar_string = true;
  aln_opts.is_clip = true;
  aln_opts.left_column_free = false;
  aln_opts.right_column_free = true;

  print_debug(_HERE_, " alignment query=", after, " ref=", ref);

  paw::pairwise_ext_alignment(after, ref, aln_opts);
  print_debug(_HERE_, " alignment resuts: ", weaver::alignment_results_to_string(*aln_opts.ar));

  assert(aln_opts.ar);
  assert(aln_opts.ar->query_begin == 0);
  assert(aln_opts.ar->query_end <= static_cast<int>(after.size()));

  int const cigar_len = static_cast<int>(after.size()) - aln_opts.ar->query_end;

  if (cigar_len > 0)
  {
    paw::Cigar new_item(cigar_len, paw::CigarOperation::SOFT_CLIP);
    aln_opts.ar->cigar_string_ptr->insert(aln_opts.ar->cigar_string_ptr->begin(), new_item);
  }

  results_ptr = std::move(aln_opts.ar); // change ownership
  assert(aln_opts.ar == nullptr);
  assert(results_ptr != nullptr);

#ifndef NDEBUG
  assert(results_ptr->aligned_strings_ptr);
  std::pair<std::string, std::string> const & aln_strings = *(results_ptr->aligned_strings_ptr);

  print_debug(_HERE_, " Aligned strings: \n", aln_strings.first, '\n', aln_strings.second);
#endif // NDEBUG
}

//! Extend and align the "end" portion of the read
void extend_and_align_end(std::unique_ptr<paw::AlignmentResults> & end_alignment,
                          std::vector<gfa_arc_t const *> & end_arcs,
                          weaver::SRSeed & seed,
                          std::string_view query_seq)
{
  using namespace weaver;

  int const read_end_pos = sketch_value_pos(seed.end_read_value);

  // +1 to change to 1-based indexing
  int const remaining_end_query_size = static_cast<int>(query_seq.size()) - (read_end_pos + 1);
  int const threshold{1};
  assert(remaining_end_query_size > 0);

  // extend from end, modifying both read_end_pos and ref_end_pos
  std::string_view after(query_seq);

  // query_seq contains the full query sequence, here we determine the "end" part of the query sequence that
  // has not been mapped to the graph. This should in normal circumstances be more than 0
  after.remove_prefix(read_end_pos + 1);
  assert(static_cast<int>(after.size()) == remaining_end_query_size);

  int mismatches{0};

  // first try to extend exactly by the remaining_end_query_size
  extend_no_indel(GFALocation(seed.end_ref_value), end_arcs, after, mismatches, threshold);

  if (mismatches <= threshold)
  {
    print_debug(_HERE_, " mismatches=", mismatches);
    int const matches = remaining_end_query_size - mismatches;
    assert(matches >= 0);
    end_alignment = std::make_unique<paw::AlignmentResults>();
    weaver::Options const & copts = *(weaver::Options::const_instance());
    end_alignment->score = matches * copts.match - copts.mismatch * mismatches;
  }
  else
  {
    // clear the outdates values
    end_arcs.clear();
    mismatches = 0;

    GFALocation end_ref(seed.end_ref_value);

    // Optimize seeding
    // if (false)
    {
      // extend until a mismatch is seen
      std::vector<gfa_arc_t const *> dummy_arcs;
      uint64_t by = extend_no_indel(end_ref, dummy_arcs, after, mismatches, 0);

      if (by > PAD_UNTIL_MISMATCH && dummy_arcs.empty())
      {
        by -= PAD_UNTIL_MISMATCH;
        assert(by > 0);

#ifndef NDEBUG
        int by2 = end_ref.advance_when_same_contig(by, end_arcs);
        assert(by2 == 0);
#else
        end_ref.advance_when_same_contig(by, end_arcs);
#endif

        seed.end_ref_value = end_ref.get_value();
        seed.end_read_value += (by << 1);
        seed.length += by;
        after.remove_prefix(by);

        if (seed.cig)
          weaver::append_cigar_back(*seed.cig, static_cast<uint32_t>(by), paw::CigarOperation::MATCH);
      }

      assert(end_ref.is_same_location(GFALocation(seed.end_ref_value)));
    }

    // do a more expensive alignment
    end_ref.advance_when_same_contig(1, end_arcs);
    std::string ref;

    {
      int constexpr MAX_EXTENSION{70}; // TODO make an option
      assert(end_ref.is_valid());

      int unaccounted = end_ref.advance_when_same_contig_and_get_sequence( //
        after.size() + MAX_EXTENSION,                                      //
        end_arcs,                                                          //
        ref);

      assert(end_ref.is_valid());
      assert(unaccounted >= 0);

      if (unaccounted > 0)
      {
        seed.end_unaccounted_read_bases = unaccounted;
        ref += std::string(unaccounted, 'X');
      }
    }

    // We need to do a more expensive alignment
    paw::AlignmentOptions<uint8_t> opts_uint8;
    paw::AlignmentOptions<uint16_t> opts_uint16;
    weaver::Options const & copts = *(weaver::Options::const_instance());
    assert(copts.mismatch > 0);

    if (static_cast<int>(after.size()) < (200 / std::max(copts.match + copts.mismatch - copts.gap_extend, 2)))
      alignment(opts_uint8, after, ref, end_alignment);
    else
      alignment(opts_uint16, after, ref, end_alignment);
  }
}

//! Extend and align the "begin" portion of the read
void extend_and_align_begin(std::unique_ptr<paw::AlignmentResults> & begin_alignment,
                            std::vector<gfa_arc_t const *> & begin_arcs,
                            weaver::SRSeed & seed,
                            std::string_view query_seq_rev)
{
  using namespace weaver;

  int const read_begin_pos = sketch_value_pos(seed.begin_read_value);
  int const remaining_begin_query_size = read_begin_pos;
  int const threshold{1};
  assert(remaining_begin_query_size > 0);

  // extend from begin, modifying both read_begin_pos and ref_begin_pos
  std::string_view before(query_seq_rev);

  // similar as above, but for the "begin" part of the query
  before.remove_prefix(query_seq_rev.size() - read_begin_pos);
  assert(static_cast<int>(before.size()) == remaining_begin_query_size);
  int mismatches{0};

  {
    GFALocation begin_ref_location(seed.begin_ref_value);
    extend_no_indel(begin_ref_location, begin_arcs, before, mismatches, threshold);
  }

  if (mismatches <= threshold)
  {
    print_debug(_HERE_, " mismatches=", mismatches);
    int const matches = remaining_begin_query_size - mismatches;
    assert(matches >= 0);
    begin_alignment = std::make_unique<paw::AlignmentResults>();
    weaver::Options const & copts = *(weaver::Options::const_instance());
    begin_alignment->score = matches * copts.match - copts.mismatch * mismatches;
  }
  else
  {
    // clear the arcs
    begin_arcs.clear();
    GFALocation begin_ref(seed.begin_ref_value);

    // Optimize seeding
    // if (false)
    {
      // extend until a mismatch is seen
      std::vector<gfa_arc_t const *> dummy_arcs;
      uint64_t by = extend_no_indel(begin_ref, dummy_arcs, before, mismatches, 0);

      if (by > PAD_UNTIL_MISMATCH && dummy_arcs.empty())
      {
        by -= PAD_UNTIL_MISMATCH;
// upate end_ref_value
#ifndef NDEBUG
        int by2 = begin_ref.advance_when_same_contig(by, begin_arcs);
        assert(by2 == 0);
#else
        begin_ref.advance_when_same_contig(by, begin_arcs);
#endif
        seed.begin_ref_value = begin_ref.get_value();
        seed.begin_read_value -= (by << 1);
        seed.length += by;
        before.remove_prefix(by);

        if (seed.cig)
          weaver::append_cigar_front(*seed.cig, static_cast<uint32_t>(by), paw::CigarOperation::MATCH);
      }

      assert(begin_ref.is_same_location(GFALocation(seed.begin_ref_value)));
    }

    begin_ref.advance_when_same_contig(1, begin_arcs);
    std::string ref;

    {
      int constexpr MAX_EXTENSION{70}; // TODO make an option
      assert(begin_ref.is_valid());

      int unaccounted = begin_ref.advance_when_same_contig_and_get_sequence(before.size() + MAX_EXTENSION, //
                                                                            begin_arcs,
                                                                            ref);

      assert(begin_ref.is_valid());
      assert(unaccounted >= 0);

      if (unaccounted > 0)
      {
        seed.begin_unaccounted_read_bases = unaccounted;
        ref += std::string(unaccounted, 'X');
      }
    }

    // We need to do a more expensive alignment
    paw::AlignmentOptions<uint8_t> opts_uint8;
    paw::AlignmentOptions<uint16_t> opts_uint16;

    weaver::Options const & copts = *(weaver::Options::const_instance());
    assert(copts.mismatch > 0);

    if (static_cast<int>(before.size()) < (200 / std::max(copts.match + copts.mismatch - copts.gap_extend, 2)))
      alignment(opts_uint8, before, ref, begin_alignment);
    else
      alignment(opts_uint16, before, ref, begin_alignment);
  }
}

} // namespace

namespace weaver
{
void SRAlignment::clear()
{
  clear_begin_alignment_extension();
  clear_end_alignment_extension();
  score = SRAlignment::MISSING_SCORE;
}

void SRAlignment::clear_begin_alignment_extension()
{
  begin_alignment_ext = nullptr;
  begin_arcs.clear();
}

void SRAlignment::clear_end_alignment_extension()
{
  end_alignment_ext = nullptr;
  end_arcs.clear();
}

int SRAlignment::get_ext_score() const
{
  if (begin_alignment_ext == nullptr || end_alignment_ext == nullptr)
    return SRAlignment::MISSING_SCORE;
  else
    return begin_alignment_ext->score + end_alignment_ext->score;
}

std::string alignment_results_to_string(paw::AlignmentResults const & ar)
{
  std::ostringstream ss;

  ss << "Alignment: score=" << ar.score                                       //
     << " query_begin/end=" << ar.query_begin << " " << ar.query_end          //
     << " database_begin/end=" << ar.database_begin << " " << ar.database_end //
     << " clip_begin/end=" << ar.clip_begin << " " << ar.clip_end;            //

  if (ar.cigar_string_ptr != nullptr)
  {
    ss << " cigar=" << cigar2string(ar.cigar_string_ptr->begin(), ar.cigar_string_ptr->end());
  }

  return ss.str();
}

SRAlignment extend_and_align(SRSeed & seed, std::string_view query_seq, std::string_view query_seq_rev)
{
  assert(seed.is_valid()); // check if seed is ok

  seed.align_to_graph(query_seq);

  // check if seed is ok again, now extra checks will be made because the cigar string is available
  assert(seed.is_valid());

  SRAlignment alignment;

  // check if it is a bad seed that had no mem
  if (seed.score == SRAlignment::MISSING_SCORE)
  {
    alignment.score = SRAlignment::MISSING_SCORE;
    return alignment;
  }

  extend_and_align_end(alignment.end_alignment_ext, alignment.end_arcs, seed, query_seq);
  assert(seed.is_valid()); // check if seed is ok
  assert(GFALocation::gfa->are_arcs_with_same_rank_as_segments(alignment.end_arcs));

  extend_and_align_begin(alignment.begin_alignment_ext, alignment.begin_arcs, seed, query_seq_rev);
  assert(seed.is_valid()); // check if seed is ok
  assert(GFALocation::gfa->are_arcs_with_same_rank_as_segments(alignment.begin_arcs));

#ifndef NDEBUG
  if (alignment.end_arcs.size() > 0 || alignment.begin_arcs.size() > 0)
  {
    print_debug(_HERE_, " n_begin_arcs=", alignment.begin_arcs.size(), " n_end_arcs=", alignment.begin_arcs.size());
  }
#endif // NDEBUG

  // check final results are in place
  assert(alignment.begin_alignment_ext);
  assert(alignment.end_alignment_ext);
  assert(alignment.begin_alignment_ext->score != SRAlignment::MISSING_SCORE);
  assert(alignment.end_alignment_ext->score != SRAlignment::MISSING_SCORE);

  alignment.score = seed.get_score() + alignment.get_ext_score();

#ifndef NDEBUG
  print_debug(_HERE_, " alignment.score=", alignment.score);
#endif // NDEBUG

  return alignment;
}

bool seed_alignment_order_gt(SRSeed const & s, SRAlignment const & sa, SRSeed const & o, SRAlignment const & oa)
{
  int const s_score = get_est_seed_alignment_score(s, sa);
  int const o_score = get_est_seed_alignment_score(o, oa);

  return s_score > o_score || (s_score == o_score && s.length > o.length);
}

int get_est_seed_alignment_score(SRSeed const & seed, SRAlignment const & sa)
{
  assert(sa.begin_alignment_ext);
  assert(sa.end_alignment_ext);
  return seed.get_est_score() + sa.get_ext_score();
}

void remove_seeds_with_no_score(std::vector<SRSeed> & seeds, std::vector<SRAlignment> & aln)
{
  int const n = static_cast<int>(seeds.size());
  int new_s{0};

  for (int s{0}; s < n; ++s)
  {
    if (seeds[s].score > SRAlignment::MISSING_SCORE)
    {
      // this seed has a score
      assert(s >= new_s);

      if (s > new_s)
      {
        seeds[new_s] = std::move(seeds[s]);
        aln[new_s] = std::move(aln[s]);
      }

      ++new_s;
    }
  }

  if (new_s < n)
  {
    seeds.resize(new_s);
    aln.resize(new_s);
  }
}

void remove_duplicates(GFA const & gfa, T_icu const & icu, std::vector<SRSeed> & seeds, std::vector<SRAlignment> & aln)
{
  if (seeds.size() <= 1)
    return;

  int const n = seeds.size();
  int new_s1{0};

  for (int s1{0}; s1 < n; ++s1)
  {
    SRSeed const & seed1 = seeds[s1];
    SRAlignment const & aln1 = aln[s1];
    auto const s1_score = seed1.get_score() + aln1.get_ext_score();
    int const pos1 = sketch_value_pos(seed1.begin_read_value);
    bool is_good{true};

    for (int s2{0}; s2 < n; ++s2)
    {
      if (s1 == s2)
        continue;

      SRSeed const & seed2 = seeds[s2];
      SRAlignment const & aln2 = aln[s2];

      if (aln2.begin_alignment_ext == nullptr)
        continue; // this has been deleted

      auto const s2_score = seed2.get_score() + aln2.get_ext_score();
      int const pos2 = sketch_value_pos(seed2.begin_read_value);

      if ((s1_score < s2_score || (s1 < s2 && s1_score == s2_score)))
      {
        if (pos1 >= pos2)
        {
          if (is_within_distance(gfa, icu, seed1.begin_ref_value, seed2.begin_ref_value, pos1 - pos2 + 100))
          {
            // print_info(_HERE_, " bad seed=", s1);
            is_good = false;
            break;
          }
        }
        else
        {
          if (is_within_distance(gfa, icu, seed2.begin_ref_value, seed1.begin_ref_value, pos2 - pos1 + 100))
          {
            // print_info(_HERE_, " bad seed=", s1);
            is_good = false;
            break;
          }
        }
      }
    } // for s2 ends

    if (is_good)
    {
      assert(s1 >= new_s1);

      if (s1 > new_s1)
      {
        seeds[new_s1] = std::move(seeds[s1]);
        aln[new_s1] = std::move(aln[s1]);
      }

      ++new_s1;
    }
  } // for s1 ends

  if (new_s1 < n)
  {
    seeds.resize(new_s1);
    aln.resize(new_s1);
  }
}

void remove_duplicates(GFA const & gfa,
                       T_icu const & icu,
                       std::vector<SRSeed> & seeds,
                       std::vector<SRAlignment> & aln,
                       std::vector<SAMRecord> & recs)
{
  if (seeds.size() <= 1)
    return;

  int const n = seeds.size();
  int new_s1{0};

  for (int s1{0}; s1 < n; ++s1)
  {
    SRSeed const & seed1 = seeds[s1];
    SRAlignment const & aln1 = aln[s1];
    auto const s1_score = seed1.get_score() + aln1.get_ext_score();
    int const pos1 = sketch_value_pos(seed1.begin_read_value);
    bool is_good{true};

    for (int s2{0}; s2 < n; ++s2)
    {
      if (s1 == s2)
        continue;

      SRSeed const & seed2 = seeds[s2];
      SRAlignment const & aln2 = aln[s2];

      if (aln2.begin_alignment_ext == nullptr)
        continue; // this has been deleted

      auto const s2_score = seed2.get_score() + aln2.get_ext_score();
      int const pos2 = sketch_value_pos(seed2.begin_read_value);

      if ((s1_score < s2_score || (s1 < s2 && s1_score == s2_score)))
      {
        if (pos1 >= pos2)
        {
          if (is_within_distance(gfa, icu, seed1.begin_ref_value, seed2.begin_ref_value, pos1 - pos2 + 100))
          {
            // print_info(_HERE_, " bad seed=", s1);
            is_good = false;
            break;
          }
        }
        else
        {
          if (is_within_distance(gfa, icu, seed2.begin_ref_value, seed1.begin_ref_value, pos2 - pos1 + 100))
          {
            // print_info(_HERE_, " bad seed=", s1);
            is_good = false;
            break;
          }
        }
      }
    } // for s2 ends

    if (is_good)
    {
      assert(s1 >= new_s1);

      if (s1 > new_s1)
      {
        seeds[new_s1] = std::move(seeds[s1]);
        aln[new_s1] = std::move(aln[s1]);
        recs[new_s1] = std::move(recs[s1]);
      }

      ++new_s1;
    }
  } // for s1 ends

  if (new_s1 < n)
  {
    seeds.resize(new_s1);
    aln.resize(new_s1);
    recs.resize(new_s1);
  }
}

} // namespace weaver
