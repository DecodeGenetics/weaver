#include "sr_alignment_sel.hpp"

#include <algorithm> // std::min
#include <cmath>     // pow
#include <cstdint>   // uint64_t
#include <numeric>
#include <vector> // std::vector

#include "alignment_utils.hpp"
#include "gfa.hpp" // GFA
#include "icu.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "sam_record.hpp" // SAMRecord
#include "sr_seed.hpp"

namespace
{
int constexpr NO_PAIR_PENALTY{60};
int constexpr PAIR_BONUS{0};

//! Calculate the read alignment identity in range [0,1]
double get_alignment_identity(int const alignment_score, int const read_length)
{
  using namespace weaver;

  assert(read_length > 0);

  double const s = static_cast<double>(alignment_score);
  double const l = static_cast<double>(read_length);

  Options const & copts = *(Options::const_instance());

  print_debug(_HERE_, " length=", l, " score=", s);

  double const difference = std::max(0.0,
                                     (l * static_cast<double>(copts.match) - s) / //
                                       l /                                        //
                                       static_cast<double>(copts.match + copts.mismatch));

  assert(difference >= 0.0);
  assert(difference <= 1.0);
  return 1.0 - difference;
}

//! Out of equal possible indexes, select one "randomly" but consistently.
int select_an_index(std::vector<int> const & p, int score)
{
  assert(not p.empty());

  if (p.size() == 1)
    return p[0];

  // select pseudo randomly (but consistently) from the indexes
  uint64_t const sum = std::accumulate(p.begin(), p.end(), static_cast<uint64_t>(score));
  return p[sum % p.size()];
}

weaver::SRAlignmentSel select_unpaired_alignment(std::vector<weaver::SAMRecord> const & records, bool const is_read1)
{
  using namespace weaver;

  Options const & copts = *(Options::const_instance());
  SRAlignmentSel sas;
  int const n_scores{static_cast<int>(records.size())};
  assert(n_scores > 0);

  if (n_scores == 0)
    return sas;

  assert(records[0].weaver_score != SRAlignmentSel::MISSING_SCORE);
  sas.best_pair_score = records[0].weaver_score - NO_PAIR_PENALTY;
  sas.log_sum_exp = add_log(sas.log_sum_exp, static_cast<double>(sas.best_pair_score) * copts.log_base);
  std::vector<int> best_scores(1, 0);

  for (int s{1}; s < n_scores; ++s)
  {
    assert(records[s].weaver_score != SRAlignmentSel::MISSING_SCORE);
    int const score{records[s].weaver_score - NO_PAIR_PENALTY};
    sas.log_sum_exp = add_log(sas.log_sum_exp, static_cast<double>(score) * copts.log_base);

    if (score >= sas.best_pair_score)
    {
      if (score > sas.best_pair_score)
        best_scores.clear();

      sas.second_best_pair_score = sas.best_pair_score;
      sas.best_pair_score = score;
      best_scores.push_back(s);
    }
  }

  sas.log_sum_exp = add_log(sas.log_sum_exp, 0.0);

#ifndef NDEBUG
  print_debug(_HERE_, " scaled best pair score=", static_cast<double>(sas.best_pair_score) * copts.log_base);
  print_debug(_HERE_, " sas.log_sum_exp=", sas.log_sum_exp);

  double mapq = -5.0 / std::log(10) *
                subtract_log(0.0, static_cast<double>(sas.best_pair_score) * copts.log_base - sas.log_sum_exp);

  if (n_scores > 1 /*&& !std::isinf(mapq)*/)
    print_debug(_HERE_, " mapq=", mapq);
#endif // NDEBUG

  assert(best_scores.size() > 0);
  int const max_score_i = select_an_index(best_scores, sas.best_pair_score);
  assert(max_score_i >= 0);
  assert(max_score_i < n_scores);

  if (is_read1)
  {
    sas.s1 = max_score_i;

    if (sas.second_best_pair_score != SRAlignmentSel::MISSING_SCORE)
      sas.other_read1_best_read_score = sas.second_best_pair_score + NO_PAIR_PENALTY;
  }
  else
  {
    sas.s2 = max_score_i;

    if (sas.second_best_pair_score != SRAlignmentSel::MISSING_SCORE)
      sas.other_read2_best_read_score = sas.second_best_pair_score + NO_PAIR_PENALTY;
  }

  return sas;
}

int get_mapq(int const best_read_score, int const read_length, int const best_pair_score, double const log_sum_exp)
{
  using namespace weaver;

  // unmapped reads get MQ=0
  if (best_read_score == SRAlignmentSel::MISSING_SCORE)
    return 0;

  double const identity = get_alignment_identity(best_read_score, read_length);
  assert(identity <= 1.0);
  assert(identity >= 0.0);

  Options const & copts = (*Options::const_instance());

  double const scaled_score = static_cast<double>(best_pair_score) * copts.log_base;
  assert(log_sum_exp >= scaled_score);

  // pre-MQ calculation. Check if infinity before determining the final MQ
  double const mapq_pre = subtract_log(0.0, scaled_score - log_sum_exp);

  if (std::isinf(mapq_pre))
    return 60;

  int mapq = std::round(-5.0 / std::log(10.0) * mapq_pre * std::pow(identity, copts.identity_penalty));

  // MQ is in [0,60]. All ambigous mappings are set to 0.
  return mapq >= 60 ? 60 : (mapq <= 3 ? 0 : mapq);
}

} // namespace

namespace weaver
{
int SRAlignmentSel::get_score_diff() const
{
  return second_best_pair_score == MISSING_SCORE ? 255 : std::min(255, best_pair_score - second_best_pair_score);
}

int SRAlignmentSel::get_mapq1(int const best_read1_score, int const read1_length) const
{
  return ::get_mapq(best_read1_score, read1_length, best_pair_score, log_sum_exp);
}

int SRAlignmentSel::get_mapq2(int const best_read2_score, int const read2_length) const
{
  return ::get_mapq(best_read2_score, read2_length, best_pair_score, log_sum_exp);
}

SRAlignmentSel select_alignment(GFA const & gfa,
                                T_icu const & icu,
                                std::vector<SRSeed> const & seeds1,
                                std::vector<SAMRecord> const & sam_records1,
                                std::vector<SRSeed> const & seeds2,
                                std::vector<SAMRecord> const & sam_records2)
{
  using namespace weaver;

  assert(seeds1.size() > 0 || seeds2.size() > 0);
  assert(seeds1.size() == sam_records1.size());
  assert(seeds2.size() == sam_records2.size());

  if (seeds1.empty())
    return select_unpaired_alignment(sam_records2, /*is_read1=*/false);

  if (seeds2.empty())
    return select_unpaired_alignment(sam_records1, /*is_read1=*/true);

  SRAlignmentSel sas;
  std::vector<int> scores;
  std::vector<int> best_pair_scores;
  int const n1 = seeds1.size();
  int const n2 = seeds2.size();
  int const n_scores = n1 * n2;
  scores.reserve(n_scores);
  Options const & copts = *(Options::const_instance());
  assert(copts.log_base >= 0);

  for (int s1{0}; s1 < n1; ++s1)
  {
    SRSeed const & seed1 = seeds1[s1];
    SAMRecord const & sam_record1 = sam_records1[s1];
    assert(sam_record1.alignment_score != SAMRecord::MISSING_TAG);
    assert(sam_record1.weaver_score != SAMRecord::MISSING_TAG);

    for (int s2{0}; s2 < n2; ++s2)
    {
      SRSeed const & seed2 = seeds2[s2];
      SAMRecord const & sam_record2 = sam_records2[s2];
      assert(sam_record2.alignment_score != SAMRecord::MISSING_TAG);
      assert(sam_record2.weaver_score != SAMRecord::MISSING_TAG);

      int score{sam_record1.weaver_score};
      assert(score != SRAlignment::MISSING_SCORE);
      static_assert(SAMRecord::MISSING_TAG == SRAlignment::MISSING_SCORE);

      if (!seed1.do_you_see_me(gfa, icu, seed2))
      {
        score += sam_record2.weaver_score - NO_PAIR_PENALTY;
      }
      else // seed1 sees seed2
      {
        assert(seed2.do_you_see_me(gfa, icu, seed1));
        score += sam_record2.weaver_score + PAIR_BONUS;
      }

      // check for underflow
      if (score < SRAlignmentSel::MISSING_SCORE)
        score = SRAlignmentSel::MISSING_SCORE;

      if (score >= sas.best_pair_score)
      {
        if (score > sas.best_pair_score)
          best_pair_scores.clear();

        sas.second_best_pair_score = sas.best_pair_score;
        sas.best_pair_score = score;
        best_pair_scores.push_back(scores.size());
      }
      else if (score > sas.second_best_pair_score)
      {
        sas.second_best_pair_score = score;
      }

      scores.push_back(score);
      sas.log_sum_exp = add_log(sas.log_sum_exp, static_cast<double>(score) * copts.log_base);
    }
  }

  // compare all scores against 0, as an approximation
  sas.log_sum_exp = add_log(sas.log_sum_exp, 0.0);

  assert(static_cast<int>(scores.size()) == n_scores);
  int const max_score_i = select_an_index(best_pair_scores, sas.best_pair_score);
  sas.s1 = max_score_i / n2;
  sas.s2 = max_score_i % n2;

#ifndef NDEBUG
  {
    double mapq = -5.0 / std::log(10) *
                  subtract_log(0.0, static_cast<double>(sas.best_pair_score) * copts.log_base - sas.log_sum_exp);

    if (not std::isinf(mapq))
      print_debug(_HERE_, " mapq=", mapq);
  }
#endif // NDEBUG

  assert(sas.s1 >= 0);
  assert(sas.s1 < n1);
  assert(sas.s2 >= 0);
  assert(sas.s2 < n2);

  if (sas.second_best_pair_score == SRAlignmentSel::MISSING_SCORE)
    return sas; // Unique match for both reads

  for (int s{0}; s < n_scores; ++s)
  {
    int const s1 = s / n2; // index for read 1 hit
    int const s2 = s % n2; // index for read 2 hit

    if (s1 != sas.s1) // read 1 is different from best
    {
      assert(s1 < static_cast<int>(sam_records1.size()));

      // update second best read 1 score
      sas.other_read1_best_read_score = std::max(sas.other_read1_best_read_score, sam_records1[s1].weaver_score);
    }

    if (s2 != sas.s2) // read 2 is different from best
    {
      assert(s2 < static_cast<int>(sam_records2.size()));

      // update second best read 2 score
      sas.other_read2_best_read_score = std::max(sas.other_read2_best_read_score, sam_records2[s2].weaver_score);
    }
  }

  return sas;
}

} // namespace weaver
