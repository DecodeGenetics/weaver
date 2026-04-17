/*!
 * @file sr_chain.hpp
 * @brief Implementation of the SRChain class.
 */

#include "sr_chain.hpp"

#include <cassert> // assert
#include <limits>  // std::numeric_limits
#include <numeric> // std::iota
#include <vector>  // std::vector

#include <paw/align/alignment_results.hpp>

#include "icu.hpp"
#include "logging.hpp"
#include "sr_alignment.hpp"
#include "sr_seed.hpp"

namespace weaver
{
void SRChain::check_score(GFA const & /*gfa*/,
                          T_icu const & /*icu*/,
                          int new_score, //
                          int new_seed_i1,
                          int new_seed_i2,
                          int new_score1,
                          int new_score2,
                          std::vector<SRSeed> const & /*seeds1*/,
                          std::vector<SRSeed> const & /*seeds2*/)
{
  int seed_i1 = hits.empty() ? -1 : hits[0].first;
  int seed_i2 = hits.empty() ? -1 : hits[0].second;

  if (new_score >= score)
  {
    // Check if the new score is the highest score
    if (new_score > score)
      hits.resize(0);

    // Check if read 1 score is now secondary
    if (seed_i1 != -1 && new_seed_i1 != seed_i1 && score1 >= score_secondary1)
    {
      if (score1 > score_secondary1)
      {
        score_secondary1 = score1; // new best secondary score for read 1
        seed_secondary1_i.resize(0);
      }

      seed_secondary1_i.emplace_back(seed_i1, seed_i2);
    }

    // Check if read 2 score is now secondary
    if (seed_i2 != -1 && new_seed_i2 != seed_i2 && score2 >= score_secondary2)
    {
      if (score2 > score_secondary2)
      {
        score_secondary2 = score2;
        seed_secondary2_i.resize(0);
      }

      seed_secondary2_i.emplace_back(seed_i1, seed_i2);
    }

    // Set a new best score
    hits.emplace_back(new_seed_i1, new_seed_i2);
    score1 = new_score1;
    score2 = new_score2;
    score = new_score;
  }
  else
  {
    // Check secondary scores of read 1
    if (new_seed_i1 != -1 && new_seed_i1 != seed_i1 && new_score1 >= score_secondary1)
    {
      if (new_score1 > score_secondary1)
      {
        score_secondary1 = new_score1; // new best secondary score
        seed_secondary1_i.resize(0);
      }

      seed_secondary1_i.emplace_back(new_seed_i1, new_seed_i2);
    }

    // Check secondary scores of read 2
    if (new_seed_i2 != -1 && new_seed_i2 != seed_i2 && new_score2 >= score_secondary2)
    {
      if (new_score2 > score_secondary2)
      {
        score_secondary2 = new_score2;
        seed_secondary2_i.resize(0);
      }

      seed_secondary2_i.emplace_back(new_seed_i1, new_seed_i2);
    }
  }
}

int SRChain::get_mapping_quality1() const
{
  if (score_secondary1 <= SRChain::MIN_SCORE)
    return 60;

  int const est_mq = ((score1 - score_secondary1) * 3) / 2;

  if (est_mq <= 0)
    return 0;

  if (est_mq >= 60)
    return 60;

  return est_mq;
}

int SRChain::get_mapping_quality2() const
{
  if (score_secondary2 <= SRChain::MIN_SCORE)
    return 60;

  int const est_mq = ((score2 - score_secondary2) * 3) / 2;

  if (est_mq <= 0)
    return 0;

  if (est_mq >= 60)
    return 60;

  return est_mq;
}

void SRChain::get_primary_and_secondary_index(std::pair<int, int> & primary_index,
                                              std::pair<int, int> & secondary_index,
                                              std::vector<SRSeed> const & s1,
                                              std::vector<SRSeed> const & s2)
{
  auto get_a_hit = [this, &s1, &s2](std::vector<std::pair<int, int>> const & p) -> std::pair<int, int>
  {
    std::size_t const n = p.size();

    if (n == 0)
      return std::pair<int, int>{-1, -1};

    if (n == 1)
      return p[0];

    // select pseudo randomly (but consistently) from the hits
    std::size_t sum = score;

    for (std::size_t h{0}; h < n; ++h)
    {
      std::pair<int, int> const & hit = p[h];

      if (hit.first != -1)
      {
        assert(hit.first < static_cast<int>(s1.size()));
        auto const & seed = s1[hit.first];
        sum += (seed.begin_ref_value + seed.end_ref_value * 2);
      }

      if (hit.second != -1)
      {
        assert(hit.second < static_cast<int>(s2.size()));
        auto const & seed = s2[hit.second];
        sum += (seed.begin_ref_value * 2 + seed.end_ref_value);
      }
    }

    return p[sum % n];
  };

  if (hits.size() > 0)
  {
    primary_index = get_a_hit(hits);
    secondary_index.first = get_a_hit(seed_secondary1_i).first;
    secondary_index.second = get_a_hit(seed_secondary2_i).second;
  }
}

std::pair<int, int> SRChain::get_primary_index(std::vector<SRSeed> const & s1, std::vector<SRSeed> const & s2) const
{
  std::size_t const n = hits.size();

  if (n == 0)
    return std::pair<int, int>(-1, -1);

  if (n == 1)
    return hits[0];

  // select pseudo randomly (but consistently) from the hits
  std::size_t sum{0};

  for (std::size_t h{0}; h < n; ++h)
  {
    std::pair<int, int> const & hit = hits[h];

    if (hit.first != -1)
    {
      assert(hit.first < static_cast<int>(s1.size()));
      auto const & seed = s1[hit.first];
      sum += (seed.begin_ref_value + seed.end_ref_value * 2);
    }

    if (hit.second != -1)
    {
      assert(hit.second < static_cast<int>(s2.size()));
      auto const & seed = s2[hit.second];
      sum += (seed.begin_ref_value * 2 + seed.end_ref_value);
    }
  }

  return hits[sum % n];
}

std::string SRChain::to_string() const
{
  std::string str;

  if (hits.empty())
    str += "i1 NA";
  else
    str += "i1 " + std::to_string(hits[0].first);

  if (hits.empty())
    str += ", i2 NA";
  else
    str += ", i2 " + std::to_string(hits[0].second);

  if (score == MIN_SCORE)
    str += ", score NA";
  else
    str += ", score " + std::to_string(score);

  if (score1 == MIN_SCORE)
    str += ", score1 NA";
  else
    str += ", score1 " + std::to_string(score1);

  if (score2 == MIN_SCORE)
    str += ", score2 NA";
  else
    str += ", score2 " + std::to_string(score2);

  if (score_secondary1 == MIN_SCORE)
    str += ", 2nd score1 NA";
  else
    str += ", 2nd score1 " + std::to_string(score_secondary1);

  if (score_secondary2 == MIN_SCORE)
    str += ", 2nd score2 NA";
  else
    str += ", 2nd score2 " + std::to_string(score_secondary2);

  return str;
}

SRChain get_chains(GFA const & gfa,
                   T_icu const & icu,
                   std::vector<SRSeed> const & seeds1,
                   std::vector<SRAlignment> const & alignments1,
                   std::vector<SRSeed> & seeds2,
                   std::vector<SRAlignment> const & alignments2)
{
  SRChain chains;

  if (seeds1.empty() && seeds2.empty())
    return chains;

  assert(seeds1.size() == alignments1.size());
  assert(seeds2.size() == alignments2.size());

  if (seeds1.empty())
  {
    for (int s2{0}; s2 < static_cast<int>(seeds2.size()); ++s2)
    {
      auto const score = alignments2[s2].score;
      assert(score != SRChain::MIN_SCORE);
      chains.check_score(gfa, icu, score, -1, s2, SRChain::MIN_SCORE, score, seeds1, seeds2);
    }

    return chains;
  }

  if (seeds2.empty())
  {
    for (int s1{0}; s1 < static_cast<int>(seeds1.size()); ++s1)
    {
      auto const score = alignments1[s1].score;
      assert(score != SRChain::MIN_SCORE);
      chains.check_score(gfa, icu, score, s1, -1, score, SRChain::MIN_SCORE, seeds1, seeds2);
    }

    return chains;
  }

  for (int s1{0}; s1 < static_cast<int>(seeds1.size()); ++s1)
  {
    auto const & seed1 = seeds1[s1];
    auto const score1 = alignments1[s1].score;
    assert(score1 != SRChain::MIN_SCORE);

    for (int s2{0}; s2 < static_cast<int>(seeds2.size()); ++s2)
    {
      auto const & seed2 = seeds2[s2];
      auto const score2 = alignments2[s2].score;
      assert(score2 != SRChain::MIN_SCORE);
      int extra_score{0};
      int extra_score_each_read{0};

      if (seed1.do_you_see_me(gfa, icu, seed2))
      {
        assert(seed2.do_you_see_me(gfa, icu, seed1));
        extra_score = 42;
        extra_score_each_read = 10;
      }

      chains.check_score(gfa, //
                         icu,
                         score1 + score2 + extra_score,
                         s1,
                         s2,
                         score1 + extra_score_each_read,
                         score2 + extra_score_each_read,
                         seeds1,
                         seeds2);
    }
  }

  assert(chains.score != SRChain::MIN_SCORE);
  assert(chains.score1 != SRChain::MIN_SCORE);
  assert(chains.score2 != SRChain::MIN_SCORE);
  return chains;
}

} // namespace weaver
