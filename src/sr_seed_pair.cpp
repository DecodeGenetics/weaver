#include "sr_seed_pair.hpp"

#include <cassert>
#include <limits>
#include <string>
#include <vector>

#include "gfa.hpp"
#include "icu.hpp"
#include "logging.hpp"

namespace
{
bool larger_chain_score(weaver::SRSeedPair const & a, weaver::SRSeedPair const & b)
{
  assert(a.score != weaver::SRSeedPair::MIN_SCORE);
  assert(b.score != weaver::SRSeedPair::MIN_SCORE);
  return a.score > b.score;
}

std::vector<weaver::SRSeedPair> get_best_seed_pairs_unpaired_read1(std::vector<weaver::SRSeed> const & seeds1)
{
  std::vector<weaver::SRSeedPair> chains;
  int const n1{static_cast<int>(seeds1.size())};

  if (n1 > 0)
  {
    std::vector<int> order1 = get_seed_est_score_sorted_order_indices(seeds1);
    int const stop_criteria = std::min(n1, weaver::SRSeedPair::MAX_CHAINS);

    for (int o1{0}; o1 < stop_criteria; ++o1)
      chains.push_back({static_cast<int>(order1[o1]), -1, seeds1[order1[o1]].get_est_score()});
  }

  return chains;
}

std::vector<weaver::SRSeedPair> get_best_seed_pairs_unpaired_read2(std::vector<weaver::SRSeed> const & seeds2)
{
  std::vector<weaver::SRSeedPair> chains;
  int const n2{static_cast<int>(seeds2.size())};

  if (n2 > 0)
  {
    std::vector<int> order2 = get_seed_est_score_sorted_order_indices(seeds2);
    int const stop_criteria = std::min(n2, weaver::SRSeedPair::MAX_CHAINS);

    for (int o2{0}; o2 < stop_criteria; ++o2)
      chains.push_back({-1, static_cast<int>(order2[o2]), seeds2[order2[o2]].get_est_score()});
  }

  return chains;
}

} // namespace

namespace weaver
{
std::string SRSeedPair::to_string() const
{
  std::string str;

  if (seed_i1 == -1)
    str += "i1 NA";
  else
    str += "i1 " + std::to_string(seed_i1);

  if (seed_i2 == -1)
    str += ", i2 NA";
  else
    str += ", i2 " + std::to_string(seed_i2);

  if (score == std::numeric_limits<int>::min())
    str += ", score NA";
  else
    str += ", score " + std::to_string(score);

  return str;
}

int get_seed_pair_score(GFA const & gfa, T_icu const & icu, SRSeed const & seed1, SRSeed const & seed2)
{
  int const score = seed1.get_est_score() + seed2.get_est_score();

  // seed1> <seed2
  if (seed1.do_you_see_me(gfa, icu, seed2))
  {
    assert(seed2.do_you_see_me(gfa, icu, seed1));
    return score + 60;
  }
  else
  {
    return score;
  }
}

std::vector<SRSeedPair> get_best_seed_pairs(GFA const & gfa,                    // graph
                                            T_icu const & icu,                  // The "I see you" index
                                            std::vector<SRSeed> const & seeds1, // seeds of first sequence
                                            std::vector<SRSeed> const & seeds2) // seeds of second sequence
{
  if (seeds2.empty())
    return get_best_seed_pairs_unpaired_read1(seeds1);

  if (seeds1.empty())
    return get_best_seed_pairs_unpaired_read2(seeds2);

  std::vector<SRSeedPair> chains;
  std::vector<int> order1 = get_seed_est_score_sorted_order_indices(seeds1);
  std::vector<int> order2 = get_seed_est_score_sorted_order_indices(seeds2);

  assert(seeds1.size() == order1.size());
  assert(seeds2.size() == order2.size());

  int const n1 = std::min(2048, static_cast<int>(seeds1.size()));
  int const n2 = std::min(2048, static_cast<int>(seeds2.size()));
  int const max_est_score1{seeds1[order1[0]].get_est_score()};
  int const max_est_score2{seeds2[order2[0]].get_est_score()};

  int best_chain_score{std::numeric_limits<int>::min()};
  int o1{0};

  // loop over best seed1 scores
  {
    int num_checks{0};

    for (/*none*/; o1 < n1; ++o1)
    {
      if (num_checks >= SRSeedPair::MAX_CHAINS_CHECKED_HARD_LIMIT ||
          (num_checks >= SRSeedPair::MAX_CHAINS_CHECKED && seeds1[order1[o1]].get_est_score() < max_est_score1))
      {
        break;
      }

      for (int o2{0}; o2 < n2; ++o2)
      {
        ++num_checks;
        int const chain_score = get_seed_pair_score(gfa, icu, seeds1[order1[o1]], seeds2[order2[o2]]);

        if ((chain_score + SRSeedPair::MAX_CHAIN_SCORE_DIFF) < best_chain_score)
          continue;
        else if (chain_score > best_chain_score)
          best_chain_score = chain_score;

        chains.emplace_back(order1[o1], order2[o2], chain_score);
      }
    }
  }

  // loop over best seed2 scores, start checking seeds1 at o1
  {
    assert(o1 > 0);
    int const init_o1{o1};
    int num_checks{init_o1};

    for (int o2{0}; o2 < n2; ++o2)
    {
      if (num_checks >= SRSeedPair::MAX_CHAINS_CHECKED_HARD_LIMIT ||
          (num_checks >= SRSeedPair::MAX_CHAINS_CHECKED && seeds2[order2[o2]].get_est_score() < max_est_score2))
      {
        break;
      }

      for (int o1_for2{init_o1}; o1_for2 < n1; ++o1_for2)
      {
        ++num_checks;
        int const chain_score = get_seed_pair_score(gfa, icu, seeds1[order1[o1_for2]], seeds2[order2[o2]]);

        if ((chain_score + SRSeedPair::MAX_CHAIN_SCORE_DIFF) < best_chain_score)
          continue;
        else if (chain_score > best_chain_score)
          best_chain_score = chain_score;

        chains.emplace_back(order1[o1_for2], order2[o2], chain_score);
      }
    }
  }

  // Sort chains, largest scores first
  std::sort(chains.begin(), chains.end(), larger_chain_score);

  if (static_cast<int>(chains.size()) > SRSeedPair::MAX_CHAINS)
    chains.resize(SRSeedPair::MAX_CHAINS);

  assert(!chains.empty());
  assert(chains[0].score == best_chain_score);

  // Then remove any seed pair hits that have too low of a score
  auto find_it = std::find_if(chains.begin(),
                              chains.end(),
                              [best_chain_score](SRSeedPair const & chain)
                              { return chain.score + SRSeedPair::MAX_CHAIN_SCORE_DIFF < best_chain_score; });

  chains.resize(std::distance(chains.begin(), find_it));
  assert(static_cast<int>(chains.size()) <= SRSeedPair::MAX_CHAINS);
  return chains;
}

void filter_sr_seeds(std::vector<SRSeed> & seeds1,
                     std::vector<SRSeed> & seeds2,
                     std::vector<SRSeedPair> const & best_seed_pairs)
{
#ifndef NDEBUG
  // for (auto const & best_seed_pair : best_seed_pairs)
  //   print_debug(_HERE_, " one of the best seed pair=", best_seed_pair.to_string());
#endif // NDEBUG

  if (seeds1.size() > 1)
  {
    int new_e{0};

    for (int e1{0}; e1 < static_cast<int>(seeds1.size()); ++e1)
    {
      // it is a nested loop but both loops are very small
      bool const found_any_e1 = std::any_of(best_seed_pairs.begin(),
                                            best_seed_pairs.end(),
                                            [e1](SRSeedPair const & sp) { return sp.seed_i1 == e1; });

      if (e1 > new_e)
        seeds1[new_e] = std::move(seeds1[e1]);

      if (found_any_e1)
        ++new_e;
    }

    seeds1.resize(new_e);
    assert(seeds1.size() > 0);
  }

  if (seeds2.size() > 1)
  {
    int new_e{0};

    for (int e2{0}; e2 < static_cast<int>(seeds2.size()); ++e2)
    {
      // it is a nested loop but both loops are very small
      bool const found_any_e2 = std::any_of(best_seed_pairs.begin(),
                                            best_seed_pairs.end(),
                                            [e2](SRSeedPair const & sp) { return sp.seed_i2 == e2; });

      if (e2 > new_e)
        seeds2[new_e] = std::move(seeds2[e2]);

      if (found_any_e2)
        ++new_e;
    }

    seeds2.resize(new_e);
    assert(seeds2.size() > 0);
  }
}

void get_sr_chains_and_filter(GFA const & gfa,
                              T_icu const & icu,
                              std::vector<SRSeed> & seeds1,
                              std::vector<SRSeed> & seeds2)
{
  std::vector<SRSeedPair> best_seed_pairs = get_best_seed_pairs(gfa, icu, seeds1, seeds2);
  filter_sr_seeds(seeds1, seeds2, best_seed_pairs);
}

} // namespace weaver
