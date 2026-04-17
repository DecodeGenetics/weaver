#include "haplotype_stats.hpp"

#include <cassert> // assert
#include <limits>  // std::numeric_limits
#include <vector>  // std::vector

#include "alignment_utils.hpp" // weaver::add_log
#include "edit.hpp"            // weaver::EditCalls
#include "edit_stats.hpp"      // weaver::EditStats
#include "logging.hpp"         // print_info

namespace weaver
{
HaplotypeStats::HaplotypeStats(const MMI::T_haplotypes * _haps_ptr) : haps_ptr(_haps_ptr)
{
  set_haps_ptr(haps_ptr);
}

void HaplotypeStats::set_haps_ptr(MMI::T_haplotypes const * new_haps_ptr)
{
  assert(new_haps_ptr != nullptr);
  haps_ptr = new_haps_ptr;

  MMI::T_haplotypes const & haps = *haps_ptr;
  snid_haps_stats.resize(haps.size());

  for (int c{0}; c < static_cast<int>(haps.size()); ++c)
  {
    std::vector<EditCalls> const & edit_calls = haps[c];
    std::vector<EditStats> & estats = snid_haps_stats[c];
    estats.reserve(edit_calls.size()); // will be the same size eventually

    for (auto const & edit_call : edit_calls)
      estats.emplace_back(edit_call);
  }
}

bool HaplotypeStats::operator==(HaplotypeStats const & other_hap_stats) const
{
  if (haps_ptr != other_hap_stats.haps_ptr)
    return false;

  if (snid_haps_stats.size() != other_hap_stats.snid_haps_stats.size())
    return false;

  int const n_snid{static_cast<int>(snid_haps_stats.size())};

  for (int c{0}; c < n_snid; ++c)
  {
    if (snid_haps_stats[c].size() != other_hap_stats.snid_haps_stats[c].size())
      return false;

    std::vector<EditStats> const & edits = snid_haps_stats[c];
    std::vector<EditStats> const & other_edits = other_hap_stats.snid_haps_stats[c];
    int const n_edits{static_cast<int>(edits.size())};

    for (int e{0}; e < n_edits; ++e)
    {
      if (edits[e] != other_edits[e])
        return false;
    }
  }

  return true;
}

void HaplotypeStats::run_hap_weighting()
{
  //
  assert(haps_ptr != nullptr);
  MMI::T_haplotypes const & haps = *haps_ptr;

  for (int c{0}; c < static_cast<int>(haps.size()); ++c)
  {
    assert(haps.size() == snid_haps_stats.size());
    assert(c < static_cast<int>(snid_haps_stats.size()));

    std::vector<EditCalls> const & edit_calls = haps[c];

    if (edit_calls.size() == 0)
      continue; // no edits on snid

    std::vector<EditStats> & edit_stats = snid_haps_stats[c];
    int const num_nonref_haplotypes{static_cast<int>(edit_calls[0].calls.size())};
    int const num_haplotypes{num_nonref_haplotypes + 1};
    assert(edit_calls.size() == edit_stats.size());
    int const n_edits{static_cast<int>(edit_calls.size())};

    // initialization
    double constexpr r = 1.26e-8;
    double constexpr N_e = 25000;
    // double constexpr eps = -20.0;

    // initialization for forward algorithm
    int pos{-1};
    // std::vector<std::vector<double>> alphas_betas(n_edits);
    std::vector<double> alpha(num_haplotypes, 0.0);

    for (int f{0}; f < n_edits; ++f)
    {
      EditCalls const & ec = edit_calls[f];
      EditStats & es = edit_stats[f];
      uint8_t const hom_call = es.get_hom_call();

      if (hom_call == Edit::MISSING_CALL)
      {
        // no new information, just copy
        es.alpha_beta = alpha;
        es.alpha_beta_log_sum = log_sum(es.alpha_beta);
        // es.alpha = alpha;
        continue;
      }

      if (pos == ec.pos)
      {
        // variant at the same pos, reuse alpha for now
        es.alpha_beta = alpha;
        es.alpha_beta_log_sum = log_sum(es.alpha_beta);
        // es.alpha = alpha;
        continue;
      }

      // print_info(_HERE_, " edit,no_edit read counts = ", es.edit_read_count, " ", es.no_edit_read_count);
      // print_info(_HERE_, " hom_call = ", static_cast<int>(hom_call), " at edit: ", ec.to_string());

      double const d = pos < 0 ? 0.0 : 4. * static_cast<double>(ec.pos - pos) * r * N_e;
      assert(d >= 0.0);
      double const interm = std::exp(-d / static_cast<double>(num_haplotypes));
      double const p = (1 - interm) / static_cast<double>(num_haplotypes);
      double const q = interm + p;
      double const lp = pos < 0 ? std::numeric_limits<double>::lowest() : std::log(p);
      double const lq = pos < 0 ? 0.0 : std::log(q);
      // print_info(_HERE_, " lq, lp = ", lq, ", ", lp);

      assert(ec.pos > pos);
      pos = ec.pos;

      double log_sum_i_call_0{std::numeric_limits<double>::lowest()};
      double log_sum_i_call_1{std::numeric_limits<double>::lowest()};
      double obs_i_call_0{0.0};
      double obs_i_call_1{0.0};

      if (hom_call == 0)
        obs_i_call_1 = es.get_hom_no_edit_eps();
      else if (hom_call == 1)
        obs_i_call_0 = es.get_hom_edit_eps();

      for (int s{0}; s < num_haplotypes; ++s)
      {
        add_log_ref(log_sum_i_call_0, alpha[s] + lp + obs_i_call_0);
        add_log_ref(log_sum_i_call_1, alpha[s] + lp + obs_i_call_1);
      }

      for (int i{0}; i < num_haplotypes; ++i)
      {
        uint8_t const i_call = i == 0 ? 0 : ec.calls[i - 1];

        if (i_call == 0)
        {
          double const new_log_sum = subtract_log(log_sum_i_call_0, alpha[i] + lp + obs_i_call_0);

          if (std::isinf(new_log_sum))
            alpha[i] = alpha[i] + lq + obs_i_call_0;
          else
            alpha[i] = add_log(new_log_sum, alpha[i] + lq + obs_i_call_0);
        }
        else if (i_call == 1)
        {
          double const new_log_sum = subtract_log(log_sum_i_call_1, alpha[i] + lp + obs_i_call_1);

          if (std::isinf(new_log_sum))
            alpha[i] = alpha[i] + lq + obs_i_call_1;
          else
            alpha[i] = add_log(new_log_sum, alpha[i] + lq + obs_i_call_1);
        }
        else
        {
          assert(i_call == Edit::MISSING_CALL);

          // take lower likelihood call
          alpha[i] = std::min(log_sum_i_call_0, log_sum_i_call_1);
        }
      }

      // normalize alpha
      auto max_it = std::max_element(alpha.begin(), alpha.end());
      assert(max_it != alpha.end());
      double const max_val = *max_it;

      for (int a{0}; a < static_cast<int>(alpha.size()); ++a)
        alpha[a] -= max_val;

      // print_info(_HERE_, " here is alpha: ");
      //
      // for (int a{0}; a < static_cast<int>(alpha.size()); ++a)
      //   print_info(a, " ", alpha[a]);

      es.alpha_beta = alpha;
      es.alpha_beta_log_sum = log_sum(es.alpha_beta);
      // es.alpha = alpha;
    }

    // initialization for reverse algorithm
    pos = std::numeric_limits<int>::max();
    // std::vector<std::vector<double>> betas(n_edits);
    std::vector<double> beta(num_haplotypes, 0.0);

    for (int f{n_edits - 2}; f >= 0; --f)
    {
      EditCalls const & ec = edit_calls[f];
      EditStats & es = edit_stats[f];
      EditStats const & es_next = edit_stats[f + 1];
      uint8_t const hom_call = es_next.get_hom_call();

      if (hom_call == Edit::MISSING_CALL)
      {
        // no new information, just copy
        for (int h{0}; h < num_haplotypes; ++h)
          es.alpha_beta[h] += beta[h];

        es.alpha_beta_log_sum = log_sum(es.alpha_beta);
        continue;
      }

      if (pos == ec.pos)
      {
        // variant at the same pos, reuse beta for now
        for (int h{0}; h < num_haplotypes; ++h)
          es.alpha_beta[h] += beta[h];

        es.alpha_beta_log_sum = log_sum(es.alpha_beta);
        continue;
      }

      // print_info(_HERE_, " edit,no_edit read counts = ", es.edit_read_count, " ", es.no_edit_read_count);
      // print_info(_HERE_, " hom_call = ", static_cast<int>(hom_call), " at edit: ", ec.to_string());

      assert(pos >= ec.pos);
      double const d = pos == std::numeric_limits<int>::max() ? 0.0 : 4.0 * static_cast<double>(pos - ec.pos) * r * N_e;
      assert(d >= 0.0);
      double const interm = std::exp(-d / static_cast<double>(num_haplotypes));
      double const p = (1 - interm) / static_cast<double>(num_haplotypes);
      double const q = interm + p;
      double const lp = pos == std::numeric_limits<int>::max() ? std::numeric_limits<double>::lowest() : std::log(p);
      double const lq = pos == std::numeric_limits<int>::max() ? 0.0 : std::log(q);
      // print_info(_HERE_, " lq, lp = ", lq, ", ", lp);

      assert(ec.pos < pos);
      pos = ec.pos;

      double log_sum_i_call_0{std::numeric_limits<double>::lowest()};
      double log_sum_i_call_1{std::numeric_limits<double>::lowest()};
      double obs_i_call_0{0.0}; // = hom_call != 1 ? 0.0 : es.get_hom_no_edit_eps();
      double obs_i_call_1{0.0}; // = hom_call != 0 ? 0.0 : es.get_hom_edit_eps();

      if (hom_call == 0)
        obs_i_call_1 = es.get_hom_no_edit_eps();
      else if (hom_call == 1)
        obs_i_call_0 = es.get_hom_edit_eps();

      for (int s{0}; s < num_haplotypes; ++s)
      {
        add_log_ref(log_sum_i_call_0, beta[s] + lp + obs_i_call_0);
        add_log_ref(log_sum_i_call_1, beta[s] + lp + obs_i_call_1);
      }

      for (int i{0}; i < num_haplotypes; ++i)
      {
        uint8_t const i_call = i == 0 ? 0 : ec.calls[i - 1];

        if (i_call == 0)
        {
          double const new_log_sum = subtract_log(log_sum_i_call_0, beta[i] + lp + obs_i_call_0);

          if (std::isinf(new_log_sum))
            beta[i] = beta[i] + lq + obs_i_call_0;
          else
            beta[i] = add_log(new_log_sum, beta[i] + lq + obs_i_call_0);
        }
        else if (i_call == 1)
        {
          double const new_log_sum = subtract_log(log_sum_i_call_1, beta[i] + lp + obs_i_call_1);

          if (std::isinf(new_log_sum))
            beta[i] = beta[i] + lq + obs_i_call_1;
          else
            beta[i] = add_log(new_log_sum, beta[i] + lq + obs_i_call_1);
        }
        else
        {
          assert(i_call == Edit::MISSING_CALL);

          // take lower likelihood call
          beta[i] = std::min(log_sum_i_call_0, log_sum_i_call_1);
        }
      }

      // normalize beta
      auto max_it = std::max_element(beta.begin(), beta.end());
      assert(max_it != beta.end());
      double const max_val = *max_it;

      for (int b{0}; b < num_haplotypes; ++b)
        beta[b] -= max_val;

      print_info(_HERE_, " here is beta: ");

      for (int b{0}; b < num_haplotypes; ++b)
        print_info(b, " ", beta[b]);

      for (int h{0}; h < num_haplotypes; ++h)
        es.alpha_beta[h] += beta[h];

      es.alpha_beta_log_sum = log_sum(es.alpha_beta);
    } // reverse algorithm ends here

    // calculate the edit probability
    /*
    for (int f{0}; f < n_edits; ++f)
    {
      EditCalls const & ec = edit_calls[f];
      EditStats & es = edit_stats[f];
      std::vector<double> const & alpha_beta = alphas_betas[f];
      es.alpha_beta_logsum = log_sum(es.alpha_beta);
      double edit_logsum{std::numeric_limits<double>::lowest()};

      for (int i{0}; i < num_haplotypes; ++i)
      {
        if (i > 0 && ec.calls[i - 1] == 1)
          add_log_ref(edit_logsum, alpha_beta[i]);
      }

      es.edit_prob = std::exp(edit_logsum - alpha_beta_logsum);
      print_debug(_HERE_, " edit: ", ec.to_string(), " has allele freq = ", es.get_edit_frequency(), " ", es.edit_prob);
    }
    */
  }

  this->is_using_stats = true;
}

/*
void HaplotypeStats::run_hmm()
{
  assert(haps_ptr != nullptr);
  MMI::T_haplotypes const & haps = *haps_ptr;

  for (int c{0}; c < static_cast<int>(haps.size()); ++c)
  {
    assert(haps.size() == snid_haps_stats.size());
    assert(c < static_cast<int>(snid_haps_stats.size()));

    std::vector<EditCalls> const & edit_calls = haps[c];

    if (edit_calls.size() == 0)
      continue; // no edits on snid

    std::vector<EditStats> & edit_stats = snid_haps_stats[c];
    int const num_nonref_haplotypes{static_cast<int>(edit_calls[0].calls.size())};
    int const num_haplotypes{num_nonref_haplotypes + 1};
    assert(edit_calls.size() == edit_stats.size());
    int const n_edits{static_cast<int>(edit_calls.size())};

    // initialization
    double constexpr r = 1.26e-8; // -8 for cM?
    double constexpr N_e = 25000;
    // double d{};
    // double p{}; // prob of transition
    // double q{}; // prob of no transition

    // initialization for forward algorithm
    int const n_pairs = (num_haplotypes * (num_haplotypes + 1)) / 2;
    std::vector<double> alpha_log(n_pairs, 0.0);
    std::vector<double> alpha_new(n_pairs, std::numeric_limits<double>::lowest());

    int pos{-1};

    // forward
    for (int f{0}; f < n_edits;) // no increment
    {
      EditCalls const & ec = edit_calls[f];
      EditStats & es = edit_stats[f];

      assert(f == 0 || pos >= 0);
      assert(ec.pos > pos);
      double const d = 4 * (ec.pos - pos) * r * N_e;
      double const interm = std::exp(-d / static_cast<double>(num_haplotypes));
      double const p = (1 - interm) / static_cast<double>(num_haplotypes);
      double const q = interm + p;
      double const lqq = f == 0 ? 0 : std::log(q * q);
      double const lpq = f == 0 ? std::numeric_limits<double>::lowest() : std::log(p * q);
      double const lpp = f == 0 ? std::numeric_limits<double>::lowest() : std::log(p * p);
      assert(lqq >= lpq);
      pos = ec.pos;
      // es.post_prob_log.resize(num_haplotypes, std::numeric_limits<double>::lowest());

      int f_end{f + 1};

      while (f_end < static_cast<int>(edit_calls.size()) && //
             ec.pos == edit_calls[f_end].pos &&             //
             ec.type == edit_calls[f_end].type)
      {
        ++f_end;
      }

      // do biallelic case first
      if (f + 1 == f_end)
      {
        print_debug(_HERE_, " forward edit with ", ec.to_string());

        std::vector<double> obs_log(n_pairs, 0.0); // observed haplotypes probability in log scale

        {
          int idx{0};

          for (int j{0}; j < num_haplotypes; ++j)
          {
            auto const call_j = j == 0 ? 0 : ec.calls[j - 1];

            for (int i{0}; i <= j; ++i)
            {
              auto const call_i = i == 0 ? 0 : ec.calls[i - 1];

              if (call_i == Edit::MISSING_CALL && call_j == Edit::MISSING_CALL)
              {
                // print_info(_HERE_, idx, " is missing call.");
                obs_log[idx++] += (es.edit_read_count + es.no_edit_read_count) * READ_COUNT_TO_LOG;
                continue;
              }

              if ((call_i == 0 && call_j != 1) || (call_i != 1 && call_j == 0))
              {
                // print_info(_HERE_,
                //            " ",
                //            idx,
                //            " is hom no edit call ",
                //            es.edit_read_count * READ_COUNT_TO_LOG,
                //            " depths=",
                //            es.no_edit_read_count,
                //            ",",
                //            es.edit_read_count);

                obs_log[idx++] += es.edit_read_count * READ_COUNT_TO_LOG; // hom no edit call
              }
              else if ((call_i == 1 && call_j != 0) || (call_i != 0 && call_j == 1))
              {
                // print_info(_HERE_,
                //            " ",
                //            idx,
                //            " is hom call ",
                //            es.no_edit_read_count * READ_COUNT_TO_LOG,
                //            " depths=",
                //            es.no_edit_read_count,
                //            ",",
                //            es.edit_read_count);

                obs_log[idx++] += es.no_edit_read_count * READ_COUNT_TO_LOG; // hom edit call
              }
              else
              {
                // print_info(_HERE_, " ", idx, " is het call 0 depths=", es.no_edit_read_count, ",",
                // es.edit_read_count);
                ++idx;
              }
            }
          }

          assert(idx == n_pairs);
        }

        {
          auto max_it = std::max_element(alpha_log.begin(), alpha_log.end());
          assert(max_it != alpha_log.end());
          double const max_alpha{*max_it};

          for (int i{0}; i < n_pairs; ++i)
            alpha_new[i] = alpha_log[i] - max_alpha;
        }

        std::vector<double> log_values_to_sum{};

        // print_info(_HERE_, " lqq, lpq, lpp = ", lqq, ",", lpq, ",", lpp);
        // print_info(_HERE_, " here's alpha before ", ec.to_string());
        //
        // for (int a{0}; a < std::min(15, n_pairs); ++a)
        //   print_info(a, " ", alpha_log[a]);
        //
        // print_info(_HERE_, " here's obs_log with ", ec.to_string());
        //
        // for (int a{0}; a < std::min(15, n_pairs); ++a)
        //   print_info(a, " ", obs_log[a]);

        // first, update all alpha_new values with values from the previous idx
        {
          int idx = 0;

          for (int j{0}; j < num_haplotypes; ++j)
          {
            for (int i{0}; i <= j; ++i)
            {
              assert(idx < static_cast<int>(alpha_new.size()));
              assert(idx < static_cast<int>(obs_log.size()));

              alpha_new[idx] += obs_log[idx] + lqq;
              ++idx;
            }
          }

          assert(idx == n_pairs);
        }

        // print_debug(_HERE_, " here's alpha new after same hap:");
        //
        // for (int a{0}; a < std::min(15, n_pairs); ++a)
        //   print_debug(a, " ", alpha_new[a]);

        {
          int idx = 0;

          // then update all cells with one difference in haplotype
          for (int j{0}; j < num_haplotypes; ++j)
          {
            // heterozygous cases
            for (int i{0}; i <= j; ++i, ++idx)
            {
              assert(idx < n_pairs);

              // k == i && l < j loop
              for (int l{0}; l < i; ++l)
              {
                int idx2 = (i + 1) * i / 2 + l;
                add_log_ref(alpha_new[idx], alpha_log[idx2] + obs_log[idx2] + lpq);
              }

              // k == i && l > j loop
              for (int l{i}; l < num_haplotypes; ++l)
              {
                int idx2 = (l + 1) * l / 2 + i;
                add_log_ref(alpha_new[idx], alpha_log[idx2] + obs_log[idx2] + lpq);
              }

              // l == j && k < i
              // idx2 = (l+1)*l/2+i
              // maybeTODO optimize by splitting in parts to avoid conditionals
              for (int l{0}; l < num_haplotypes; ++l)
              {
                for (int k{0}; k <= l; ++k, ++idx2)
                {
                  assert(idx2 < n_pairs);

                  if (idx == idx2)
                  {
                    assert(i == k && j == l);
                    continue;
                    // add_log_ref(alpha_new[idx], alpha_log[idx2] + obs_log[idx2] + lqq);
                  }

                  if (i == k || j == l || i == l || j == k)
                  {
                    // print_info("i,j,k,l ", i, ",", j, ",", k, ",", l, " exp ", alpha_log[idx2] + obs_log[idx2] +
                    // lpq);
                    add_log_ref(alpha_new[idx], alpha_log[idx2] + obs_log[idx2] + lpq);
                  }
                  else
                  {
                    // print_info("i,j,k,l ", i, ",", j, ",", k, ",", l, " exp ", alpha_log[idx2] + obs_log[idx2] +
                    // lpp);
                    add_log_ref(alpha_new[idx], alpha_log[idx2] + obs_log[idx2] + lpp);
                  }
                }
              }

if ((j > 0 && ec.calls[j - 1] == 1) || (i > 0 && ec.calls[i - 1] == 1))
  add_log_ref(es.alpha_edit_prob_log, alpha_new[idx]);
}
}

assert(idx == n_pairs);
}

es.alpha_beta_log_sum = log_sum(alpha_new);

// print_info(_HERE_, " alpha_edit_prob_log ", es.alpha_edit_prob_log);
// print_info(_HERE_, " alpha_beta_log_sum = ", es.alpha_beta_log_sum);
// print_info(_HERE_, " ", std::exp(es.alpha_edit_prob_log - es.alpha_beta_log_sum));

std::swap(alpha_log, alpha_new);

f = f_end;
} // biallelic case ends
}
}

this->is_using_stats = true;
}
*/

// void update_pair_stats(HaplotypeStats & haplotype_stats,
//                        std::vector<int> const & edits1,
//                        std::vector<int> const & edits2,
//                        int const snid)
// {
//   if (snid >= static_cast<int>(haplotype_stats.snid_haps_stats.size()))
//     return;
//
//   std::vector<EditStats> & estats = haplotype_stats.snid_haps_stats[snid];
//   std::vector<int> pair_edits(edits1);
//   pair_edits.insert(pair_edits.end(), edits2.begin(), edits2.end());
//   std::sort(pair_edits.begin(), pair_edits.end(), is_less_edit_index);
//   int const n_edits{static_cast<int>(pair_edits.size())};
//
//   for (int i{0}; i < n_edits; ++i)
//   {
//     int const e = pair_edits[i];
//
//     if (e >= 0)
//     {
//       assert(e < static_cast<int>(estats.size()));
//       print_debug(_HERE_, " YES edit at ", e);
//       ++(estats[e].edit_read_count);
//     }
//     else // e < 0
//     {
//       int const no_e = -e - 1;
//       assert(e < static_cast<int>(estats.size()));
//       print_debug(_HERE_, " NO edit at ", no_e);
//       ++(estats[no_e].no_edit_read_count);
//     }
//   }
// }

void update_stats(HaplotypeStats & haplotype_stats, std::vector<int> const & edits, int const snid)
{
  if (snid >= static_cast<int>(haplotype_stats.snid_haps_stats.size()))
    return;

  std::vector<EditStats> & estats = haplotype_stats.snid_haps_stats[snid];
  int const n_edits{static_cast<int>(edits.size())};

  for (int i{0}; i < n_edits; ++i)
  {
    int const e = edits[i];

    if (e >= 0)
    {
      assert(e < static_cast<int>(estats.size()));
      print_debug(_HERE_, " YES edit at ", e);
      ++(estats[e].edit_read_count);
    }
    else // e < 0
    {
      int const no_e = -e - 1;
      assert(e < static_cast<int>(estats.size()));
      print_debug(_HERE_, " NO edit at ", no_e);
      ++(estats[no_e].no_edit_read_count);
    }
  }
}

HaplotypeStats merge_haplotype_stats(std::vector<HaplotypeStats> & p_hap_stats)
{
  print_debug(_HERE_, " merging haplotype stats from ", p_hap_stats.size(), " haplotypes.");

  assert(p_hap_stats.size() > 0);
  HaplotypeStats merged_stats(std::move(p_hap_stats[0]));
  int const n_threads{static_cast<int>(p_hap_stats.size())};

  for (int t{1}; t < n_threads; ++t)
  {
    HaplotypeStats const & thread_stats = p_hap_stats[t];
    assert(merged_stats.snid_haps_stats.size() == p_hap_stats[t].snid_haps_stats.size());
    int const n_snid{static_cast<int>(merged_stats.snid_haps_stats.size())};

    for (int c{0}; c < n_snid; ++c)
    {
      std::vector<EditStats> & merged_into = merged_stats.snid_haps_stats[c];
      std::vector<EditStats> const & merged_from = thread_stats.snid_haps_stats[c];
      assert(merged_into.size() == merged_from.size());
      int const n_edits{static_cast<int>(merged_into.size())};

      for (int e{0}; e < n_edits; ++e)
        merged_into[e].merge_with(merged_from[e]);
    }
  }

  return merged_stats;
}

void merge_and_then_scatter_haplotype_stats(std::vector<HaplotypeStats> & p_hap_stats)
{
  int const threads{static_cast<int>(p_hap_stats.size())};
  assert(threads > 0);
  HaplotypeStats new_hap_stats = merge_haplotype_stats(p_hap_stats);

  // Run the HMM after merging the haplotype stats
  // new_hap_stats.run_hmm();
  new_hap_stats.run_hap_weighting();

  assert(static_cast<int>(p_hap_stats.size()) == threads);
  // Give every thread the results
  p_hap_stats.assign(threads, new_hap_stats);
  assert(static_cast<int>(p_hap_stats.size()) == threads); // it replaces the old values

  // all values are equal
  assert(std::equal(p_hap_stats.begin(),     //
                    p_hap_stats.end() - 1,   //
                    p_hap_stats.begin() + 1, //
                    p_hap_stats.end()));     //
}

} // namespace weaver
