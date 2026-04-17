#pragma once

#include <string>
#include <vector>

#include "edit.hpp"
#include "edit_stats.hpp"
#include "mmi.hpp"
#include "sam_record.hpp"

namespace weaver
{
//! Statistics for the haplotypes in the graph.
class HaplotypeStats
{
public:
  //!
  double static constexpr PHRED_TO_LOG{-0.23025850929940456840179914546843642076011014886287729760333279009675726};

  /*!
   * @brief Transforms read count to log
   *
   * @details
   * Calculated here:
   * https://www.wolframalpha.com/input?i2d=true&i=Divide%5B-1%2Clog10%5C%2840%29exp%5C%2840%291%5C%2841%29%5C%2841%29%5D
   */
  double static constexpr READ_COUNT_TO_LOG{25.0 * PHRED_TO_LOG}; // eps -28, value includes the -3 from the het case

  // -log10(0.5) * PHRED_TO_LOG
  double static constexpr HET_COUNT_TO_LOG{3.01029995663981195213738894724493026768189881462108541 * PHRED_TO_LOG};

  //! Pointer to the associated haplotypes.
  MMI::T_haplotypes const * haps_ptr{nullptr};

  //! All haplotype edit statistics are recorded here for each snid
  std::vector<std::vector<EditStats>> snid_haps_stats{};

  //! Set to true if the haplotype stats should be used/it is running the ambigous reads
  bool is_using_stats{false};

  void set_haps_ptr(MMI::T_haplotypes const * new_haps_ptr);

  //
  void run_hap_weighting();

  //
  void run_hmm();

  explicit HaplotypeStats(const MMI::T_haplotypes * _haps);

  HaplotypeStats() = default;
  HaplotypeStats(HaplotypeStats const &) = default;
  HaplotypeStats(HaplotypeStats &&) = default;
  HaplotypeStats & operator=(HaplotypeStats const &) = default;
  HaplotypeStats & operator=(HaplotypeStats &&) = default;

  bool operator==(HaplotypeStats const & other_hap_stats) const;
};

void update_stats(HaplotypeStats & haplotype_stats, std::vector<int> const & edits, int const snid);

HaplotypeStats merge_haplotype_stats(std::vector<HaplotypeStats> && p_hap_stats);
void merge_and_then_scatter_haplotype_stats(std::vector<HaplotypeStats> & p_hap_stats);

} // namespace weaver
