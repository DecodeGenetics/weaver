#include "edit_stats.hpp"

#include <cassert>
#include <cmath>

#include "edit.hpp"
#include "haplotype_stats.hpp"
#include "logging.hpp"

namespace
{
double constexpr EPS{25.0};

} // namespace

namespace weaver
{
EditStats::EditStats(EditCalls const & edit_calls) :
  non_missing_haplotypes(1), edit_read_count(0), no_edit_read_count(0), count(0)
{
  // There is no "++non_missing_haplotypes" here because it is initialized with the reference allele added
  int const num_alt_haplotypes = edit_calls.calls.size(); // does not include the reference haplotype

  for (int c{0}; c < num_alt_haplotypes; ++c)
  {
    auto const call = edit_calls.calls[c];

    if (call != Edit::MISSING_CALL)
      ++non_missing_haplotypes;

    if (call == 1)
      ++count;
  }
}

bool EditStats::operator==(EditStats const & o) const
{
  return non_missing_haplotypes == o.non_missing_haplotypes && edit_read_count == o.edit_read_count &&
         no_edit_read_count == o.no_edit_read_count && count == o.count;
}

bool EditStats::operator!=(EditStats const & o) const
{
  return !(*this == o);
}

double EditStats::get_hom_no_edit_eps() const
{
  return no_edit_read_count * HaplotypeStats::HET_COUNT_TO_LOG - //
         edit_read_count * HaplotypeStats::READ_COUNT_TO_LOG;
}

double EditStats::get_hom_edit_eps() const
{
  return edit_read_count * HaplotypeStats::HET_COUNT_TO_LOG - //
         no_edit_read_count * HaplotypeStats::READ_COUNT_TO_LOG;
}

uint8_t EditStats::get_hom_call() const
{
  if (edit_read_count + no_edit_read_count <= 2)
    return Edit::MISSING_CALL;

  if (no_edit_read_count == 0)
    return 1;

  if (edit_read_count == 0)
    return 0;

  if ((no_edit_read_count * 20) < edit_read_count)
    return 1;
  else if ((edit_read_count * 20) < no_edit_read_count)
    return 0;

  // het calls go here
  return Edit::MISSING_CALL;
}

void EditStats::merge_with(EditStats const & other_edit_stats)
{
  assert(this->non_missing_haplotypes == other_edit_stats.non_missing_haplotypes);
  edit_read_count += other_edit_stats.edit_read_count;
  no_edit_read_count += other_edit_stats.no_edit_read_count;
  assert(this->count == other_edit_stats.count);
}

double EditStats::get_hom_edit_phred() const
{
  return edit_read_count * 10.0 * std::log10(2) + no_edit_read_count * EPS;
}

double EditStats::get_het_phred() const
{
  return (edit_read_count + no_edit_read_count) * 10.0 * std::log10(2);
}

double EditStats::get_hom_no_edit_phred() const
{
  return no_edit_read_count * 10.0 * std::log10(2) + edit_read_count * EPS;
}

double EditStats::get_edit_frequency() const
{
  assert(non_missing_haplotypes > 0);
  return static_cast<double>(count) / static_cast<double>(non_missing_haplotypes);
}

double EditStats::get_no_edit_frequency() const
{
  assert(non_missing_haplotypes > 0);
  return 1. - static_cast<double>(count) / static_cast<double>(non_missing_haplotypes);
}

bool constexpr IS_HAPLOID{true};

double EditStats::get_likelihood_of_edit() const
{
  double const e_freq = get_edit_frequency();
  return IS_HAPLOID ? e_freq : 2.0 * e_freq - e_freq * e_freq;
}

double EditStats::get_likelihood_of_no_edit() const
{
  double const no_e_freq = get_no_edit_frequency();
  return IS_HAPLOID ? no_e_freq : 2.0 * no_e_freq - no_e_freq * no_e_freq;
}

} // namespace weaver
