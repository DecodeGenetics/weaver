#include "region.hpp"

#include <algorithm>
#include <charconv>
#include <string>

#include "gfa.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "sketch_value.hpp"
#include "sr_seed.hpp"
#include "stable_contigs.hpp"

namespace weaver
{
Region::Region(std::string const & region)
{
  auto find_colon_it = std::find(region.cbegin(), region.cend(), ':');
  this->chr = std::string(region.cbegin(), find_colon_it);

  if (find_colon_it == region.cend())
  {
    // chrN
    // then begin=0, end=max
  }
  else
  {
    auto find_dash_it = std::find(find_colon_it, region.cend(), '-');

    if (find_dash_it != region.cend())
    {
      // chrN:A-B
      auto ret = std::from_chars(&*std::next(find_colon_it), &*find_dash_it, this->begin);
      --this->begin;

      if (ret.ec != std::errc())
      {
        print_warning(_HERE_, " failed to parse region begin pos: ", region);
        assert(ret.ec == std::errc()); // will trigger
        this->begin = -1;
      }

      ret = std::from_chars(&*std::next(find_dash_it), &*region.cend(), this->end);

      if (ret.ec != std::errc())
      {
        print_warning(_HERE_, " failed to parse region end pos: ", region);
        assert(ret.ec == std::errc()); // will trigger
        this->end = -1;
      }
    }
    else
    {
      // chrN:A
      auto ret = std::from_chars(&*std::next(find_colon_it), &*find_dash_it, this->begin);
      --this->begin;

      if (ret.ec != std::errc())
      {
        print_warning(_HERE_, " failed to parse region pos: ", region);
        assert(ret.ec == std::errc()); // will trigger
        this->begin = -1;
      }

      this->end = this->begin + 1;
    }
  }
}

bool Region::check() const
{
  return chr.size() > 0 && end >= begin && begin >= 0 && end >= 0;
}

int Region::size() const
{
  return end - begin;
}

void parse_target_region_option(GFA const & /*gfa*/)
{
  Options & opts = *(Options::instance());

  if (opts.region.size() == 0)
    return;

  // Target region set
  Region target_region(opts.region);

  if (!target_region.check())
  {
    print_error(_HERE_, " Failed to parse region: '", opts.region, "'");
    std::exit(1);
  }

  auto find_it = stable_contigs.name2sfa_idx.find(target_region.chr);

  if (find_it == stable_contigs.name2sfa_idx.end())
  {
    print_error(_HERE_, " Did not find chromosome name in graph: '", target_region.chr, "'");
    std::exit(1);
  }

  assert(find_it->second >= 0);

  int const sfa_idx = find_it->second;
  assert(sfa_idx < static_cast<int>(stable_contigs.contigs.size()));
  Contig const & contig = stable_contigs.contigs[sfa_idx];
  target_region.begin = std::min(target_region.begin, contig.max);
  target_region.end = std::min(target_region.end, contig.max);

  // set region options
  opts.region_sfa_idx = sfa_idx;
  opts.region_rank = contig.rank;
  opts.region_lower_pos = target_region.begin;
  opts.region_upper_pos = target_region.end;
  // opts.region_lower_bucket_index = stable_contigs.get_bucket_index(order_begin);
  // opts.region_upper_bucket_index = stable_contigs.get_bucket_index(order_end);
  // assert(opts.region_upper_bucket_index >= opts.region_lower_bucket_index);
}

bool is_seed_in_region(GFA const & gfa, SRSeed const & seed)
{
  if (seed.is_empty())
    return false;

  Options const & copts = *(Options::const_instance());
  gfa_seg_t const & begin_segment = gfa.get_segment(sketch_value_rid(seed.begin_ref_value));
  int const begin_sfa_idx = stable_contigs.get_sfa_idx(begin_segment.snid, begin_segment.soff);

  if (/*begin_segment.rank != copts.region_rank ||*/ begin_sfa_idx != copts.region_sfa_idx)
    return false;

  gfa_seg_t const & end_segment = gfa.get_segment(sketch_value_rid(seed.end_ref_value));
  int const end_sfa_idx = stable_contigs.get_sfa_idx(end_segment.snid, end_segment.soff);

  if (/*end_segment.rank != copts.region_rank ||*/ end_sfa_idx != copts.region_sfa_idx)
    return false;

  Contig const & begin_contig = stable_contigs.contigs[begin_sfa_idx];
  Contig const & end_contig = stable_contigs.contigs[end_sfa_idx];
  int const begin_pos = sketch_value_pos(seed.begin_ref_value) + begin_segment.soff - begin_contig.min;
  int const end_pos = sketch_value_pos(seed.end_ref_value) + end_segment.soff - end_contig.min;

  return begin_pos <= end_pos ? begin_pos <= copts.region_upper_pos && end_pos >= copts.region_lower_pos
                              : end_pos <= copts.region_upper_pos && begin_pos >= copts.region_lower_pos;
}

bool are_all_seeds_outside_region(GFA const & gfa, std::vector<SRSeed> const & seeds)
{
  for (SRSeed const & seed : seeds)
  {
    if (is_seed_in_region(gfa, seed))
      return false;
  }

  return true;
}

} // namespace weaver
