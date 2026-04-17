#include "stable_contigs.hpp"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "gfa.hpp"
#include "hashmap.hpp"
#include "logging.hpp"

#include <gfa-priv.h>

namespace weaver
{
StableContigs stable_contigs; // global instance

StableContigs::StableContigs(int const _bucket_size) : bucket_size(_bucket_size)
{
}

void StableContigs::add_contig(std::string const & name, int const snid, int const min, int const max, int const rank)
{
  print_debug(_HERE_,
              " adding stable contig with name=",
              name,
              " snid=",
              snid,
              " min=",
              min,
              " max=",
              max,
              " rank=",
              rank,
              " sfa_idx=",
              contigs.size());

  if (name2sfa_idx.count(name) == 1)
  {
    print_error(_HERE_, " graph contains duplicated contig name=", name);
    std::exit(1);
    return;
  }

  int const sfa_idx{static_cast<int>(contigs.size())};
  name2sfa_idx[name] = sfa_idx;
  snid2sfa_idx[snid].push_back(sfa_idx);
  assert(max > 0);
  contigs.emplace_back(name, min, max, rank, get_num_buckets() + (max - 1) / bucket_size + 1);
}

void StableContigs::clear()
{
  contigs.clear();
  name2sfa_idx.clear();
  snid2sfa_idx.clear();
}

int StableContigs::get_bucket_index(uint64_t order) const
{
  if (order == std::numeric_limits<uint64_t>::max())
    return get_num_buckets();

  int const sfa_idx = static_cast<int>(order >> 32ull);
  int const pos = static_cast<uint32_t>(order); // get 32 least significant bits

#ifndef NDEBUG
  if (sfa_idx >= static_cast<int>(contigs.size()))
    print_warning(_HERE_, " sfa_idx pointing to invalid contig Unable to find in graph the sfa_idx=", sfa_idx);
#endif

  // Checks how many buckets have been created from all previous stable FASTA contigs
  int const prev_accumulated_num_buckets = (sfa_idx == 0) ? 0 : contigs[sfa_idx - 1].accumulated_num_buckets;
  int const bucket_index = prev_accumulated_num_buckets + pos / bucket_size;
  return bucket_index;
}

int StableContigs::get_num_buckets() const
{
  return contigs.empty() ? 0 : contigs[contigs.size() - 1].accumulated_num_buckets;
}

int StableContigs::get_num_contigs() const
{
  return contigs.size();
}

Contig const & StableContigs::get_contig_reference(std::string const & contig_name) const
{
  auto find_sfa_idx_it = name2sfa_idx.find(contig_name);

  if (find_sfa_idx_it == name2sfa_idx.end())
  {
    print_warning(_HERE_, " Unable to find stable contig name from name=", contig_name);
    return contigs[contigs.size() - 1];
  }
  else
  {
    return contigs[find_sfa_idx_it->second];
  }
}

int StableContigs::get_sfa_idx(int const snid, int const pos) const
{
  auto get_dist = [](Contig const & contig, int const pos) -> int
  {
    if (pos > contig.max)
      return pos - contig.max;
    else if (pos < contig.min)
      return contig.min - pos;
    else
      return 0;
  };

  auto find_sfa_idx_it = snid2sfa_idx.find(snid);

  if (find_sfa_idx_it == snid2sfa_idx.end())
  {
    print_warning(_HERE_, " Unable to find stable contig name from snid=", snid);
    return contigs.size() - 1;
  }

  // TODO binary search
  int const n_sfa = find_sfa_idx_it->second.size();

  for (int s{n_sfa - 1}; s > 0; --s)
  {
    int const c = find_sfa_idx_it->second[s];
    assert(c > 0);
    assert(c < static_cast<int>(contigs.size()));
    int const dist_this = get_dist(contigs[c], pos);

    if (dist_this == 0 || dist_this <= get_dist(contigs[c - 1], pos))
      return c;
  }

  return find_sfa_idx_it->second[0];
}

void set_stable_contigs_from_graph(GFA const & graph)
{
  stable_contigs.clear();
  gfa_t const & g = graph.get_graph();

  // Get all the stable contigs
  int n_sfa;
  int write_seq{0};
  gfa_sfa_t * r = gfa_gfa2sfa(&g, &n_sfa, write_seq); // from gfa-priv.h in gfatools

  for (int sfa_idx{0}; sfa_idx < n_sfa; ++sfa_idx)
  {
    gfa_sfa_t const & s = r[sfa_idx];
    std::string contig_name(g.sseq[s.snid].name);

    if (s.rank != 0)
      contig_name += "_" + std::to_string(s.soff) + "_" + std::to_string(s.soff + s.len);

    stable_contigs.add_contig(contig_name, s.snid, s.soff, s.soff + s.len, s.rank);
  }

  free(r);
}

} // namespace weaver
