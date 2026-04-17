#include "make_mmi.hpp"

#include <cstdint>
#include <fstream>
#include <mutex>
#include <string>
#include <vector>

#include <paw/station.hpp>

#include "edit.hpp"
#include "gfa.hpp"     // GFA
#include "gfa_arc.hpp" // arc_begin
#include "hashmap.hpp"
#include "index_io.hpp"
#include "io.hpp" // hts_file_ptr
#include "logging.hpp"
#include "options.hpp"
#include "segment.hpp"
#include "seq_view.hpp"
#include "sequence_utils.hpp"
#include "sketch.hpp"
#include "sketch_cache.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp"
#include "stable_contigs.hpp"
#include "variant.hpp"

namespace
{
void read_alts_file(phmap::flat_hash_set<std::string> & alts, std::string const & graph_fn)
{
  using namespace weaver;

  std::string alts_fn = graph_fn + ".alts";

  if (!weaver::filesystem::exists(alts_fn))
    return;

  std::ifstream alts_in;
  alts_in.open(alts_fn.c_str());

  if (!alts_in.is_open())
  {
    print_warning(_HERE_, " Unable to open alt segment file '", alts_fn, "'");
    return;
  }

  std::string segment_name;

  while (getline(alts_in, segment_name))
  {
    if (not segment_name.empty())
      alts.insert(segment_name);
  }

  print_info("Read ", alts.size(), " alt segments from file '", alts_fn, "'");
}

std::mutex mmi_uniq_mutex; // protects mmi_uniq

std::vector<uint8_t> get_edit_calls(std::vector<uint16_t> const & hap_calls, uint16_t allele)
{
  auto const num_hap_calls = hap_calls.size();
  std::vector<uint8_t> edit_calls(hap_calls.size(), 0);

  assert(allele < num_hap_calls);

  for (size_t c{0}; c < num_hap_calls; ++c)
  {
    auto const hap_call = hap_calls[c];

    if (hap_call == allele)
    {
      edit_calls[c] = 1;
    }
    else if (hap_call == weaver::Variant::MISSING_CALL)
    {
      edit_calls[c] = weaver::Edit::MISSING_CALL;
    }
#ifndef NDEBUG
    else
    {
      assert(edit_calls[c] == 0);
    }
#endif

    // else do nothing, call is already 0
  }

  return edit_calls;
}

void update_edit_calls(std::vector<uint8_t> & edit_calls, std::vector<uint16_t> const & hap_calls, uint16_t allele)
{
  auto const num_hap_calls = hap_calls.size();
  assert(allele < num_hap_calls);

  for (size_t c{0}; c < num_hap_calls; ++c)
  {
    if (hap_calls[c] == allele)
      edit_calls[c] = 1;
  }
}

/*
void add_snp_edits(phmap::flat_hash_set<weaver::Edit, weaver::EditHash> & edits,
                   int pos,
                   std::string_view ref,
                   std::string_view alt)
{
  using namespace weaver;
  assert(ref.size() == alt.size());

  for (size_t p{0}; p < ref.size(); ++p)
  {
    if (ref[p] != alt[p])
      edits.emplace(pos + p, 'X', alt[p]);
  }
}
*/

void add_snp_edits(phmap::flat_hash_map<weaver::Edit, std::vector<uint8_t>, weaver::EditHash> & edits,
                   weaver::Variant const & var,
                   uint16_t a)
{
  using namespace weaver;
  assert(var.seqs.size() >= 2);
  assert(a > 0);
  std::vector<char> const & ref = var.seqs[0];
  std::vector<char> const & alt = var.seqs[a];

  assert(ref.size() == alt.size());

  for (size_t p{0}; p < ref.size(); ++p)
  {
    if (ref[p] != alt[p])
    {
      // edits use 0-based indexing but variants 1-based, so we need to do "- 1"
      Edit new_edit(var.pos + p - 1, 'X', alt[p]);
      auto it = edits.find(new_edit);

      if (it == edits.end())
        edits[new_edit] = get_edit_calls(var.calls, a); // insert new
      else
        update_edit_calls(/*edit_calls=*/it->second, /*hap_calls=*/var.calls, a); // update existing values
    }
  }
}

void add_indel_edits(phmap::flat_hash_map<weaver::Edit, std::vector<uint8_t>, weaver::EditHash> & edits,
                     weaver::Variant const & var,
                     uint16_t a)
{
  using namespace weaver;

  // TODO test before enabling
  // return;

  assert(var.seqs.size() >= 2);
  assert(a > 0);
  assert(var.seqs[0].size() > 0);
  assert(var.seqs[a].size() > 0);
  assert(var.seqs[0].size() != var.seqs[a].size());

  // expect the padded base
  if (var.seqs[0][0] == var.seqs[a][0])
  {
    // int const prefix_size = get_prefix_size(var.seqs[0], var.seqs[a]);
    int const suffix_size = get_suffix_size(var.seqs[0], var.seqs[a]);

    print_debug(_HERE_,
                " indel ",
                std::string_view(var.seqs[0].data(), var.seqs[0].size()),
                " ",
                std::string_view(var.seqs[a].data(), var.seqs[a].size()),
                " suffix size = ",
                suffix_size);

    if (static_cast<int>(var.seqs[0].size()) > static_cast<int>(var.seqs[a].size()) &&
        static_cast<int>(var.seqs[a].size()) == suffix_size + 1)
    {
      // simple deletion
      Edit new_edit(var.pos,
                    'D',
                    std::vector<char>(std::next(var.seqs[0].begin(), 1), std::next(var.seqs[0].end(), -suffix_size)));

      auto it = edits.find(new_edit);

      if (it == edits.end())
        edits[new_edit] = get_edit_calls(var.calls, a); // insert new
      else
        update_edit_calls(/*edit_calls=*/it->second, /*hap_calls=*/var.calls, a); // update existing values
    }
    else if (static_cast<int>(var.seqs[a].size()) > static_cast<int>(var.seqs[0].size()) &&
             static_cast<int>(var.seqs[0].size()) == suffix_size + 1)
    {
      // simple insertion
      Edit new_edit(var.pos,
                    'I',
                    std::vector<char>(std::next(var.seqs[a].begin(), 1), std::next(var.seqs[a].end(), -suffix_size)));

      auto it = edits.find(new_edit);

      if (it == edits.end())
        edits[new_edit] = get_edit_calls(var.calls, a); // insert new
      else
        update_edit_calls(/*edit_calls=*/it->second, /*hap_calls=*/var.calls, a); // update existing values
    }
  }

  // // TODO complex variant
  // std::string_view ref(var.seqs[0].data(), var.seqs[0].size());
  // std::string_view alt(var.seqs[a].data(), var.seqs[a].size());
}

std::vector<weaver::EditCalls> get_edits(weaver::GFA const & /*gfa*/, std::vector<weaver::Variant> const & hap_variants)
{
  using namespace weaver;

  phmap::flat_hash_map<Edit, std::vector<uint8_t>, EditHash> edits;
  int const num_hap_variants = static_cast<int>(hap_variants.size());
  // int const num_calls = num_hap_variants == 0 ? 0 : hap_variants[0].calls.size();

  for (int v{0}; v < num_hap_variants; ++v)
  {
    Variant const & var = hap_variants[v];
    assert(var.seqs.size() >= 2);
    std::vector<char> const & ref = var.seqs[0];

    for (int a{1}; a < static_cast<int>(var.seqs.size()); ++a)
    {
      std::vector<char> const & alt = var.seqs[a];

      if (ref.size() == alt.size())
        add_snp_edits(edits, var, a);
      else
        add_indel_edits(edits, var, a);
    }
  }

  int const num_edits{static_cast<int>(edits.size())};
  std::vector<EditCalls> contig_edit_calls;
  contig_edit_calls.reserve(num_edits);

  Options const & copts = *(Options::const_instance());
  int const N_DIV = copts.min_edit_as_effect < 1e-10
                    ? std::numeric_limits<int>::max()
                    : std::round(static_cast<double>(copts.match + copts.mismatch) / copts.min_edit_as_effect);

  // store in a vector (which will be sorted later)
  for (auto it = edits.begin(); it != edits.end(); ++it)
  {
    std::vector<uint8_t> const & calls = it->second;
    int n{static_cast<int>(calls.size())};

    if (n >= N_DIV)
    {
      // for graphs with a lot of haplotypes, remove variants which are so rare that they cannot practically affect
      // results to reduce memory/cpu time
      int const MIN_ALT_COUNT{n / N_DIV + 1};
      int c{0};
      int alt_count{0};

      while (alt_count < MIN_ALT_COUNT && c < n)
      {
        if (calls[c] == 1)
          ++alt_count;

        ++c;
      }

      if (alt_count >= MIN_ALT_COUNT)
        contig_edit_calls.emplace_back(it->first.pos, it->first.type, it->first.seq, std::move(it->second));
    }
    else
    {
      contig_edit_calls.emplace_back(it->first.pos, it->first.type, it->first.seq, std::move(it->second));
    }
  }

  std::sort(contig_edit_calls.begin(), contig_edit_calls.end());
  assert(std::unique(contig_edit_calls.begin(), contig_edit_calls.end()) == contig_edit_calls.end());
  return contig_edit_calls;
}

void sort_hits_and_remove_keys_with_too_many_hits(weaver::MMI::T_multi_value_map & mmi_multi_map,
                                                  weaver::MMI::T_uniq_value_map & mmi_uniq_map,
                                                  weaver::GFA const & gfa)
{
  size_t constexpr MAX_MINIMIZER_HITS{262144}; // 2^18
  std::vector<uint64_t> keys_to_remove;        // stores keys with extreme amount of hits
  auto order_mmi = [&gfa](uint64_t const a, uint64_t const b) -> bool { return gfa.mmi_value_order(a, b); };

  for (auto & mmi_pair : mmi_multi_map)
  {
    if (mmi_pair.second.size() > MAX_MINIMIZER_HITS) // 2^18
      keys_to_remove.push_back(mmi_pair.first);
    else
      std::sort(mmi_pair.second.begin(), mmi_pair.second.end(), order_mmi);
  }

#ifndef NDEBUG
  auto compare_approx_value = [&gfa](uint64_t a, uint64_t b) { //
    return gfa.get_approximate_stable_position(a) < gfa.get_approximate_stable_position(b);
  };

  for (auto const & mmi_pair : mmi_multi_map)
  {
    auto const & v = mmi_pair.second; // references the vector containing all values for this minimizer key
    assert(std::find(v.begin(), v.end(), std::numeric_limits<uint64_t>::max()) == v.end());
    assert(std::is_sorted(v.begin(), v.end(), compare_approx_value));
  }
#endif // NDEBUG

  if (keys_to_remove.size() > 0)
  {
    print_info("Removing ", keys_to_remove.size(), " keys with > ", MAX_MINIMIZER_HITS, " hits.");

    // remove keys with extreme amount of hits (see above)
    for (auto key_to_remove : keys_to_remove)
    {
      auto find_erase_it = mmi_multi_map.find(key_to_remove);
      assert(find_erase_it != mmi_multi_map.end());
      auto find_erase_uniq_it = mmi_uniq_map.find(key_to_remove);
      assert(find_erase_uniq_it != mmi_uniq_map.end());

      if (find_erase_it == mmi_multi_map.end()) // should never happen but is very cheap to check for
        continue;

      if (find_erase_uniq_it == mmi_uniq_map.end()) // should never happen but is very cheap to check for
        continue;

      mmi_multi_map.erase(find_erase_it);
      mmi_uniq_map.erase(find_erase_uniq_it);
    }
  }
}

void sketch_rid(std::vector<std::pair<uint64_t, uint64_t>> & new_sketches, //
                weaver::GFA const & gfa,
                uint32_t const rid,
                int const k,
                int const w)
{
  using namespace weaver;

  using u64_pair = std::pair<uint64_t, uint64_t>;
  std::string_view const segment_sequence = gfa.view_segment(rid);
  SketchCache<u64_pair> cache{}; // cache for sketches

  sketch(new_sketches,
         cache,
         segment_sequence,
         w,
         k,
         rid,
         /*is_reading_forward=*/true,
         /*max_read_size=*/-1,
         /*read_init=*/0,
         /*add_pos=*/0);

  uint32_t const fwd_v = forward_segment_id2vertex_id(rid); // forward vertex id
  uint32_t const rev_v = fwd_v + 1;                         // reverse vertex id

  gfa_arc_t const * rev_arc_ptr = arc_begin(gfa, rev_v);
  gfa_arc_t const * rev_end_ptr = arc_end(gfa, rev_v, rev_arc_ptr); // arc pointer behind the reverse vertex id

  /* Check arcs from a reverse oriented segment
   *   s1 - s2 +  (or s2<s1) In this case add iff s1 rid is smaller than s2 rid.
   * We should not add s1 - s2 - since we will add those the other way around.
   */
  for (/*no initial*/; rev_arc_ptr != rev_end_ptr; ++rev_arc_ptr)
  {
    uint32_t const next_rid = static_cast<uint32_t>(rev_arc_ptr->w) >> 1; // rid of segment on the other side of the arc
    bool const is_reading_forward{(rev_arc_ptr->w & 1) == 0u}; // true iff we are reading from the beginning of s2

    if (!is_reading_forward || next_rid < rid) // such that links are only crossed once
      continue;

    print_debug(_HERE_, " sketching arc ", arc_to_string(*rev_arc_ptr));
    gfa_seg_t const & segment = gfa.get_segment(rid);
    gfa_seg_t const & next_segment = gfa.get_segment(next_rid);
    std::string segment_prefix = get_segment_prefix<std::string>(segment, k + 3 * w);
    std::string next_segment_prefix = get_segment_prefix<std::string>(next_segment, k + 3 * w);
    int const next_segment_prefix_size{static_cast<int>(next_segment_prefix.size())};
    std::string next_segment_seq_rev = get_reverse_complement(next_segment_prefix);

    print_debug(_HERE_, " next_segment_str_rev=", next_segment_seq_rev, " (length=", next_segment_prefix_size, ")");

    next_segment_seq_rev.append(segment_prefix);
    print_debug(_HERE_, " after append ", next_segment_seq_rev, " (length=", next_segment_seq_rev.size(), ")");
    int first_pos{-1};

    if (new_sketches.empty())
    {
      std::vector<u64_pair> new_sketches_cp; // new sketches for the copy
      SketchCache cache_cp(cache);           // don't modify cache, so we need to copy
      push_last_sketch(new_sketches_cp, cache_cp, -1, k);

      if (new_sketches_cp.empty())
      {
        print_debug(_HERE_, " No sketches found in ", rid);
        continue; // if still nothing then give up
      }

      first_pos = sketch_value_pos(new_sketches_cp[0].second);
    }
    else
    {
      first_pos = sketch_value_pos(new_sketches[0].second);
    }

    print_debug(_HERE_, " first_pos=", first_pos);

    // time for sketching
    if (first_pos >= 0)
    {
      SketchCache<u64_pair> empty;
      std::vector<u64_pair> potential_sketches;

      // sketch with the previous cache
      sketch(potential_sketches,
             empty,
             next_segment_seq_rev,
             w,
             k,
             rid,
             /*is_reading_forward=*/true,
             /*max_read_size=*/-1,
             /*read_init=*/0,
             /*add_pos=*/0);
      push_last_sketch(potential_sketches, empty, -1, k);

      for (auto & sketch : potential_sketches)
      {
        int const sketch_pos = sketch_value_pos(sketch.second);

        // Add sketches that have a position greater or equal to
        // "next_segment_prefix_size" (usually k+3w)
        // and strictly smaller than the first sketch of the rid
        if ((sketch_pos + k / 2 + w) >= next_segment_prefix_size)
        {
          if ((sketch_pos - next_segment_prefix_size) >= first_pos)
            break; // we can break since sketch_pos will only grow bigger

          // reduce pos to account for sequence from previous segment
          if (sketch_pos < next_segment_prefix_size)
          {
            // sketch is on the other RID
            uint64_t const new_pos = next_segment_prefix_size - sketch_pos - 1;

            sketch.second = (static_cast<uint64_t>(next_rid) << 32) |    // new RID
                            (new_pos << 1) |                             // new pos
                            static_cast<uint64_t>(!(sketch.second & 1)); // flipped strand
          }
          else
          {
            sketch.second -= static_cast<uint64_t>(next_segment_prefix_size << 1); // sketch is on this RID
          }

          print_debug(_HERE_, " adding sketch ", sketch_to_string(sketch));
          assert(sketch_value_rid(sketch.second) < gfa.get_num_segments());
          new_sketches.push_back(sketch); // add sketch
        }
      }
    }
  }

  gfa_arc_t const * fwd_arc_ptr = arc_begin(gfa, fwd_v);
  gfa_arc_t const * fwd_end_ptr = arc_end(gfa, fwd_v, fwd_arc_ptr); // arc pointer behind the forward vertex id

  /* Check arcs from a forward oriented segment
   *   s1 + s2 + Add all such arcs.
   * or
   *   s1 + s2 - with s1<s2. In this case add arc iff s1 rid is smaller than s2 rid.
   */
  for (/*no initial*/; fwd_arc_ptr != fwd_end_ptr; ++fwd_arc_ptr)
  {
    // Sketch the sequence across an arc
    uint32_t const next_rid = static_cast<uint32_t>(fwd_arc_ptr->w) >> 1; // rid of segment on the other side of the arc
    bool const is_next_forward{(fwd_arc_ptr->w & 1) == 0u};               // true if s1>s2, false otherwise

    if (next_rid < rid && !is_next_forward)
      continue; // such that links are only crossed once

    // Sketch the next segment sequence to find the position of the first sketch.
    // print_debug(_HERE_, " sketching arc ", arc_to_string(*fwd_arc_ptr));
    std::string_view const next_sequence = gfa.view_segment(next_rid); // sequence of next rid
    int const max_read_size = get_position_of_the_first_sketch(next_sequence, next_rid, w, k, is_next_forward);
    // print_debug(_HERE_, " Position of the first sketch is ", max_read_size);

    // Resketch the next segment sequence bases using the cache from previous vertex, up to max_read_size base-pairs.
    SketchCache temp_cache(cache); // copy the cache for each arc

    sketch(new_sketches,
           temp_cache,
           next_sequence,
           w,
           k,
           next_rid,
           /*is_reading_forward=*/is_next_forward,
           /*max_read_size=*/max_read_size,
           /*read_init=*/0,
           /*add_pos=*/0);

    push_last_sketch(new_sketches, temp_cache, max_read_size, k);
  }

  push_last_sketch(new_sketches, cache, -1, k); // push the very last sketch of the current segment
}

void populate_uniq_and_multi_mmi_maps(weaver::MMI::T_uniq_value_map & uniq_map,
                                      weaver::MMI::T_multi_value_map & multi_map,
                                      weaver::MMI::T_unique_set const & mmi_uniq)
{
  using namespace weaver;
  uint64_t constexpr MISSING_VALUE{std::numeric_limits<uint64_t>::max()};

  for (auto sketch : mmi_uniq)
  {
    auto find_uniq_it = uniq_map.find(sketch.first);

    if (find_uniq_it == uniq_map.end())
    {
      // never before seen sketch key
      uniq_map[sketch.first] = sketch.second;
    }
    else if (find_uniq_it->second == MISSING_VALUE)
    {
      // found a previous value which is missing
      assert(sketch.second != std::numeric_limits<uint64_t>::max());
      assert(multi_map.find(sketch.first) != multi_map.end());
      multi_map[sketch.first].push_back(sketch.second);
    }
    else
    {
      // found a previous value which was not missing
      assert(find_uniq_it->second != MISSING_VALUE);
      assert(sketch.second != MISSING_VALUE);
      assert(multi_map.find(sketch.first) == multi_map.end()); // first and only time we see it

      // make a new entry in multi_map with both previous value and new value
      std::vector<uint64_t> hits{find_uniq_it->second, sketch.second};

      assert(hits.size() == 2);
      assert(hits[0] != MISSING_VALUE);
      assert(hits[1] != MISSING_VALUE);

      find_uniq_it->second = MISSING_VALUE;
      multi_map[sketch.first] = hits; // add entry
    }
  }
}

bool is_sketch_with_unique_context(std::vector<std::pair<uint64_t, uint64_t>>::const_iterator const & new_sketch_begin,
                                   std::vector<std::pair<uint64_t, uint64_t>>::const_iterator const & new_sketch_end,
                                   std::vector<std::pair<uint64_t, uint64_t>>::const_iterator const & new_sketch_it,
                                   weaver::MMI::T_uniq_value_map const & pre_uniq_map)
{
  if (pre_uniq_map.count(new_sketch_it->first) == 0)
    return true; // if the sketch itself is unique then

  auto const rid = weaver::sketch_value_rid(new_sketch_it->second);
  int constexpr SIZE_OF_CONTEXT{3};
  bool is_fwd_unique{false};

  {
    auto new_sketch_fwd = std::next(new_sketch_it);

    for (int i{0}; i < SIZE_OF_CONTEXT; ++i)
    {
      // rid change means context is no longer relevant, we assume the context is unique in this case
      if (new_sketch_fwd != new_sketch_end || weaver::sketch_value_rid(new_sketch_fwd->second) != rid ||
          pre_uniq_map.count(new_sketch_fwd->first) == 0)
      {
        is_fwd_unique = true;
        break;
      }

      ++new_sketch_fwd;
    }
  }

  if (is_fwd_unique)
  {
    auto new_sketch_rev = new_sketch_it;

    for (int i{0}; i < SIZE_OF_CONTEXT && new_sketch_rev != new_sketch_begin; ++i)
    {
      if (new_sketch_it == new_sketch_begin)
        return true; // same as change of rid

      --new_sketch_rev;

      // rid change means context is no longer relevant, we assume the context is unique in this case
      if (weaver::sketch_value_rid(new_sketch_rev->second) != rid || pre_uniq_map.count(new_sketch_rev->first) == 0)
        return true;
    }
  }

  return false;
}

void get_minimizers_of_a_haplotype_in_variant_data(int const /*haplotype index=*/c,
                                                   weaver::MMI::T_unique_set * mmi_uniq_ptr,
                                                   weaver::MMI::T_uniq_value_map const * const pre_uniq_map_ptr,
                                                   weaver::GFA const * const gfa_ptr,
                                                   std::vector<weaver::Variant> const * const small_variants_ptr,
                                                   int const rid,
                                                   int const k,
                                                   int const w)
{
  using namespace weaver;
  using u64_pair = std::pair<uint64_t, uint64_t>;

  assert(mmi_uniq_ptr != nullptr);
  assert(pre_uniq_map_ptr != nullptr);
  assert(gfa_ptr != nullptr);
  assert(small_variants_ptr != nullptr);

  GFA const & gfa = *gfa_ptr;
  std::vector<Variant> const & small_variants = *small_variants_ptr;

  std::vector<u64_pair> new_sketches;
  weaver::SketchCache<u64_pair> cache{}; // cache for sketches
  int const num_sites{static_cast<int>(small_variants.size())};
  assert(num_sites > 0);
  // std::vector<SeqView> seq_views;

  int curr_rid = rid;
  int curr_rid_pos{0};

  for (int s{0}; s <= num_sites; ++s)
  {
    // loop over all small variant sites
    // the final iteration has s==num_sites to add the remaining sequence after the last variant
    gfa_seg_t const * segment_ptr = &gfa.get_segment(curr_rid);
    assert(segment_ptr->rank == 0);

    // variant site is nullptr in the last iteration only
    weaver::Variant const * site = s == num_sites ? nullptr : &small_variants[s];

    if (site != nullptr)
    {
      // print_debug(_HERE_, " s=", s, " seqs.size()=", site->seqs.size(), " site.pos=", site->pos, " rid=",
      // curr_rid);
      assert(c < static_cast<int>(site->calls.size()));
      auto const call = site->calls[c];

      // reference or a missing call, ignore
      if (call == 0 || call == weaver::Variant::MISSING_CALL)
        continue;
    }

    if (site == nullptr || (site->pos - 1 - segment_ptr->soff) >= segment_ptr->len)
    {
      // print_debug(_HERE_, " site not on segment. rid pos=", curr_rid_pos, " segment_ptr->len=",
      // segment_ptr->len);
      assert(curr_rid_pos <= segment_ptr->len);

      std::string_view const rid_segment_seq(segment_ptr->seq + curr_rid_pos, segment_ptr->len - curr_rid_pos);
      // seq_views.emplace_back(rid_segment_seq, curr_rid, curr_rid_pos);

      // the next variant is not on this segment
      sketch_boundary(new_sketches,
                      cache,
                      rid_segment_seq,
                      w,
                      k,
                      curr_rid,
                      /*is_reading_forward=*/true,
                      /*max_read_size=*/-1,
                      /*read_init=*/0,
                      /*add_pos=*/curr_rid_pos);

      // go to the next rid
      gfa_arc_t const * arc = gfa.get_arc_to_next_stable_sequence(forward_segment_id2vertex_id(curr_rid));

      if (arc == nullptr)
      {
        push_last_sketch(new_sketches, cache, -1, k);
        break;
      }

      // set a new segment
      curr_rid = (arc->w >> 1);
      curr_rid_pos = 0;
      --s; // no site was processed
      continue;
    }
    else
    {
      assert(site != nullptr);
      auto const call = site->calls[c];
      assert(call != Variant::MISSING_CALL);
      assert(call < site->seqs.size());
      assert(call != 0);
      std::vector<char> const & ref_seq = site->seqs[0];
      std::vector<char> const & call_seq = site->seqs[call]; // should not be the same as ref_seq
      assert(call_seq != ref_seq);
      int const segment_site_pos = site->pos - 1 - segment_ptr->soff; // -1 to switch to 0-based indexing
      assert(segment_site_pos >= curr_rid_pos);

      if (segment_site_pos + static_cast<int>(ref_seq.size()) >= segment_ptr->len)
      {
        print_debug(_HERE_,
                    " incompatible call @ ",
                    site->pos,
                    " ",
                    std::string_view(call_seq.data(), call_seq.size()),
                    " segment_pos=",
                    segment_site_pos,
                    " segment_len=",
                    segment_ptr->len);

        // the next variant is not on this segment
        std::string_view const rid_segment_seq(segment_ptr->seq + curr_rid_pos, segment_ptr->len - curr_rid_pos);
        // seq_views.emplace_back(rid_segment_seq, curr_rid, curr_rid_pos);

        sketch_boundary(new_sketches,
                        cache,
                        rid_segment_seq,
                        w,
                        k,
                        curr_rid,
                        /*is_reading_forward=*/true,
                        /*max_read_size=*/-1,
                        /*read_init=*/0,
                        /*add_pos=*/curr_rid_pos);

        // go to the next rid
        gfa_arc_t const * arc = gfa.get_arc_to_next_stable_sequence(forward_segment_id2vertex_id(curr_rid));

        if (arc == nullptr)
        {
          push_last_sketch(new_sketches, cache, -1, k);
          break;
        }

        // set a new segment
        curr_rid = (arc->w >> 1);
        curr_rid_pos = 0;
      }
      else
      {
        // the next variant is on this segment
        if (segment_site_pos > curr_rid_pos)
        {
          std::string_view const rid_segment_seq(segment_ptr->seq + curr_rid_pos, segment_site_pos - curr_rid_pos);
          // seq_views.emplace_back(rid_segment_seq, curr_rid, curr_rid_pos);

          sketch_boundary(new_sketches,
                          cache,
                          rid_segment_seq,
                          w,
                          k,
                          curr_rid,
                          /*is_reading_forward=*/true,
                          /*max_read_size=*/-1,
                          /*read_init=*/0,
                          /*add_pos=*/curr_rid_pos);

          curr_rid_pos = segment_site_pos;
        }

        std::string_view const ref_view(ref_seq.data(), ref_seq.size());
        std::string_view const call_view(call_seq.data(), call_seq.size());
        // seq_views.emplace_back(ref_view, call_view, curr_rid, curr_rid_pos);

        if (call_seq.size() <= ref_seq.size())
        {
          // deletion or substitution call
          sketch(new_sketches,
                 cache,
                 call_view,
                 w,
                 k,
                 curr_rid,
                 /*is_reading_forward=*/true,
                 /*max_read_size=*/-1,
                 /*read_init=*/0,
                 /*add_pos=*/curr_rid_pos);
        }
        else // call_seq.size() > ref_seq.size()
        {
          // insertion call
          std::vector<u64_pair> tmp_sketches;

          sketch(tmp_sketches,
                 cache,
                 call_view,
                 w,
                 k,
                 curr_rid,
                 /*is_reading_forward=*/true,
                 /*max_read_size=*/-1,
                 /*read_init=*/0,
                 /*add_pos=*/curr_rid_pos + (call_seq.size() - ref_seq.size() + 1) / 2);

          for (u64_pair const & tmp_sketch : tmp_sketches)
          {
            if (sketch_value_pos(tmp_sketch.second) <= curr_rid_pos)
              new_sketches.push_back(tmp_sketch);
            else
              break;
          }

          if (cache.min.second != std::numeric_limits<uint64_t>::max() && //
              sketch_value_rid(cache.min.second) == curr_rid)
          {
            // fix position of previous minimizer, if it exists
            int const insert_size = call_seq.size() - ref_seq.size();
            assert(insert_size > 0);

            if (sketch_value_pos(cache.min.second) > insert_size)
              cache.min.second -= (insert_size << 1);
          }
        }

        curr_rid_pos += ref_seq.size();
        assert(curr_rid_pos < segment_ptr->len);
      }
    }
  }

  // Insert new sketches, if any
  if (new_sketches.size() > 0)
  {
    std::lock_guard<std::mutex> lock(mmi_uniq_mutex); // insert sketches into mmi_uniq with thread safety
    auto const new_sketch_begin = new_sketches.cbegin();
    auto const new_sketch_end = new_sketches.cend();
    MMI::T_unique_set & mmi_uniq = *mmi_uniq_ptr;
    MMI::T_uniq_value_map const & pre_uniq_map = *pre_uniq_map_ptr;
    int num_invalid_sketches{0};

    for (auto new_sketch_it = new_sketch_begin; new_sketch_it != new_sketch_end; ++new_sketch_it)
    {
      if (is_sketch_value_valid(gfa, new_sketch_it->second))
      {
        if (mmi_uniq.count(*new_sketch_it) == 0 &&
            is_sketch_with_unique_context(new_sketch_begin, new_sketch_end, new_sketch_it, pre_uniq_map))
        {
          mmi_uniq.insert(mmi_uniq.end(), *new_sketch_it);
        }
      }
      else
      {
        ++num_invalid_sketches;
      }
    }

    if (num_invalid_sketches > 0)
    {
      print_debug(_HERE_, " num invalid sketches in rid = ", rid, " were ", num_invalid_sketches);
      num_invalid_sketches = 0;
    }
  }
}

void get_minimizer_of_variant_data(weaver::MMI::T_unique_set & mmi_uniq,
                                   weaver::MMI::T_uniq_value_map const & pre_uniq_map,
                                   weaver::MMI::T_haplotypes & haplotypes,
                                   weaver::GFA const & gfa,
                                   std::string const & vcf_fn,
                                   int k,
                                   int w)
{
  using namespace weaver;

  int const num_segments{gfa.get_num_segments()};

  // Input streams
  weaver::hts_file_ptr in_vcf(nullptr, weaver::close_hts_file);
  weaver::tbx_t_ptr in_tbx(nullptr, weaver::close_tbx_t);
  weaver::hts_itr_t_ptr in_it(nullptr, weaver::close_hts_itr_t);

  // go through the base again and add small variants
  print_info(_HERE_, " Reading small variants from ", vcf_fn);

  in_vcf = weaver::open_hts_file(vcf_fn.c_str(), "r"); // open vcf.gz
  in_tbx = weaver::open_tbx_t(vcf_fn.c_str());         // open vcf.gz.tbi

  for (int rid{0}; rid < num_segments; ++rid)
  {
    gfa_seg_t const & first_segment = gfa.get_segment(rid);

    if (first_segment.rank != 0)
      continue;

    // Get stable sequence name/contig name
    gfa_sseq_t const & stable_sequence = gfa.get_stable_sequence(first_segment);

    if (first_segment.soff != stable_sequence.min)
      continue;

    std::string const contig_name(stable_sequence.name);

    print_debug(_HERE_, " Processing rid=", rid, " contig name=", contig_name, " with snid=", first_segment.snid);
    in_it = weaver::open_hts_itr_t(in_tbx.get(), contig_name.c_str(), 1, std::numeric_limits<int>::max());

    assert(first_segment.snid == static_cast<int>(haplotypes.size()));

    if (first_segment.snid >= static_cast<int>(haplotypes.size()))
      haplotypes.resize(first_segment.snid + 1);

    assert(first_segment.snid < static_cast<int>(haplotypes.size()));
    std::vector<Variant> small_variants = get_variants_in_a_region(in_vcf, in_tbx, in_it);

    print_debug(_HERE_, " parsed ", small_variants.size(), " small variants.");

#ifndef NDEBUG
    if (small_variants.size() < 500)
    {
      for (Variant const & small_variant : small_variants)
        print_debug(_HERE_, " ", small_variant.to_string());
    }
#endif // NDEBUG

    haplotypes[first_segment.snid] = get_edits(gfa, small_variants);

    assert(std::is_sorted(small_variants.begin(), small_variants.end()));
    int const num_calls = small_variants.size() == 0 ? 0 : small_variants[0].calls.size();

    Options const & copts = *(Options::const_instance());
    int const threads = copts.threads <= 1 ? 1 : copts.threads;
    paw::Station station(threads);

    for (int c{0}; c < num_calls; ++c)
    {
      if (c == (num_calls - 1))
      {
        station.add_to_thread(static_cast<std::size_t>(threads - 1),
                              get_minimizers_of_a_haplotype_in_variant_data,
                              c,
                              &mmi_uniq,
                              &pre_uniq_map,
                              &gfa,
                              &small_variants,
                              rid,
                              k,
                              w);
      }
      else
      {
        station.add(get_minimizers_of_a_haplotype_in_variant_data,
                    c,
                    &mmi_uniq,
                    &pre_uniq_map,
                    &gfa,
                    &small_variants,
                    rid,
                    k,
                    w);
      }
    }

    {
      std::string thread_info = station.join();
      print_info("Alternative indexing chunks per thread: ", thread_info);
    }
  }
}

} // namespace

namespace weaver
{
MMI make_mmi_index(GFA const & gfa, std::string const & vcf_fn, int k, int w)
{
  if (k <= 0)
    k = K_DEFAULT;

  if (w <= 0)
    w = W_DEFAULT;

  using u64_pair = std::pair<uint64_t, uint64_t>;
  print_info("Indexing graph with k = ", k, " and w = ", w);
  MMI::T_unique_set mmi_uniq;
  int const num_segments{gfa.get_num_segments()};
  int64_t num_invalid_sketches{0};    // for debugging only
  bool are_any_non_zero_ranks{false}; // check if any sequence has a rank which is not 0

  // This loop sketches/finds minimizers on the reference sequence (sequence with rank=0)
  for (int rid{0}; rid < num_segments; ++rid)
  {
    gfa_seg_t const & segment = gfa.get_segment(rid);

    if (segment.rank != 0)
    {
      are_any_non_zero_ranks = true;
      continue; // we will sketch those later
    }

    std::vector<u64_pair> new_sketches{};
    sketch_rid(new_sketches, gfa, rid, k, w);

    for (auto const & new_sketch : new_sketches)
    {
      if (is_sketch_value_valid(gfa, new_sketch.second))
        mmi_uniq.insert(mmi_uniq.end(), new_sketch);
      else
        ++num_invalid_sketches;
    }

    if (num_invalid_sketches > 0)
    {
      print_debug(_HERE_, " num invalid sketches in rid = ", rid, " were ", num_invalid_sketches);
      num_invalid_sketches = 0;
    }
  }

  print_debug(_HERE_, " are_any_non_zero_ranks contigs = ", are_any_non_zero_ranks);
  print_debug(_HERE_, " num mmi keys = ", mmi_uniq.size());

  MMI::T_multi_value_map pre_multi_map; // multi map: minimizer -> vector of 2 or more locations
  MMI::T_uniq_value_map pre_uniq_map;   // uniq map : minimizer -> one unique location

  assert(is_mmi_uniq_map_valid(mmi_uniq, num_segments));

  populate_uniq_and_multi_mmi_maps(pre_uniq_map, pre_multi_map, mmi_uniq);

  print_debug(_HERE_, " GFA fn = ", gfa.fn);
  phmap::flat_hash_set<std::string> alts;
  ::read_alts_file(alts, gfa.fn);

  if (are_any_non_zero_ranks && alts.size() > 0)
  {
    // Make a set of the alt segments we want to index
    std::vector<std::string> alt_segment_names;

    for (auto it = alt_segment_names.begin(); it != alt_segment_names.end(); ++it)
      alts.insert(*it);

    for (int rid{0}; rid < num_segments; ++rid)
    {
      gfa_seg_t const & segment = gfa.get_segment(rid);

      if (segment.rank == 0)
        continue;

      if (alts.find(segment.name) == alts.end())
      {
        print_debug(_HERE_, " skipping alt ", segment.name);
        continue;
      }

      print_debug(_HERE_, " sketching alt ", segment.name);
      std::vector<u64_pair> new_sketches{};
      sketch_rid(new_sketches, gfa, rid, k, w);
      auto const new_sketch_begin = new_sketches.cbegin();
      auto const new_sketch_end = new_sketches.cend();

      for (auto new_sketch_it = new_sketch_begin; new_sketch_it != new_sketch_end; ++new_sketch_it)
      {
        if (is_sketch_value_valid(gfa, new_sketch_it->second))
        {
          if (is_sketch_with_unique_context(new_sketch_begin, new_sketch_end, new_sketch_it, pre_uniq_map))
          {
            mmi_uniq.insert(mmi_uniq.end(), *new_sketch_it);
          }
        }
        else
        {
          ++num_invalid_sketches;
        }
      }

      if (num_invalid_sketches > 0)
      {
        print_debug(_HERE_, " num invalid sketches in rid = ", rid, " were ", num_invalid_sketches);
        num_invalid_sketches = 0;
      }
    }
  }

  print_debug(_HERE_, " num mmi keys = ", mmi_uniq.size());
  MMI::T_haplotypes haplotypes; // Store the set of haplotypes from the input VCF file

  if (vcf_fn.size() > 0)
  {
    int64_t const num_sketches_before{static_cast<int64_t>(mmi_uniq.size())};
    are_any_non_zero_ranks = true;
    get_minimizer_of_variant_data(mmi_uniq, pre_uniq_map, haplotypes, gfa, vcf_fn, k, w);

    // Add additional haplotypes on stable segments
    print_info(_HERE_,
               " num additional sketches from other haplotypes = ",
               static_cast<int64_t>(mmi_uniq.size()) - num_sketches_before);
  }

  // Reuse allocated memory
  MMI::T_multi_value_map mmi_multi_map = std::move(pre_multi_map);
  MMI::T_uniq_value_map mmi_uniq_map = std::move(pre_uniq_map);

  if (are_any_non_zero_ranks)
  {
    mmi_multi_map.clear();
    mmi_uniq_map.clear();

    assert(is_mmi_uniq_map_valid(mmi_uniq, num_segments));
    populate_uniq_and_multi_mmi_maps(mmi_uniq_map, mmi_multi_map, mmi_uniq);
    mmi_uniq.clear();
  }

  // All hits are sorted such that binary search is possible
  sort_hits_and_remove_keys_with_too_many_hits(mmi_multi_map, mmi_uniq_map, gfa);

  // Generate the MiniMizer Index (MMI)
  MMI mmi;
  mmi.haplotypes = std::move(haplotypes);

  for (auto it = mmi_multi_map.begin(); it != mmi_multi_map.end(); ++it)
  {
    assert(it->second.size() > 1);

    // non-unique hits only
    uint64_t const new_vals_start_i = mmi.values.size(); // start index for the new values
    uint64_t const new_vals_amount = it->second.size();  // the amount of new values
    mmi.values.insert(mmi.values.end(), it->second.begin(), it->second.end());

    assert(static_cast<double>(new_vals_start_i) < pow(2, 35));
    assert(static_cast<double>(new_vals_amount) < pow(2, 28));

    // 1 bit (non-unique?) | 35 bits (start index) | 28 bits (amount)
    uint64_t const val = (1ull << 63) | (new_vals_start_i << 28) | new_vals_amount;
    assert((val << 1) >> 29 == new_vals_start_i);
    assert((val << 36) >> 36 == new_vals_amount);

    auto find_it = mmi_uniq_map.find(it->first);
    assert(find_it != mmi_uniq_map.end());
    assert(find_it->second == std::numeric_limits<uint64_t>::max());
    find_it->second = val; // replace missing value with a range
  }

  mmi.map = std::move(mmi_uniq_map);

#ifndef NDEBUG
  for (auto key_val : mmi.map)
  {
    assert(key_val.second != std::numeric_limits<uint64_t>::max());
  }
#endif // NDEBUG

  return mmi;
}

} // namespace weaver
