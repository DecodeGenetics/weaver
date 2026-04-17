/*!
 * @file sr_seed.cpp
 * @brief Implements the SRSeed class.
 */

#include "sr_seed.hpp"

#include <array>   // std::array
#include <cassert> // assert
#include <numeric> // std::iota
#include <sstream> // std::ostringstream
#include <string>  // std::string

#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/pairwise_alignment.hpp>

#include <weaver/constants.hpp>

#include "gfa_arc.hpp"
#include "gfa_location.hpp"
#include "hashmap.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "read_sketch.hpp"
#include "segment.hpp"
#include "sequence_utils.hpp"
#include "sketch.hpp"
#include "sketch_cache.hpp"
#include "sketch_to_string.hpp"
#include "sketch_value.hpp"
#include "sr_alignment.hpp"
#include "sr_cigar.hpp"
#include "sr_seed_extend.hpp"
#include "sr_seed_pair.hpp"

namespace
{
/*
void filter_seeds(std::vector<weaver::SRSeed> & seeds, int const highest_score_seed)
{
  int seed_size_filter{60}; // TODO make an option

  if (seeds.size() <= 12 || (seeds.size() <= 42 && highest_score_seed <= seed_size_filter))
    return;

  auto seed_filter = [highest_score_seed, &seed_size_filter](weaver::SRSeed const & s) -> bool { //
    return s.get_est_score() < (highest_score_seed - seed_size_filter);
  };

  seeds.erase(std::remove_if(seeds.begin(), seeds.end(), seed_filter), seeds.end());

  while (seeds.size() > 42 && seed_size_filter >= 15) // TODO make an option
  {
    seed_size_filter -= 10;
    seeds.erase(std::remove_if(seeds.begin(), seeds.end(), seed_filter), seeds.end());
  }

  while (seeds.size() > 42 && seed_size_filter >= 1) // TODO make an option
  {
    --seed_size_filter;
    seeds.erase(std::remove_if(seeds.begin(), seeds.end(), seed_filter), seeds.end());
  }
}
*/

template <typename Tit, typename T>
Tit sketch_lower_bound(weaver::GFA const & gfa, Tit first, Tit last, T const val)
{
  Tit it;
  typename std::iterator_traits<Tit>::difference_type count;
  typename std::iterator_traits<Tit>::difference_type step;
  count = std::distance(first, last);

  while (count > 0)
  {
    it = first;
    step = count / 2;
    std::advance(it, step);

    if (val > static_cast<T>(gfa.get_approximate_stable_position(*it)))
    {
      first = ++it;
      count -= step + 1;
    }
    else
    {
      count = step;
    }
  }

  return first;
}

template <typename Tit, typename T>
Tit sketch_upper_bound(weaver::GFA const & gfa, Tit first, Tit last, T const val)
{
  Tit it;
  typename std::iterator_traits<Tit>::difference_type count;
  typename std::iterator_traits<Tit>::difference_type step;
  count = std::distance(first, last);

  while (count > 0)
  {
    it = first;
    step = count / 2;
    std::advance(it, step);

    if (val >= static_cast<T>(gfa.get_approximate_stable_position(*it)))
    {
      first = ++it;
      count -= step + 1;
    }
    else
    {
      count = step;
    }
  }

  return first;
}

template <typename Tit, typename T>
Tit sketch_linear_upper_bound(weaver::GFA const & gfa, Tit first, Tit last, T const val)
{
  Tit it = first;

  while (it != last)
  {
    if (val >= static_cast<T>(gfa.get_approximate_stable_position(*it)))
      ++it;
    else
      break;
  }

  return it;
}

template <std::size_t N, typename Tit>
inline std::array<Tit, N> max_n_elements(Tit begin, Tit end)
{
  std::array<Tit, N> max_its;
  max_its.fill(end);

  auto check_index = [&](int i, Tit it) -> bool
  {
    if (max_its[i] == end || it->second > max_its[i]->second)
    {
      for (int j{i + 1}; j < static_cast<int>(N); ++j)
        max_its[j] = max_its[j - 1];

      max_its[i] = it;
      return true;
    }

    return false;
  };

  while (begin != end)
  {
    for (int i{0}; i < static_cast<int>(N); ++i)
    {
      if (check_index(i, begin))
        break;
    }

    ++begin;
  }

  return max_its;
}

void get_seeds_with_approx_pos(std::vector<weaver::SRSeed> & seeds, //
                               weaver::ReadSketch const & read_sketch,
                               bool const is_any_extending_init,
                               weaver::GFA const & gfa,
                               weaver::T_icu const & icu)
{
  using namespace weaver;

  std::vector<uint64_t> const & mmi_results = read_sketch.mmi_results;
  uint64_t const read_value = read_sketch.second;
  int const read_pos = sketch_value_pos(read_value);
  int constexpr NO_HIT{-3};
  int constexpr EXACT_HIT{-2};
  int constexpr APPROXIMATE_HIT{-1};
  int const number_of_seeds{static_cast<int>(seeds.size())};
  std::vector<int> hits(number_of_seeds, NO_HIT);

  for (auto it = mmi_results.begin(); it != mmi_results.end(); ++it)
  {
    uint64_t ref_value = *it;
    assert(sketch_value_rid(ref_value) < gfa.get_num_segments());
    assert(sketch_value_rid(read_value) == gfa.get_num_segments());

    if ((read_value & 1) == 1)
      ref_value ^= 1; // toggle least significant bit when read_value is reversed

    bool is_any_extending{is_any_extending_init};

    // Check if this ref_value extends a previous seed
    for (int s{0}; s < number_of_seeds; ++s)
    {
      auto & hit = hits[s];

      if (hit != NO_HIT) // hit already found, only one is expected
        continue;

      SRSeed & seed = seeds[s];

#ifndef NDEBUG
      // seed should be valid before extending
      {
        bool const is_seed_valid = seed.is_valid();

        if (not is_seed_valid)
        {
          print_error(_HERE_, " before extending, the seed was not valid\n", seed.to_string());
          assert(false);
        }
      }

#endif // NDEBUG

      int const seed_read_pos = sketch_value_pos(seed.end_read_value);
      int const read_distance{read_pos - seed_read_pos};
      assert(read_distance > 0);

      // Can the seed be extended exactly?
      if (seed_extend_if_exact_distance(gfa, seed, icu, ref_value, read_value, read_distance))
      {
        hit = EXACT_HIT;
        is_any_extending = true;
        // print_debug(_HERE_, " !!! Extending !!! new length=", seed.length, " ", read_distance, " ",
        // seed.to_string()); highest_score_seed = std::max(highest_score_seed, seed.get_est_score());

#ifndef NDEBUG
        // seed should be valid before extending
        bool const is_seed_valid = seed.is_valid();

        if (not is_seed_valid)
        {
          print_error(_HERE_, " after extending, the seed was not valid\n", seed.to_string());
          assert(false);
        }

#endif // NDEBUG
      }
    }

    if (not is_any_extending && /*approximate matching=*/true)
    {
      for (int s{0}; s < number_of_seeds; ++s)
      {
        auto & hit = hits[s];

        if (hit != NO_HIT)
        {
          // this seed already has found an extension/a hit
          continue;
        }

        SRSeed & seed = seeds[s];

        if (seed.is_extending_cut == false)
          continue; // seed has not been cut, do not check for approximate match until it has

        assert(seed.is_valid()); // seed should be valid before extending

        int const seed_read_pos = sketch_value_pos(seed.end_read_value);
        int const read_distance{read_pos - seed_read_pos};
        assert(read_distance > 0);

        // Can the seed be extended approximately?
        if (seed_extend_if_approximate_distance(gfa, seed, icu, ref_value, read_value, read_distance))
        {
          hit = APPROXIMATE_HIT;
          is_any_extending = true;

          print_debug(_HERE_,
                      " !!! Approximate extending !!! new length=",
                      seed.length,
                      " ",
                      read_distance,
                      " ",
                      seed.to_string());

          // highest_score_seed = std::max(highest_score_seed, seed.get_est_score());

          // Check if this now has the largest seed score
          assert(seed.is_valid()); // seed should be valid after extending
        }
      }
    }

    // After going through all the current seeds and we never found an extension, we add this as a new seed
    if (not is_any_extending)
    {
      SRSeed new_seed(ref_value, read_value & ~1ull); // unset strand in read_value
      // print_debug(_HERE_, " !!! New seed !!! ", new_seed.to_string());
      seeds.push_back(std::move(new_seed));
    }
  }
  // }

  for (int s{0}; s < number_of_seeds; ++s)
  {
    if (hits[s] == NO_HIT)
      seeds[s].is_extending_cut = true;
  }
}

std::vector<weaver::ReadSketch> filter_read_sketches(std::vector<weaver::ReadSketch> const & read_sketches,
                                                     std::vector<std::pair<int64_t, int64_t>> const & intervals,
                                                     weaver::GFA const & gfa)
{
  using namespace weaver;
  std::vector<weaver::ReadSketch> filtered_read_sketches;

  for (ReadSketch const & read_sketch : read_sketches)
  {
    ReadSketch new_read_sketch(read_sketch.first, read_sketch.second);
    std::vector<uint64_t> const & mmi_results = read_sketch.mmi_results;
    std::vector<uint64_t> new_mmi_results;
    std::vector<uint64_t>::const_iterator const first = mmi_results.begin();
    std::vector<uint64_t>::const_iterator const last = mmi_results.end();

    for (auto const & interval : intervals)
    {
      assert(interval.first != interval.second);
      assert(interval.first < interval.second);
      assert(interval.first + 1 < interval.second);

      auto begin_it = sketch_lower_bound(gfa, first, last, interval.first);
      auto end_it = sketch_linear_upper_bound(gfa, begin_it, last, interval.second - 1);
      new_read_sketch.mmi_results.insert(new_read_sketch.mmi_results.end(), begin_it, end_it);
    }

    if (new_read_sketch.mmi_results.size() > 0)
      filtered_read_sketches.push_back(new_read_sketch);
  }

  return filtered_read_sketches;
}

void call_get_seeds_with_approx_pos(std::vector<weaver::ReadSketch> const & read_sketches,
                                    std::vector<weaver::SRSeed> & seeds,
                                    std::vector<std::pair<int64_t, int64_t>> const & intervals,
                                    weaver::GFA const & gfa,
                                    weaver::T_icu const & icu)
{
  // Extracts all read sketches that have matches within the selected buckets
  std::vector<weaver::ReadSketch> filtered_read_sketches = filter_read_sketches(read_sketches, intervals, gfa);

  for (weaver::ReadSketch const & read_sketch : filtered_read_sketches)
  {
    get_seeds_with_approx_pos(seeds,
                              read_sketch,
                              /*is_any_extending_init=*/false,
                              gfa,
                              icu);
  }
}

} // namespace

namespace weaver
{
using Tseen_hits = phmap::flat_hash_map<int64_t, int>;

SRSeed::SRSeed(uint64_t ref_val, uint64_t read_val) noexcept :
  begin_ref_value(ref_val ^ 1ull), // store the begin/original values with
                                   // the other strand
  end_ref_value(ref_val),          // store most recent values of ref and read
  begin_read_value(read_val ^ 1ull),
  end_read_value(read_val),
  length(1)
{
}

SRSeed::SRSeed(SRSeed const & o) noexcept :
  begin_ref_value(o.begin_ref_value),
  end_ref_value(o.end_ref_value),
  begin_read_value(o.begin_read_value),
  end_read_value(o.end_read_value),
  score(o.score),
  length(o.length),
  begin_unaccounted_read_bases(o.begin_unaccounted_read_bases),
  end_unaccounted_read_bases(o.end_unaccounted_read_bases),
  num_cuts(o.num_cuts),
  is_extending_cut(o.is_extending_cut),
  arcs(o.arcs),
  cig(nullptr)
{
  if (o.cig != nullptr)
    cig = std::make_unique<std::vector<paw::Cigar>>(*o.cig);
}

SRSeed & SRSeed::operator=(SRSeed const & o) noexcept
{
  begin_ref_value = o.begin_ref_value;
  end_ref_value = o.end_ref_value;
  begin_read_value = o.begin_read_value;
  end_read_value = o.end_read_value;
  score = o.score;
  length = o.length;
  begin_unaccounted_read_bases = o.begin_unaccounted_read_bases;
  end_unaccounted_read_bases = o.end_unaccounted_read_bases;
  num_cuts = o.num_cuts;
  is_extending_cut = o.is_extending_cut;
  arcs = o.arcs;

  if (o.cig != nullptr)
    cig = std::make_unique<std::vector<paw::Cigar>>(*o.cig);

  return *this;
}

void SRSeed::extend_end(int l, uint64_t new_end_ref_value, uint64_t new_end_read_value)
{
  assert(l >= 0);
  length += l;
  end_ref_value = new_end_ref_value;
  end_read_value = new_end_read_value;
  // assert(get_read_length() == length);
}

void SRSeed::clear_essentials()
{
  begin_ref_value = 0;
  end_ref_value = 0;
  begin_read_value = 0;
  end_read_value = 0;
  score = 0;
  length = 0;
}

bool SRSeed::is_begin_reversed() const
{
  return (begin_ref_value & 1) != (begin_read_value & 1);
}

bool SRSeed::is_end_reversed() const
{
  return (end_ref_value & 1) != (end_read_value & 1);
}

bool SRSeed::is_empty() const
{
  return length == 0;
}

std::string SRSeed::to_string() const
{
  std::ostringstream ss;

  ss << " begin ref=" << sketch_value_to_string(begin_ref_value)   //
     << " end_ref=" << sketch_value_to_string(end_ref_value)       //
     << " begin read=" << sketch_value_to_string(begin_read_value) //
     << " end_read=" << sketch_value_to_string(end_read_value)     //
     << " length=" << length                                       //
     << " is_begin_reversed=" << is_begin_reversed()               //
     << " is_end_reversed=" << is_end_reversed()                   //
     << " num_cuts=" << num_cuts;                                  //

  for (int l{0}; l < static_cast<int>(arcs.size()); ++l)
  {
    ss << " arcs[" << l << "] v->w=" << (arcs[l]->v_lv >> 33) << "|" << ((arcs[l]->v_lv >> 32) & 1);
    ss << " -> " << (arcs[l]->w >> 1) << "|" << (arcs[l]->w & 1);
  }

  if (cig)
    ss << " cigar=" << inv_cigar2string(cig->rbegin(), cig->rend());

  return ss.str();
}

int SRSeed::get_read_length() const
{
  return static_cast<int>(end_read_value > 0) + static_cast<int>((static_cast<uint32_t>(end_read_value) >> 1) -
                                                                 (static_cast<uint32_t>(begin_read_value) >> 1));
}

int SRSeed::get_est_score() const
{
  return static_cast<int>((static_cast<uint32_t>(end_read_value) >> 1) -
                          (static_cast<uint32_t>(begin_read_value) >> 1)) -
         4 * num_cuts + 1;
}

int SRSeed::get_score() const
{
  return score;
}

void SRSeed::align_to_graph(std::string_view query_seq)
{
  Options const & copts = *(Options::const_instance());

  std::string ref = this->get_ref_sequence();
  std::string read = this->get_read_sequence(query_seq);

  int mem_read_start{0};
  int mem_ref_start{0};
  int mem_length{0};

  find_mem(mem_read_start, mem_ref_start, mem_length, read, ref);

  if (mem_length == 0)
  {
    print_debug(_HERE_, " zero length mem, read=", read, " ref=", ref, " seed=", to_string());
    return;
  }

  score = copts.match * mem_length;

  // shrink the seed
  {
    int not_shrinkable_begin{0};

    if (mem_ref_start > 0)
      not_shrinkable_begin = shrink_reference_begin(mem_ref_start);

    mem_read_start -= not_shrinkable_begin;

    if (mem_read_start > 0)
      shrink_read_begin(mem_read_start);

    int const end_ref_diff = static_cast<int>(ref.size()) - (mem_ref_start + mem_length);
    assert(end_ref_diff >= 0);

    int not_shrinkable_end{0};

    if (end_ref_diff > 0)
      not_shrinkable_end = shrink_reference_end(end_ref_diff);

    int end_read_diff = static_cast<int>(read.size()) - (mem_read_start + mem_length);
    assert(end_read_diff >= 0);
    end_read_diff -= not_shrinkable_end;

    if (end_read_diff > 0)
      shrink_read_end(end_read_diff);

#ifndef NDEBUG
    // Make sure that the mem is correct
    std::string new_ref = this->get_ref_sequence();
    std::string new_read = this->get_read_sequence(query_seq);

    if (not_shrinkable_end == 0 && new_ref != new_read)
    {
      print_info(_HERE_, " new_ref != new_read ", new_ref, " != ", new_read);
      print_info(_HERE_, " ref=", ref);
      print_info(_HERE_, " read=", read);
      print_info(_HERE_, " mem readcut=", mem_read_start, " ", end_read_diff);
      print_info(_HERE_, " mem refcut=", mem_ref_start, " ", end_ref_diff);
      assert(ref == read);
    }
#endif // NDEBUG
  }
}

bool SRSeed::has_indel_in_cigar() const
{
  if (cig == nullptr)
    return false;

  for (auto it = cig->begin(); it != cig->end(); ++it)
  {
    if (it->operation == paw::CigarOperation::INSERTION || it->operation == paw::CigarOperation::DELETION)
      return true;
  }

  return false;
}

bool SRSeed::has_deletion_and_insertion_in_cigar() const
{
  if (cig == nullptr)
    return false;

  bool has_del{false};
  bool has_ins{false};

  for (auto it = cig->begin(); it != cig->end(); ++it)
  {
    has_del |= it->operation == paw::CigarOperation::DELETION;
    has_ins |= it->operation == paw::CigarOperation::INSERTION;
  }

  return has_del && has_ins;
}

int SRSeed::shrink_reference_begin(int by)
{
  GFALocation begin_ref_location(this->begin_ref_value ^ 1ull); // begin ref value with flipped strand

  std::vector<gfa_arc_t const *> new_arcs;
  int new_by = begin_ref_location.advance_when_same_contig(by, new_arcs);
  begin_ref_location.flip_strand(); // flip strand again

  // update to new value
  this->begin_ref_value = begin_ref_location.get_value();
  this->length -= (by - new_by);

  // remove walked arcs
  for (auto new_arc_ptr : new_arcs)
  {
    print_debug(_HERE_, " new arc=", arc_to_string(*new_arc_ptr));

    for (auto old_arc_ptr_it = this->arcs.begin(); old_arc_ptr_it != this->arcs.end(); ++old_arc_ptr_it)
    {
      if (is_complement_arc(new_arc_ptr, *old_arc_ptr_it) || is_same_arc(new_arc_ptr, *old_arc_ptr_it))
      {
        print_debug(_HERE_,
                    " same as old arc=",
                    arc_to_string(*(*old_arc_ptr_it)),
                    " strand=",
                    !(begin_ref_location.strand));

        this->arcs.erase(old_arc_ptr_it);
        break;
      }
    }
  }

  if (this->cig)
  {
    assert(this->cig->size() > 0);
    assert(this->cig->back().operation == paw::CigarOperation::INSERTION);

    if (new_by == 0)
      this->cig->pop_back();
    else
      this->cig->back().count -= (new_by - by);
  }

  return new_by;
}

int SRSeed::shrink_reference_end(int by)
{
  GFALocation end_ref_location(this->end_ref_value ^ 1ull); // end ref value with flipped strand

  std::vector<gfa_arc_t const *> new_arcs;
  int new_by = end_ref_location.advance_when_same_contig(by, new_arcs);
  end_ref_location.flip_strand(); // flip strand again

  // update to new value
  this->end_ref_value = end_ref_location.get_value();
  this->length -= (by - new_by);

  // remove walked arcs
  for (auto new_arc_ptr : new_arcs)
  {
    print_debug(_HERE_, " new arc=", arc_to_string(*new_arc_ptr));

    for (auto old_arc_ptr_it = this->arcs.begin(); old_arc_ptr_it != this->arcs.end(); ++old_arc_ptr_it)
    {
      if (is_complement_arc(new_arc_ptr, *old_arc_ptr_it) || is_same_arc(new_arc_ptr, *old_arc_ptr_it))
      {
        print_debug(_HERE_,
                    " same as old arc=",
                    arc_to_string(*(*old_arc_ptr_it)),
                    " strand=",
                    end_ref_location.strand);

        this->arcs.erase(old_arc_ptr_it);
        break;
      }
    }
  }

  if (this->cig)
  {
    assert(this->cig->size() > 0);
    assert(this->cig->front().operation == paw::CigarOperation::INSERTION);

    if (new_by == 0)
      this->cig->erase(this->cig->begin());
    else
      this->cig->front().count -= (by - new_by);
  }

  return new_by;
}

void SRSeed::shrink_read_begin(int by)
{
  assert(by >= 0);
  this->begin_read_value += (static_cast<uint64_t>(by) << 1ull);

  if (this->cig)
  {
    assert(this->cig);
    assert(this->cig->size() > 0);
    assert(this->cig->back().operation == paw::CigarOperation::DELETION);

    this->begin_read_value += (static_cast<uint64_t>(by) << 1);
    this->cig->pop_back();
  }
}

void SRSeed::shrink_read_end(int by)
{
  assert(by >= 0);
  this->end_read_value -= (static_cast<uint64_t>(by) << 1ull);

  if (this->cig)
  {
    assert(this->cig->size() > 0);
    assert(this->cig->front().operation == paw::CigarOperation::DELETION);

    this->cig->erase(this->cig->begin());
  }
}

std::string SRSeed::get_ref_sequence() const
{
  GFALocation begin(begin_ref_value ^ 1ull);
  GFALocation end(end_ref_value);
  std::string ref_seq;
  begin.advance_until_and_get_sequence(end, arcs, ref_seq);
  return ref_seq;
}

std::string SRSeed::get_read_sequence(std::string_view query_seq) const
{
  assert(query_seq.size() > 0);
  int const remaining_end_query_size = query_seq.size() - (sketch_value_pos(end_read_value) + 1);
  int const remaining_begin_query_size = sketch_value_pos(begin_read_value);

  query_seq.remove_suffix(remaining_end_query_size);
  query_seq.remove_prefix(remaining_begin_query_size);
  assert(query_seq.size() > 0);
  std::string query_seq_str(query_seq);
  return query_seq_str;
}

bool SRSeed::do_you_see_me(GFA const & gfa, T_icu const & icu, SRSeed const & to, int const max_dist) const
{
  SRSeed const & from = *this;
  GFALocation read1(from.begin_ref_value);
  read1.advance_on_segment(sketch_value_pos(from.begin_read_value));
  read1.flip_strand();
  GFALocation read2(to.begin_ref_value);
  read2.advance_on_segment(sketch_value_pos(to.begin_read_value));

  return is_within_distance(gfa, icu, read1.get_value(), read2.get_value(), max_dist);
}

bool SRSeed::is_valid() const
{
  if (GFALocation::gfa == nullptr)
  {
    print_warning(_HERE_, " GFALocation::gfa not set");
    return false;
  }

  GFALocation begin(begin_ref_value ^ 1ull);
  GFALocation end(end_ref_value);

  if (!begin.is_valid())
  {
    print_warning(_HERE_, " bad begin location in seed.");
    return false;
  }

  if (!end.is_valid())
  {
    print_warning(_HERE_, " bad end location in seed.");
    return false;
  }

  if (end.rid == begin.rid &&                                 //
      static_cast<int>(1l + end.pos - begin.pos) != length && //
      static_cast<int>(1l + begin.pos - end.pos) != length)   //
  {
    print_warning(_HERE_, " bad length ", length, " end.pos=", end.pos, " begin.pos=", begin.pos);
    return false;
  }

  // check path
  int rid = begin.rid;
  bool strand = begin.strand;

  for (gfa_arc_t const * arc_ptr : arcs)
  {
    gfa_arc_t const & arc = *(arc_ptr);

    if (static_cast<int>(arc.v_lv >> 33ull) != rid)
    {
      print_warning(_HERE_, " bad vertex id in arc. ", to_string());
      return false;
    }

    if (((arc.v_lv >> 32ull) & 1) != strand)
    {
      print_warning(_HERE_, " bad vertex strand in arc. ", to_string());
      return false;
    }

    rid = arc.w >> 1;
    strand = arc.w & 1;
  }

  if (rid != end.rid)
  {
    print_warning(_HERE_, " we didn't find the end rid. ", to_string());
    return false;
  }

  if (strand != end.strand)
  {
    print_warning(_HERE_, " we didn't find the end strand. ", to_string());
    return false;
  }

  std::string seq;
  begin.advance_until_and_get_sequence(end, arcs, seq);

  if (static_cast<int>(seq.size()) != length)
  {
    print_warning(_HERE_, " incorrect seed length ", (seq.size()), " != ", length);
    return false;
  }

  return true;
}

bool seed_est_score_order_gt(SRSeed const & s, SRSeed const & o)
{
  return s.get_est_score() > o.get_est_score() ||
         (s.get_est_score() == o.get_est_score() && s.end_read_value > o.end_read_value);
}

std::vector<int> get_seed_est_score_sorted_order_indices(std::vector<SRSeed> const & seeds)
{
  std::vector<int> indices(seeds.size());
  std::iota(indices.begin(), indices.end(), 0);
  auto seed_order_lambda = [&seeds](int i, int j) -> bool { return seed_est_score_order_gt(seeds[i], seeds[j]); };
  std::sort(indices.begin(), indices.end(), seed_order_lambda);
  return indices;
}

void get_seen_hits(Tseen_hits & seen_hits,
                   GFA const & gfa,
                   std::vector<ReadSketch> const & read_sketches,
                   std::size_t const max_hits = std::numeric_limits<std::size_t>::max())
{
  bool is_too_rare_result{false};

  for (auto const & read_sketch : read_sketches)
  {
    uint32_t const num_results = read_sketch.mmi_results.size();

    if (num_results == 0)
      continue; // no hits

    if (num_results > max_hits)
    {
      is_too_rare_result = true;
      continue; // too many hits
    }

    // TODO make fallback when no clz available
    int count = static_cast<int>(__builtin_clz(num_results)) - 20;

    if (count <= 0)
      count = 1;

    int64_t prev_approx_stable_pos{std::numeric_limits<int64_t>::min()};

    for (uint64_t mmi_value : read_sketch.mmi_results)
    {
      int64_t const approx_stable_pos = gfa.get_approximate_stable_position(mmi_value);

      if (approx_stable_pos != prev_approx_stable_pos)
      {
        assert(approx_stable_pos > prev_approx_stable_pos);
        seen_hits[approx_stable_pos] += count;
        prev_approx_stable_pos = approx_stable_pos;
      }
    }
  }

  if (is_too_rare_result && seen_hits.size() == 0)
  {
    // scan all hits
    for (auto const & read_sketch : read_sketches)
    {
      uint32_t const num_results = read_sketch.mmi_results.size();

      if (num_results == 0 || num_results > 65536)
        continue; // no hits or too many hits

      // TODO make fallback when no clz available
      int count = static_cast<int>(__builtin_clz(num_results)) - 20;

      if (count <= 0)
        count = 1;

      int64_t prev_approx_stable_pos{std::numeric_limits<int64_t>::min()};

      for (uint64_t mmi_value : read_sketch.mmi_results)
      {
        int64_t const approx_stable_pos = gfa.get_approximate_stable_position(mmi_value);

        if (approx_stable_pos != prev_approx_stable_pos)
        {
          assert(approx_stable_pos > prev_approx_stable_pos);
          seen_hits[approx_stable_pos] += count;
          prev_approx_stable_pos = approx_stable_pos;
        }
      }
    }
  }
}

void add_intervals(std::vector<std::pair<int64_t, int64_t>> & intervals, Tseen_hits const & seen_hits)
{
  Tseen_hits::const_iterator begin = seen_hits.cbegin();
  Tseen_hits::const_iterator end = seen_hits.cend();
  auto max_n_its = max_n_elements<NUM_HIT_INTERVALS>(begin, end);
  // int constexpr SCORE_THRESHOLD = 100;

  for (int i{0}; i < static_cast<int>(max_n_its.size()); ++i)
  {
    if (max_n_its[i] == end)
      break;

    std::pair<int64_t, int64_t> new_interval{max_n_its[i]->first - 1, max_n_its[i]->first + 2};

    // only add if it is not the same as the last one
    if (intervals.size() == 0 || intervals[intervals.size() - 1] != new_interval)
      intervals.push_back(new_interval);
  }
}

void merge_intervals(std::vector<std::pair<int64_t, int64_t>> & intervals)
{
  assert(intervals.size() > 0);
  assert(std::is_sorted(intervals.begin(), intervals.end()));
  std::vector<std::pair<int64_t, int64_t>> merged_intervals;
  merged_intervals.push_back(intervals[0]);
  int m{0}; // index to the merged interval

  for (int i{1}; i < static_cast<int>(intervals.size()); ++i)
  {
    auto const & interval = intervals[i];
    auto & merged_interval = merged_intervals[m];

    if (merged_interval.second >= interval.first)
    {
      if (interval.second > merged_interval.second)
        merged_interval.second = interval.second;
    }
    else
    {
      merged_intervals.push_back(interval);
      ++m;
    }
  }

  intervals = std::move(merged_intervals);
}

Tseen_hits const & merge_seen_hits(Tseen_hits & seen_hits1, Tseen_hits const & seen_hits2)
{
  // Add all hits from seen_hits2 into seen_hits1
  for (auto const & hit2 : seen_hits2)
  {
    auto insert_it = seen_hits1.insert(hit2);

    if (!insert_it.second)
      insert_it.first->second += hit2.second;
  }

  return seen_hits1;
}

#ifndef NDEBUG
void get_seeds_with_pair(std::vector<SRSeed> & seeds1,
                         std::vector<SRSeed> & seeds2,
                         std::vector<ReadSketch> const & read_sketches1, //
                         std::vector<ReadSketch> const & read_sketches2,
                         GFA const & gfa,
                         T_icu const & icu,
                         bool const is_debug)
#else
void get_seeds_with_pair(std::vector<SRSeed> & seeds1,
                         std::vector<SRSeed> & seeds2,
                         std::vector<ReadSketch> const & read_sketches1, //
                         std::vector<ReadSketch> const & read_sketches2,
                         GFA const & gfa,
                         T_icu const & icu,
                         bool const /*is_debug*/)
#endif
{
  Tseen_hits seen_hits1;
  Tseen_hits seen_hits2;

  get_seen_hits(seen_hits1, gfa, read_sketches1, 4096);
  get_seen_hits(seen_hits2, gfa, read_sketches2, 4096);

  if (seen_hits1.empty() && seen_hits2.empty())
  {
#ifndef NDEBUG
    if (is_debug)
      print_info(_HERE_, " No hits on either read.");
#endif // NDEBUG

    return; // no matches
  }

  std::vector<std::pair<int64_t, int64_t>> intervals;
  add_intervals(intervals, seen_hits1);
  add_intervals(intervals, seen_hits2);
  assert(intervals.size() > 0);

  if (intervals.size() > 1) // trivial case if there is only one interval, no need to merge hits or sort
  {
    // add_intervals(intervals, seen_hits1, seen_hits2);
    auto const & all_seen_hits = merge_seen_hits(seen_hits1, seen_hits2);
    add_intervals(intervals, all_seen_hits);
    std::sort(intervals.begin(), intervals.end());
    merge_intervals(intervals);
  }

  assert(intervals.size() > 0);

#ifndef NDEBUG
  if (is_debug)
  {
    print_info(_HERE_, " num merged intervals=", intervals.size());

    for (auto const & interval : intervals)
      print_info(_HERE_,
                 " ",
                 (interval.first >> 32ul),
                 "|",
                 static_cast<uint32_t>(interval.first),
                 " ",
                 (interval.second >> 32ul),
                 "|",
                 static_cast<uint32_t>(interval.second));
  }
#endif // NDEBUG

  call_get_seeds_with_approx_pos(read_sketches1, seeds1, intervals, gfa, icu);
  call_get_seeds_with_approx_pos(read_sketches2, seeds2, intervals, gfa, icu);
}

} // namespace weaver
