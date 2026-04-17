#include "variant.hpp"

#include <algorithm> // std::find
#include <charconv>  // std::from_chars
#include <vector>    // std::vector

#include "io.hpp"
#include "logging.hpp"
#include "sequence_utils.hpp"

namespace
{
void parse_genotype_field(std::vector<uint16_t> & calls, const char * seq, int const b, int const i)
{
  auto find_colon_ptr = std::find(seq + b, seq + i, ':');
  auto find_sep_ptr = std::find(seq + b, find_colon_ptr, '|');

  // try to find a '/' if no '|'
  if (find_sep_ptr == find_colon_ptr)
    find_sep_ptr = std::find(seq + b, find_colon_ptr, '/');

  // parse the first genotype
  {
    uint16_t gt{};
    auto ret = std::from_chars(seq + b, find_sep_ptr, gt);

    if (ret.ec != std::errc())
    {
      calls.push_back(weaver::Variant::MISSING_CALL); // Unable to parse genotype => Missing genotype
    }
    else
    {
      // print_info(_HERE_, " call ", gt);
      calls.push_back(gt);
    }
  }

  // parse the second genotype (if there is one)
  if (find_sep_ptr != find_colon_ptr)
  {
    uint16_t gt{};
    auto ret = std::from_chars(find_sep_ptr + 1, find_colon_ptr, gt);

    if (ret.ec != std::errc())
    {
      calls.push_back(weaver::Variant::MISSING_CALL);
    }
    else
    {
      // print_info(_HERE_, " call ", gt);
      calls.push_back(gt);
    }
  }
}

void read_small_variants(std::vector<weaver::Variant> & small_variants, const char * seq, int const length)
{
  weaver::Variant new_site;
  int b{0};
  int field{0};
  assert(length > 0);
  assert(seq[0] != '\t');
  // auto & new_site = small_variants.new_site;

  for (int i{1}; i < length; ++i)
  {
    if (seq[i] != '\t')
      continue;

    if (field == 1)
    {
      auto ret = std::from_chars(seq + b, seq + i, new_site.pos);

      if (ret.ec != std::errc())
      {
        print_warning(_HERE_, " Unable to parse position from VCF: ", std::string(seq + b, seq + i));
        return;
      }
    }
    else if (field == 3)
    {
      new_site.seqs.emplace_back(seq + b, seq + i);
      // new_site.ref.assign(seq + b, seq + i);

      if ((i - b) > weaver::Variant::MAX_SIZE)
        return;
    }
    else if (field == 4)
    {
      if (std::find(seq + b, seq + i, ',') != seq + i)
      {
        // multi allelic
        std::string alts(seq + b, seq + i);
        std::vector<std::string_view> spl_alts = weaver::split_string(alts, ',');
        int a{0};

        for (std::string_view const spl_alt : spl_alts)
        {
          if (spl_alt.size() > weaver::Variant::MAX_SIZE)
            return;

          assert((a + spl_alt.size()) <= alts.size());

          if (!std::all_of(alts.data() + a, alts.data() + a + spl_alt.size(), weaver::isACGT))
            return;

          new_site.seqs.emplace_back(alts.data() + a, alts.data() + a + spl_alt.size());
          // new_site.alts.emplace_back(alts.data() + a, alts.data() + a + spl_alt.size());
          a += spl_alt.size() + 1;
        }
      }
      else
      {
        // biallelic
        new_site.seqs.emplace_back(seq + b, seq + i);

        if (!std::all_of(seq + b, seq + i, weaver::isACGT))
          return;

        // new_site.alts.emplace_back(seq + b, seq + i);
      }
    }
    else if (field >= 9)
    {
      parse_genotype_field(new_site.calls, seq, b, i);
    }

    ++field;
    ++i;
    b = i;
  }

  parse_genotype_field(new_site.calls, seq, b, length);
  auto const num_sites{small_variants.size()};

  if (num_sites == 0)
  {
    small_variants.emplace_back(std::move(new_site));
  }
  else
  {
    auto const & last_site = small_variants[num_sites - 1u];
    bool const is_overlapping_previous = (last_site.pos + static_cast<int>(last_site.seqs[0].size())) > new_site.pos;

    if (not is_overlapping_previous)
      small_variants.emplace_back(std::move(new_site));
  }
}

} // namespace

namespace weaver
{
std::vector<Variant> get_variants_in_a_region(hts_file_ptr const & in_vcf,
                                              tbx_t_ptr const & in_tbx,
                                              hts_itr_t_ptr const & in_it)
{
  std::vector<Variant> small_variants;

  if (in_it == nullptr)
    return small_variants; // no variants in region

  if (in_vcf == nullptr || in_tbx == nullptr)
  {
    print_warning(_HERE_, " Unable to get small variants in region.");
    return small_variants;
  }

  // read small variants in the region
  {
    kstring_t str = {0, 0, 0};
    int ret = tbx_itr_next(in_vcf.get(), in_tbx.get(), in_it.get(), &str);

    while (ret > 0)
    {
      read_small_variants(small_variants, str.s, str.l);
      ret = tbx_itr_next(in_vcf.get(), in_tbx.get(), in_it.get(), &str);
    }

    free(str.s);
  }

  return small_variants;
}

} // namespace weaver
