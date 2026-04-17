#include "edit.hpp"

#include <cassert>     // assert
#include <numeric>     // std::accumulate
#include <string>      // std::string
#include <string_view> // std::string_view
#include <utility>     // std::hash

namespace weaver
{
std::string Edit::to_string() const
{
  std::string str;

  str += "pos=";
  str += std::to_string(this->pos + 1); // to 1-based indexing
  str += " type=";
  str += this->type;
  str += " sequence=";
  str += std::string(this->seq.begin(), this->seq.end());

  return str;
}

std::string EditCalls::to_string() const
{
  std::string str;

  str += "pos=";
  str += std::to_string(this->pos + 1); // to 1-based indexing
  str += " type=";
  str += this->type;
  str += " sequence=";
  str += std::string(this->seq.begin(), this->seq.end());
  str += " num_calls=";
  str += std::to_string(this->calls.size());

  return str;
}

std::size_t EditHash::operator()(Edit const & e) const
{
  std::size_t h1 = std::hash<int>()(e.pos);
  std::size_t h2 = std::hash<char>()(e.type);
  std::size_t h3 = 42 + std::hash<std::string_view>{}(std::string_view{e.seq.data(), e.seq.size()});
  return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
}

void get_hap_match(std::vector<double> & alpha, std::vector<EditCalls> const & ec, int const edit)
{
  int const n = static_cast<int>(alpha.size());
  assert(n > 0);

  if (edit >= 0)
  {
    int const e = edit;
    assert(e >= 0);
    assert(e < static_cast<int>(ec.size()));
    std::vector<uint8_t> const & calls_to_add = ec[e].calls;
    assert(alpha.size() == calls_to_add.size() + 1);

    for (int c{1}; c < n; ++c)
      alpha[c] += (calls_to_add[c - 1] == 1);
  }
  else
  {
    alpha[0] += 1; // reference match
    int const no_e = -edit - 1;
    assert(no_e >= 0);
    assert(no_e < static_cast<int>(ec.size()));
    std::vector<uint8_t> const & calls_to_add = ec[no_e].calls;
    assert(alpha.size() == calls_to_add.size() + 1);

    for (int c{1}; c < n; ++c)
      alpha[c] += (calls_to_add[c - 1] == 0);
  }
}

void add_hap_matches(std::vector<double> & hap_matches,
                     std::vector<double> const & snp_hap_matches,
                     std::vector<double> const & nonsnp_hap_matches)
{
  if (nonsnp_hap_matches.size() == 0)
  {
    hap_matches = snp_hap_matches;
  }
  else if (snp_hap_matches.size() == 0)
  {
    hap_matches = nonsnp_hap_matches;
  }
  else
  {
    assert(nonsnp_hap_matches.size() == snp_hap_matches.size());
    hap_matches.resize(snp_hap_matches.size());

    for (int h{0}; h < static_cast<int>(snp_hap_matches.size()); ++h)
      hap_matches[h] = snp_hap_matches[h] + nonsnp_hap_matches[h];
  }
}

int get_hap_matches(std::vector<double> & hap_matches,
                    std::vector<EditCalls> const & ec,
                    std::vector<int> const & edits)
{
  if (edits.size() == 0)
    return 0;

  assert(ec.size() > 0);

  if (hap_matches.empty())
    hap_matches.resize(ec[0].calls.size(), 0);

  assert(hap_matches.size() == ec[0].calls.size());
  int const n = static_cast<int>(hap_matches.size());
  int ref_matches{0};

  for (int const edit : edits)
  {
    if (edit >= 0)
    {
      int const e = edit;
      assert(e >= 0);
      assert(e < static_cast<int>(ec.size()));
      std::vector<uint8_t> const & calls_to_add = ec[e].calls;
      assert(hap_matches.size() == calls_to_add.size());

      for (int c{0}; c < n; ++c)
        hap_matches[c] += (calls_to_add[c] == 1);
    }
    else
    {
      ++ref_matches;
      int const no_e = -edit - 1;
      assert(no_e >= 0);
      assert(no_e < static_cast<int>(ec.size()));
      std::vector<uint8_t> const & calls_to_add = ec[no_e].calls;
      assert(hap_matches.size() == calls_to_add.size());

      for (int c{0}; c < n; ++c)
        hap_matches[c] += (calls_to_add[c] == 0);
    }
  }

  return ref_matches;
}

int get_ref_matches(std::vector<int> const & edits)
{
  return std::accumulate(edits.begin(), edits.end(), 0, [](int a, int v) { return a + (v < 0); });
}

} // namespace weaver
