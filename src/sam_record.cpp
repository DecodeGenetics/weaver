#include "sam_record.hpp"

#include <paw/align/cigar.hpp>

#include "fastq_data.hpp"
#include "logging.hpp"
#include "sr_cigar.hpp"
#include "stable_contigs.hpp"

namespace weaver
{
SAMRecord::SAMRecord(FastqData const & fastq_data) :
  qname(fastq_data.name), extra_tags(fastq_data.comment), seq(fastq_data.seq), qual(fastq_data.qual)
{
  // Make sure there is no suffix '/1' or '/2' in the query read name.
  assert(qname.size() <= 2 || qname[qname.size() - 2] != '/');
}

void SAMRecord::append_cigar(int const count, paw::CigarOperation op)
{
  assert(count >= 0);
  weaver::append_cigar_back(cig, static_cast<uint32_t>(count), op);
}

void SAMRecord::append_cigar(paw::Cigar const & cigar)
{
  weaver::append_cigar_back(cig, cigar.count, cigar.operation);
}

int SAMRecord::soft_clip_cigar_begin_based_on_ref(int num_ref_bases)
{
  std::vector<paw::Cigar> cig_new;
  int num_query_bases_removed{0};
  auto cig_it = cig.begin();

  for (; cig_it != cig.end(); ++cig_it)
  {
    bool is_advancing_ref = paw::advances_ref(cig_it->operation);
    bool is_advancing_query = paw::advances_query(cig_it->operation);

    if (!is_advancing_ref)
    {
      if (is_advancing_query)
        num_query_bases_removed += cig_it->count;

      continue;
    }

    // advancing ref is true
    if (static_cast<int>(cig_it->count) >= num_ref_bases)
    {
      if (is_advancing_query)
        num_query_bases_removed += num_ref_bases;

      // add clipped read bases
      cig_new.emplace_back(num_query_bases_removed, paw::CigarOperation::SOFT_CLIP);

      if (static_cast<int>(cig_it->count) > num_ref_bases)
        weaver::append_cigar_back(cig_new, cig_it->count - num_ref_bases, cig_it->operation);

      ++cig_it;
      break;
    }

    if (is_advancing_query)
      num_query_bases_removed += cig_it->count;

    num_ref_bases -= cig_it->count;
  }

  // append all remaining cigar operations
  for (; cig_it != cig.end(); ++cig_it)
    weaver::append_cigar_back(cig_new, cig_it->count, cig_it->operation);

  print_debug(_HERE_,
              " old_cigar=",
              cigar2string(cig.begin(), cig.end()),
              " new_cigar=",
              cigar2string(cig_new.begin(), cig_new.end()));

  cig = std::move(cig_new);
  pos += num_ref_bases; // adjust position
  return num_query_bases_removed;
}

int SAMRecord::soft_clip_cigar_end_based_on_ref(int num_ref_bases)
{
  int num_query_bases_removed{0};
  auto cig_it = cig.rbegin();

  for (; cig_it != cig.rend(); ++cig_it)
  {
    bool is_advancing_ref = paw::advances_ref(cig_it->operation);
    bool is_advancing_query = paw::advances_query(cig_it->operation);

    if (!is_advancing_ref)
    {
      if (is_advancing_query)
        num_query_bases_removed += cig_it->count;

      continue;
    }

    // advancing ref is true
    if (static_cast<int>(cig_it->count) > num_ref_bases)
    {
      if (is_advancing_query)
        num_query_bases_removed += num_ref_bases;

      cig_it->count -= num_ref_bases;
      break;
    }
    else
    {
      if (is_advancing_query)
        num_query_bases_removed += cig_it->count;

      num_ref_bases -= cig_it->count;
    }
  }

  int const num_remove_operations = std::distance(cig.rbegin(), cig_it);
  assert(num_remove_operations >= 0);
  int const new_size = static_cast<int>(cig.size()) - num_remove_operations;

  cig.resize(new_size);
  weaver::append_cigar_back(cig, num_query_bases_removed, paw::CigarOperation::SOFT_CLIP);
  return num_query_bases_removed;
}

void SAMRecord::make_unmapped(SAMRecord & other_sam_record)
{
  // clear other record
  other_sam_record.flags |= SAMFlags::IS_MATE_UNMAPPED;
  // other_sam_record.mpos = MISSING_POS;

  // clear this record
  flags |= SAMFlags::IS_UNMAPPED;
  sfa_idx = -1;
  snid = -1;
  pos = MISSING_POS;
  mapq = 0;
  cig.clear();
  tlen = 0;
}

std::string SAMRecord::get_cigar() const
{
  if (cig.empty())
    return std::string(1, '*');

  return cigar2string(cig.begin(), cig.end());
}

int SAMRecord::get_cigar_query_length() const
{
  return weaver::get_cigar_query_length(cig);
}

int SAMRecord::get_cigar_reference_length() const
{
  return weaver::get_cigar_reference_length(cig);
}

uint64_t SAMRecord::get_sam_order() const
{
  assert(sfa_idx >= 0 || pos < 0);
  assert(sfa_idx < 0 || pos >= 0);

  if (sfa_idx < 0 || pos < 0)
    return std::numeric_limits<uint64_t>::max();

  assert(sfa_idx < static_cast<int>(stable_contigs.contigs.size()));
  assert(pos < stable_contigs.contigs[sfa_idx].max);
  return (static_cast<uint64_t>(sfa_idx) << 32) | static_cast<uint64_t>(pos);
}

int SAMRecord::get_reference_reach() const
{
  return pos + get_cigar_reference_length();
}

bool SAMRecord::is_double_clipped(int const min_clip_length) const
{
  int const n{static_cast<int>(cig.size())};

  if (n <= 2)
    return false;

  auto const & first = cig[0];
  auto const & last = cig[n - 1];

  return first.operation == paw::CigarOperation::SOFT_CLIP && //
         last.operation == paw::CigarOperation::SOFT_CLIP &&  //
         static_cast<int>(first.count + last.count) >= min_clip_length;
}

bool SAMRecord::is_unmapped_with_mapped_mate() const
{
  return ((flags & SAMFlags::IS_UNMAPPED) != 0u) && ((flags & SAMFlags::IS_MATE_UNMAPPED) == 0u);
}

bool SAMRecord::is_cigar_valid() const
{
  if ((flags & SAMFlags::IS_UNMAPPED) == 0)
  {
    for (auto c : cig)
    {
      if (c.count == 0)
      {
        print_warning(_HERE_, "CIGAR element with 0 count found: ", cigar2string(cig.begin(), cig.end()));
        return false;
      }
    }

    int const cigar_query_len = get_cigar_query_length();

    if (cigar_query_len != static_cast<int>(seq.size()))
    {
      print_warning(_HERE_,
                    " CIGAR string broken: ",
                    cigar2string(cig.begin(), cig.end()),
                    " (query_len=",
                    cigar_query_len,
                    ") with seq size=",
                    seq.size());

      return false;
    }
  }

  return true;
}

bool SAMRecord::is_valid() const
{
  if (qname.empty())
  {
    print_warning(_HERE_, " query name empty");
    return false;
  }

  if ((flags & SAMFlags::IS_UNMAPPED) == 0)
  {
    if (sfa_idx < 0)
    {
      print_warning(_HERE_, " stable FASTA index not available on a mapped read.");
      return false;
    }

    if (snid < 0)
    {
      print_warning(_HERE_, " stable contig name ID (snid) not available on a mapped read.");
      return false;
    }

    if (pos < 0)
    {
      print_warning(_HERE_, " pos<0 on a mapped read. qname=", qname);
      return false;
    }

    if (cig.empty())
    {
      print_warning(_HERE_, " cigar string empty.");
      return false;
    }

    if (!is_cigar_valid())
    {
      print_warning(_HERE_, " cigar not valid.");
      return false;
    }
  }

  return true;
}

bool SAMRecord::is_pair_valid(SAMRecord const & /*other*/) const
{
  if (!is_valid())
    return false;

  // if (pos != other.mpos || mpos != other.pos)
  // {
  //   print_warning(_HERE_, " mate pos and pos are different.", pos, " ", other.mpos, " ", other.pos, " ", mpos);
  //   return false;
  // }

  return true;
}

} // namespace weaver
