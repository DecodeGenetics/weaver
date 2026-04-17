/*!
 * @file sam_writer.cpp
 * @brief Implements the SAMWriter class.
 */
#include "sam_writer.hpp"

#include <cstdio>
#include <fstream>
#include <gfa.h>
#include <iostream>
#include <memory>
#include <numeric>
#include <ostream>

#include <paw/align/cigar.hpp>

#include "constants.hpp"
#include "gfa.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "sam_order_line.hpp"
#include "sam_record.hpp"
#include "sequence_utils.hpp"
#include "stable_contigs.hpp"

namespace
{
// Stream deleter that does nothing (no ownership assumed).
static void stream_deleter_noop(std::ostream *)
{
}

// Stream deleter with default behaviour (ownership assumed).
static void stream_deleter_default(std::ostream * ptr)
{
  delete ptr;
}

// From samtools/bam_mate.c:226
int calc_mate_score(std::string_view qual)
{
  int constexpr MD_MIN_QUALITY{15 + 33}; // +33 because of ascii offset
  int score{0};

  for (char q : qual)
  {
    if (static_cast<int>(q) >= MD_MIN_QUALITY)
      score += (static_cast<int>(q) - 33);
  }

  return score;
}

} // namespace

namespace weaver
{
SAMWriter::SAMWriter(std::string const & fn) :
  sink{new std::ofstream{fn, std::ios::binary | std::ios::trunc}, ::stream_deleter_default}
{
}

SAMWriter::SAMWriter(std::ostream & other_stream) : sink{&other_stream, ::stream_deleter_noop}
{
}

void SAMWriter::close()
{
  sink = nullptr;
}

void SAMWriter::open(std::string const & fn)
{
  assert(sink == nullptr);
  sink.reset(new std::ofstream(fn, std::ios::binary | std::ios::trunc));
  sink.get_deleter() = stream_deleter_default; // set the default deleter
}

void SAMWriter::flush()
{
  assert(sink);
  sink->flush();
}

void SAMWriter::write_header(GFA const & /*gfa*/)
{
  assert(sink);
  // gfa_t const & g = gfa.get_graph();
  // int const num_stable_seqs = g.n_sseq;

  *sink << "@HD\tVN:1.6\tSO:coordinate\n";

  // Write contig names and their lengths to header
  for (Contig const & contig : stable_contigs.contigs)
    *sink << "@SQ\tSN:" << contig.name << "\tLN:" << contig.get_length() << "\n";

  /*
  for (int c{0}; c < num_stable_seqs; ++c)
  {
    auto const & sseq = g.sseq[c];
    std::string contig_name = sseq.name;
    int const contig_length = sseq.max - sseq.min;

    if (sseq.rank == 0)
    {
      *sink << "@SQ\tSN:" << contig_name         // contig name
            << "\tLN:" << contig_length << "\n"; // contig length
    }
    else
    {
      auto find_it = stable_contigs.contig2index.find(contig_name);

      if (find_it == stable_contigs.contig2index.end())
      {
        print_warning(_HERE_, " ignoring stable FASTA sequence with an unexpected contig name: ' ", contig_name,
  "'"); continue;
      }

      assert(find_it->second < static_cast<int>(stable_contigs.contigs.size()));

      // for (std::pair<int, int> const & partition : stable_contigs.contigs[find_it->second].partitions)
      // {
      //   *sink << "@SQ\tSN:" << contig_name                                // contig name
      //         << '_' << partition.first << '_' << partition.second        // rank>0 extra contig name
      //         << "\tLN:" << (partition.second - partition.first) << "\n"; // contig length
      // }
    }
  }
  */

  // read group
  Options const & copts = *(Options::const_instance());

  if (copts.read_group_header_line.size() > 0)
  {
    *sink << "@RG\tID:";

    std::string rg = copts.read_group_header_line.substr(8);
    std::size_t old_pos{0};
    std::size_t pos = rg.find("\\t");

    // store ID
    read_group_id = rg.substr(0, pos);

    print_debug(_HERE_, " read group id=", read_group_id);

    while (pos != std::string::npos)
    {
      *sink << rg.substr(old_pos, pos - old_pos) << '\t';
      old_pos = pos + 2;
      pos = rg.find("\\t", pos + 2);
    }

    *sink << rg.substr(old_pos) << '\n';
  }

  // User can specify extra header line to be inserted in the header
  if (copts.extra_header_lines.size() > 0)
  {
    if (copts.extra_header_lines[0] == '@')
    {
      auto spl_line = split_string(copts.extra_header_lines, '\\');
      int const n_fields = spl_line.size();
      assert(n_fields > 1);

      if (n_fields > 0)
        *sink << spl_line[0];

      for (int f{1}; f < n_fields; ++f)
      {
        if (spl_line[f].size() <= 1)
          continue;

        char const first = spl_line[f][0];
        spl_line[f].remove_prefix(1);

        if (first == 't')
          *sink << '\t' << spl_line[f];
        else if (first == 'n')
          *sink << '\n' << spl_line[f];
        else
          *sink << ' ' << spl_line[f];
      }

      *sink << '\n';
    }
    else
    {
      std::ifstream extra_header_lines_f(copts.extra_header_lines);

      if (not extra_header_lines_f.fail())
        *sink << extra_header_lines_f.rdbuf();
      else
        print_warning("Failed reading extra header lines from file: ", copts.extra_header_lines);
    }
  }

  // program line
  if (!copts.no_PG)
  {
    *sink << "@PG\tID:weaver\tPN:weaver\tVN:" //
          << weaver_VERSION_MAJOR << '.'      //
          << weaver_VERSION_MINOR << '.'      //
          << weaver_VERSION_PATCH << '-'      //
          << GIT_COMMIT_SHORT_HASH;

    if (copts.command_line.size() > 0)
      *sink << "\tCL:" << copts.command_line;

    *sink << '\n';
  }
}

void SAMWriter::write_line(std::string const & line)
{
  assert(sink);
  *sink << line;

  if (not read_group_id.empty())
    *sink << "\tRG:Z:" << read_group_id;

  int constexpr qual_field = 10;
  auto const spl_line = split_string(line, '\t', qual_field + 1);
  assert(qual_field < static_cast<int>(spl_line.size()));
  *sink << '\n';
}

void SAMWriter::write_line(std::vector<char> & data, std::string const & line) const
{
  assert(sink);
  data.insert(data.end(), line.begin(), line.end());

  if (not read_group_id.empty())
  {
    std::string rg = "\tRG:Z:" + read_group_id;
    data.insert(data.end(), rg.begin(), rg.end());
  }

  int constexpr qual_field = 10;
  auto const spl_line = split_string(line, '\t', qual_field + 1);
  assert(qual_field < static_cast<int>(spl_line.size()));
  data.push_back('\n');
}

void SAMWriter::write_sam_order_line(SAMOrderLine const & line)
{
  assert(sink);
  *sink << std::to_string(line.order) << "\t" << line.sam_line << '\n';
}

void SAMWriter::write_file(std::string const & fn)
{
  std::ifstream if_sam(fn);

  for (std::string line; std::getline(if_sam, line);)
    *sink << line << '\n';
}

void SAMWriter::write_sr_record(SAMRecord const & sam_record,
                                SAMRecord const & other_sam_record,
                                std::string const & cigar,
                                std::string const & other_cigar)
{
  assert(sink);
  assert(sam_record.seq.size() == sam_record.qual.size());
  std::string sam_line = get_sam_string(sam_record, other_sam_record, cigar, other_cigar);
  *sink << sam_line << '\n';
}

std::string get_sam_string(SAMRecord const & sam_record,
                           SAMRecord const & other_sam_record,
                           std::string const & cigar,
                           std::string const & other_cigar) noexcept
{
  // fields are appended to this sam_line
  std::string sam_line(sam_record.qname); // append qname

  // most lines are ~400 bytes if read size is 151
  sam_line.reserve(456);

  // write record
  sam_line += '\t';
  sam_line += std::to_string(sam_record.flags); // append SAM flags
  sam_line += '\t';

  // append contig name (which is called 'rname' in SAM) and position
  if (sam_record.sfa_idx < 0)
  {
    sam_line += SAMWriter::MISSING_SAM_FIELD; // rname
    sam_line += "\t0\t";                      // position==0 on unmapped
  }
  else
  {
    assert(sam_record.sfa_idx < static_cast<int>(stable_contigs.contigs.size()));
    sam_line += stable_contigs.contigs[sam_record.sfa_idx].name; // rname
    sam_line += '\t';
    sam_line += std::to_string(sam_record.pos + 1); // 1-based position
    sam_line += '\t';
  }

  sam_line += std::to_string(sam_record.mapq);
  sam_line += '\t';
  sam_line += cigar;
  sam_line += '\t';

  // append mate contig name ('rnext' in SAM) and the mate's position
  if (other_sam_record.sfa_idx < 0)
  {
    sam_line += SAMWriter::MISSING_SAM_FIELD; // rname
    sam_line += "\t0\t";                      // position==0 on unmapped
  }
  else
  {
    assert(other_sam_record.sfa_idx < static_cast<int>(stable_contigs.contigs.size()));
    sam_line += stable_contigs.contigs[other_sam_record.sfa_idx].name; // rname
    sam_line += '\t';
    sam_line += std::to_string(other_sam_record.pos + 1); // 1-based position
    sam_line += '\t';
  }

  sam_line += std::to_string(sam_record.tlen);
  sam_line += '\t';

  // write sequence and quality
  if (sam_record.seq.empty())
  {
    sam_line += SAMWriter::MISSING_SAM_FIELD;
    sam_line += '\t';
    sam_line += SAMWriter::MISSING_SAM_FIELD;
  }
  else
  {
    bool const is_seq_reversed = (sam_record.flags & SAMFlags::IS_SEQ_REVERSED) != 0u;

    if (is_seq_reversed)
    {
      for (auto seq_rev_it = sam_record.seq.rbegin(); seq_rev_it != sam_record.seq.rend(); ++seq_rev_it)
        sam_line += weaver::complement(*seq_rev_it);

      sam_line += '\t';
      sam_line.insert(sam_line.end(), sam_record.qual.rbegin(), sam_record.qual.rend());
    }
    else
    {
      sam_line += sam_record.seq;
      sam_line += '\t';
      sam_line += sam_record.qual;
    }
  }

  // write tags
  if (sam_record.num_edits != SAMRecord::MISSING_TAG)
  {
    sam_line.append("\tNM:i:", 6);
    sam_line += std::to_string(sam_record.num_edits);
  }

  if (sam_record.alignment_score != SAMRecord::MISSING_TAG)
  {
    sam_line.append("\tAS:i:", 6);
    sam_line += std::to_string(sam_record.alignment_score);
  }

  if (sam_record.weaver_score != SAMRecord::MISSING_TAG)
  {
    sam_line.append("\tWS:i:", 6);
    sam_line += std::to_string(sam_record.weaver_score);
  }

  if (sam_record.second_alignment_score != SAMRecord::MISSING_TAG)
  {
    if (sam_record.second_alignment_score < (SAMRecord::MISSING_TAG / 2))
    {
      print_warning(_HERE_, " extremely small secondary score = ", sam_record.second_alignment_score);
    }

    sam_line.append("\tXS:i:", 6);
    sam_line += std::to_string(sam_record.second_alignment_score);
  }

  // mate MQ
  sam_line.append("\tMQ:i:", 6);
  sam_line += std::to_string(other_sam_record.mapq);

  if (other_cigar.size() > 1)
  {
    sam_line.append("\tMC:Z:", 6);
    sam_line += other_cigar;
  }

  if (sam_record.extra_tags.size() > 5)
  {
    sam_line += '\t';
    sam_line += sam_record.extra_tags;
  }

  // write ms
  {
    sam_line.append("\tms:i:", 6);
    sam_line += std::to_string(calc_mate_score(other_sam_record.qual));
  }

  return sam_line;
}

} // namespace weaver
