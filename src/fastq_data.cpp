

#include "fastq_data.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <kseq.h>
#include <vector>
#include <zlib.h>

#include "logging.hpp"
#include "options.hpp"

namespace
{
void are_qnames_valid(std::string const & last, std::string const & curr)
{
  if (last.empty() || curr.empty())
  {
    print_error(" Encountered a FASTQ record with empty name.");
    std::exit(1);
  }

  int l = static_cast<int>(curr.size()) - 1;

  if (last.size() != curr.size() || (not std::equal(last.begin(), last.begin() + l, curr.begin(), curr.begin() + l)))
  {
    print_error(" Encountered a pair of FASTQ records that appear to come from different read pairs.",
                "FASTQ is likelily malformed. read1.name=",
                last,
                " read2.name=",
                curr);

    std::exit(1);
  }
}

void are_sequences_and_qual_valid(weaver::FastqData const & fastq_data)
{
  if (fastq_data.seq.empty())
  {
    print_error("Encountered a FASTQ record with empty sequence.");
    std::exit(1);
  }

  if (fastq_data.seq.size() != fastq_data.qual.size())
  {
    print_error("Encountered a FASTQ record with different sequence and quality lengths.");
    print_error("Problematic record has name = ", fastq_data.name);
    std::exit(1);
  }
}

/*!
 * @brief Transform quality like they do in NovaSeq (RTA3)
 *
 * @details
 * Based on:
 * https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf
 *
 * Bins:
 *    0 -  7: 2
 *    8 - 14: 12 (11)
 *   15 - 30: 23 (25)
 *   31+    : 37
 */
std::array<char, 64> constexpr rta3_transform_qual = {
  '#', // 0
  '#', // 1
  '#', // 2
  '#', // 3
  '#', // 4
  '#', // 5
  '#', // 6
  '#', // 7
  ',', // 8
  ',', // 9
  ',', // 10
  ',', // 11
  ',', // 12
  ',', // 13
  ',', // 14
  ':', // 15
  ':', // 16
  ':', // 17
  ':', // 18
  ':', // 19
  ':', // 20
  ':', // 21
  ':', // 22
  ':', // 23
  ':', // 24
  ':', // 25
  ':', // 26
  ':', // 27
  ':', // 28
  ':', // 29
  ':', // 30
  'F', // 31
  'F', // 32
  'F', // 33
  'F', // 34
  'F', // 35
  'F', // 36
  'F', // 37
  'F', // 38
  'F', // 39
  'F', // 40
  'F', // 41
  'F', // 42
  'F', // 43
  'F', // 44
  'F', // 45
  'F', // 46
  'F', // 47
  'F', // 48
  'F', // 49
  'F', // 50
  'F', // 51
  'F', // 52
  'F', // 53
  'F', // 54
  'F', // 55
  'F', // 56
  'F', // 57
  'F', // 58
  'F', // 59
  'F', // 60
  'F', // 61
  'F', // 62
  'F'  // 63
};

} // namespace

namespace weaver
{
FastqData::FastqData(kseq_t const & kseq) :
  name(kseq.name.s, kseq.name.l),
  comment(kseq.comment.s, kseq.comment.l),
  seq(kseq.seq.s, kseq.seq.l),
  qual(kseq.qual.s, kseq.qual.l)
{
  assert(kseq.seq.l == kseq.qual.l);

  // Expect the RG tag to appear first in the FQ comment, otherwise clear the comment
  if (comment.size() <= 5 || comment[0] != 'R' || comment[1] != 'G' || comment[2] != ':' || comment[3] != 'Z' ||
      comment[4] != ':')
  {
    comment.clear();
  }
}

std::vector<FastqData> * seq_store_get(kseq_t * seq)
{
  assert(seq != nullptr);
  std::vector<FastqData> * seqs = new std::vector<FastqData>();
  int const fastq_data_buffer_size = Options::const_instance()->fastq_data_buffer_size;
  assert(fastq_data_buffer_size % 2 == 0);
  seqs->reserve(fastq_data_buffer_size);

  // set new records
  for (int i{0}; i < fastq_data_buffer_size; i++)
  {
    if (kseq_read(seq) < 0)
      break;

    assert(seq);
    seqs->emplace_back(*seq);

#ifndef NDEBUG
    if ((i % 2) == 1)
    {
      FastqData const & last = (*seqs)[i - 1];
      FastqData const & curr = (*seqs)[i];
      ::are_qnames_valid(last.name, curr.name);
      ::are_sequences_and_qual_valid(last);
      ::are_sequences_and_qual_valid(curr);
    }
#endif // NDEBUG
  }

  return seqs;
}

std::vector<FastqData> * seq_store_get(kseq_t * seq1, kseq_t * seq2)
{
  assert(seq1 != nullptr);

  if (seq2 == nullptr)
    return seq_store_get(seq1);

  std::vector<FastqData> * seqs = new std::vector<FastqData>();
  int const fastq_data_buffer_size = Options::const_instance()->fastq_data_buffer_size;
  assert(fastq_data_buffer_size % 2 == 0);
  seqs->reserve(fastq_data_buffer_size);

  // set new records
  for (int i{0}; (i + 1) < fastq_data_buffer_size; i += 2)
  {
    if (kseq_read(seq1) < 0)
    {
      bool const no_seq2 = kseq_read(seq2) < 0;

      if (not no_seq2)
      {
        print_error(_HERE_, " end of fastq1 but fastq2 has more data.");
        FastqData fastq_data(*seq2);
        print_error(_HERE_, " fastq2 read name = ", fastq_data.name);
        std::exit(1);
      }
      else
      {
        break;
      }
    }
    else
    {
      bool const no_seq2 = kseq_read(seq2) < 0;

      if (no_seq2)
      {
        print_error(_HERE_, " end of fastq2 but fastq1 has more data.");
        FastqData fastq_data(*seq1);
        print_error(_HERE_, " fastq1 read name = ", fastq_data.name);
        std::exit(1);
      }
    }

    assert(seq1);
    assert(seq2);
    assert((seqs->size() % 2) == 0);
    seqs->emplace_back(*seq1);
    seqs->emplace_back(*seq2);

#ifndef NDEBUG
    // Sanity check fastq data
    assert(static_cast<int>(seqs->size()) == i + 2);
    FastqData const & last = (*seqs)[i];
    FastqData const & curr = (*seqs)[i + 1];
    ::are_qnames_valid(last.name, curr.name);
    ::are_sequences_and_qual_valid(last);
    ::are_sequences_and_qual_valid(curr);
#endif // NDEBUG
  }

  return seqs;
}

std::vector<std::vector<FastqData> *> get_new_fastq_buffers(
  std::vector<std::vector<FastqData>> && thread_reads_to_remap)
{
  std::vector<std::vector<FastqData> *> new_buffers;
  int const fastq_data_buffer_size = 2 * (Options::const_instance()->fastq_data_buffer_size / 8 + 1);
  assert(fastq_data_buffer_size % 2 == 0);

  // Make new buffers of reads from the reads to remap
  std::vector<FastqData> * buf = new std::vector<FastqData>();

  for (int t{0}; t < static_cast<int>(thread_reads_to_remap.size()); ++t)
  {
    std::vector<FastqData> && reads = std::move(thread_reads_to_remap[t]);

    if (reads.size() == 0)
      continue;

    assert(reads.size() % 2 == 0);
    print_info(_HERE_, " buffer has ", reads.size(), " reads.");
    assert(static_cast<int>(buf->size()) < fastq_data_buffer_size);

    if (static_cast<int>(reads.size() + buf->size()) < fastq_data_buffer_size)
    {
      // everything fits into a single buffer, and does not fill it
      buf->insert(buf->end(), //
                  std::make_move_iterator(reads.begin()),
                  std::make_move_iterator(reads.end()));

      assert(buf->size() % 2 == 0);
      continue;
    }

    // then add the remaining buffers
    for (int r{0}; r < static_cast<int>(reads.size()); /*no increment*/)
    {
      int const top_off{fastq_data_buffer_size - static_cast<int>(buf->size())};
      assert(r % 2 == 0);
      assert(top_off % 2 == 0);

      if ((r + top_off) < static_cast<int>(reads.size()))
      {
        buf->insert(buf->end(),
                    std::make_move_iterator(reads.begin() + r),
                    std::make_move_iterator(reads.begin() + r + top_off));

        assert(static_cast<int>(buf->size()) == fastq_data_buffer_size);
      }
      else
      {
        buf->insert(buf->end(), //
                    std::make_move_iterator(reads.begin() + r),
                    std::make_move_iterator(reads.end()));
      }

      new_buffers.push_back(buf);
      buf = new std::vector<FastqData>();
      r += top_off;
    }
  }

  if (buf->size() > 0)
    new_buffers.push_back(buf);
  else
    delete buf;

  return new_buffers;
}

void FastqData::remove_slash_from_name()
{
  size_t const n = name.size();

  if (n > 2 && name[n - 2] == '/')
    name.resize(n - 2);
}

void FastqData::rta3_bin_quals()
{
  for (char & c : qual)
  {
    int const q = c - '!';
    assert(q >= 0);
    assert(q < 64);
    c = rta3_transform_qual[q];
  }
}

} // namespace weaver
