#pragma once

#include <string>
#include <vector>
#include <zlib.h> // gzFile

#include "htslib/kseq.h"

// ignore unused function warnings in kseq.h
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic push
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

namespace weaver
{
class FastqData
{
public:
  /*!
   * @name Public instance variables
   * @{
   */
  std::string name;     //!< Name of the sequence
  std::string comment;  //!< Fastq comments
  std::string seq;      //!< Sequence
  std::string qual;     //!< Quality string for the sequence
  int read_group_id{0}; //!< Read group ID from comment

  /*!
   * @}
   * @name Constructors and destructor
   * @{
   */
  FastqData() = default;
  explicit FastqData(kseq_t const & kseq);

  FastqData(FastqData const &) = default;
  FastqData(FastqData &&) = default;
  FastqData & operator=(FastqData const &) = default;
  FastqData & operator=(FastqData &&) = default;
  ~FastqData() = default;

  /*!
   * @}
   */

  //! Remove the slash and read number ('/{1,2}') from the read name
  void remove_slash_from_name();

  //! Call to bin the values of the quality string
  void rta3_bin_quals();
};

std::vector<FastqData> * seq_store_get(kseq_t * seq);
std::vector<FastqData> * seq_store_get(kseq_t * seq1, kseq_t * seq2);

std::vector<std::vector<FastqData> *> get_new_fastq_buffers(
  std::vector<std::vector<FastqData>> && thread_reads_to_remap);

} // namespace weaver
