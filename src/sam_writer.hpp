#pragma once
/*!
 * @file sam_writer.hpp
 * @brief Defines the SAMWriter class.
 */

#include <functional> // std::function
#include <memory>     // std::unique_ptr
#include <ostream>    // std::ostream
#include <string>     // std::string
#include <vector>     // std::vector

namespace weaver
{
class GFA;
class SAMRecord;
class SAMOrderLine;

//! Handles writing SAM/BAM/CRAM files
class SAMWriter
{
public:
  std::string_view static constexpr MISSING_SAM_FIELD = "*";

private:
  //! Pointer with customisable delete behaviour.
  using sam_writer_ptr_t = std::unique_ptr<std::ostream, std::function<void(std::ostream *)>>;

public:
  //! Read group ID string. The ID will be written to each record (if non-empty).
  std::string read_group_id{};

  //! Stream data to sink.
  sam_writer_ptr_t sink;

  //! Open sam filename at "fn".
  explicit SAMWriter(std::string const & fn);

  //! Open sam filename and write to other stream, i.e. std::cout
  explicit SAMWriter(std::ostream & other_stream);

  //! Closes the stream
  /*! @see SAMWriter::open(std::string const & fn) for opening. */
  void close();

  //! Opens the stream. Should be closed before calling open.
  /*! @see SAMWriter::close() for closing.*/
  void open(std::string const & fn);

  //! Flush the stream.
  void flush();

  //! Write the SAM header.
  void write_header(GFA const & gfa);

  //! Open a file and write to stream.
  void write_file(std::string const & fn);

  //! Write a ready SAM line to stream.
  void write_line(std::string const & line);

  //! Write a SAM line to a vector.
  void write_line(std::vector<char> & data, std::string const & line) const;

  //! Write a ready SAM order line
  void write_sam_order_line(SAMOrderLine const & line);

  //! Write a SAM short read record.
  void write_sr_record(SAMRecord const & sam_record,
                       SAMRecord const & other_sam_record,
                       std::string const & cigar,
                       std::string const & other_cigar);
};

//! Generate a complete SAM line in a string from a pair of SAMRecord.
std::string get_sam_string(SAMRecord const & sam_record,
                           SAMRecord const & other_sam_record,
                           std::string const & cigar,
                           std::string const & other_cigar) noexcept;

} // namespace weaver
