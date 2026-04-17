#pragma once

/*!
 * @file sam_zstd_writer.hpp
 * @brief Definition of the SAMZstdWriter class.
 */

#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <zstd.h>

#include "sam_order_line.hpp"

namespace weaver
{
//! Write SAM data into a zstd compressed file.
class SAMZstdWriter
{
public:
  std::size_t static constexpr MAX_SAM_LINE_SIZE = 1024; //!< Each SAM line if garantueed to be at most this large.

  //! File pointer to raw file
  FILE * fp{nullptr};
  ZSTD_CCtx * cctx{nullptr};

  std::size_t buffer_in_size{0};  //!< Size of the zstd input buffer
  std::size_t buffer_out_size{0}; //!< Size of the zstd output buffer

  void * buffer_in{nullptr};  //! zstd input buffer
  void * buffer_out{nullptr}; //! zstd output buffer

  /*!
   * @brief Append data to this string.
   *
   * @details
   * Once the string is almost as large as the input buffer size, set that buffer this string buffer.
   */
  std::string str;

  /***********
   * METHODS *
   ***********/
  SAMZstdWriter() = delete;
  explicit SAMZstdWriter(std::string const & fn);
  ~SAMZstdWriter();

  //! Closes the file
  void close();

  //! Opens a file
  void open(std::string const & fn);

  //! Flushes the file. A hint should be given to zstd wether or not this is the last chunk or not.
  void flush(bool is_last_chunk = true);

  //! Writes an SAM order line
  void write_sam_order_line(SAMOrderLine const & sam_order_line);

  //! Writes a pair containing the order and sam line
  void write_line(std::pair<uint64_t, std::string> const & sam_order_line);
};

} // namespace weaver
