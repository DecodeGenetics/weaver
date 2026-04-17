#pragma once
/*!
 * @file sam_zstd_reader.hpp
 * @brief Defines the SAMZstdReader and ParallelSAMZstdReader classes.
 */

#include <cstdio>  // FILE
#include <string>  // std::string
#include <utility> // std::pair
#include <vector>  // std::vector
#include <zstd.h>  // ZSTD_*

#include "filesystem.hpp"     // filesystem::path
#include "sam_order_line.hpp" // SAMOrderLine

namespace weaver
{
class SAMZstdReader
{
private:
  FILE * fp{nullptr}; // file pointer
  std::size_t const buffer_size{0};
  std::size_t const buffer_output_size{0};
  std::size_t buffer_output_current_pos{0};
  SAMOrderLine order_line{};
  void * buffer{nullptr};
  void * buffer_output{nullptr};

  std::unique_ptr<ZSTD_DCtx, void (*)(ZSTD_DCtx *)> dctx;
  // ZSTD_DCtx * dctx{nullptr};
  ZSTD_outBuffer output{nullptr, 0, 0};
  ZSTD_inBuffer input{nullptr, 0, 0};

public:
  SAMZstdReader() = delete;                                  //!< The reader must have a file_index associated with it.
  explicit SAMZstdReader(int file_index);                    //!< Only allowed construction.
  SAMZstdReader(SAMZstdReader const &) = delete;             //!< No copying should be allowed.
  SAMZstdReader(SAMZstdReader &&) = delete;                  //!< No moving should be allowed.
  SAMZstdReader & operator=(SAMZstdReader const &) = delete; //!< No copying should be allowed.
  SAMZstdReader & operator=(SAMZstdReader &&) = delete;      //!< No moving should be allowed.
  ~SAMZstdReader();                                          //!< Destructor that frees buffers.

  //! Closes the underlying and clears all fields such that a new file can be opened.
  void close();

  //! Opens a file at the given path.
  /*! The provided path must exist and refer to a file. */
  void open(filesystem::path const & sam_zstd_path);

  //! Returns a pointer to the next line in the file. If no line could be read then a nullptr will be returned.
  SAMOrderLine const * read_line();

  //! Returns the file index associated with this reader.
  int get_file_index() const;
};

class ParallelSAMZstdReader
{
private:
  std::vector<SAMOrderLine const *> heap;
  std::vector<std::unique_ptr<SAMZstdReader>> readers;

public:
  ParallelSAMZstdReader() = default;
  explicit ParallelSAMZstdReader(int num_chunks);
  ~ParallelSAMZstdReader() = default;

  //! Adds a SAMZstdReader reader with file_index
  void add_reader();

  //! Open zstd compressed SAM order lined file
  void open(filesystem::path const & sam_zstd_path, int file_index);

  //! Assigns sam_line to the smallest line found in the heap.
  /*! Returns false if no assignment was made, then we are finished processing the files. */
  bool read_next_line(std::string & sam_line);

  //! Assigns sam_line to the smallest sam order line in the heap.
  bool read_next_sam_order_line(SAMOrderLine & sam_order_line);
};

std::vector<filesystem::path> get_sam_zstd_paths(std::vector<std::string> const & sam_lines_fn,
                                                 std::vector<uint32_t> const & chunk_counter);

} // namespace weaver
