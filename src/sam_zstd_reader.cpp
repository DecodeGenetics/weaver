#include "sam_zstd_reader.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <zstd.h>

#include "filesystem.hpp"
#include "logging.hpp"

namespace
{
//! Reassigns a sam order line to a line starting at \a begin and ends at \a newline
void get_sam_order_line(weaver::SAMOrderLine & out, char const * begin, char const * newline)
{
  assert(begin);
  assert(newline);
  assert(begin != newline);

#ifndef NDEBUG
  char const * tab = std::find(begin, newline, '\t');
#endif // NDEBUG
  char * end_ptr;
  assert(tab != newline);

  out.order = strtoull(begin, &end_ptr, 10);

  assert(static_cast<char const *>(end_ptr) == tab);
  char const * c_end_ptr = end_ptr + 1;
  out.sam_line.assign(c_end_ptr, newline);
}

//! Returns the pointer that points to the SAMOrderLine with the greater order.
bool cmp_sam_order_line_gt_ptrs(weaver::SAMOrderLine const * a, weaver::SAMOrderLine const * b)
{
  assert(a != nullptr);
  assert(b != nullptr);

  return a->order > b->order;
}

inline void destroy_dctx(ZSTD_DCtx * dctx)
{
  if (dctx != nullptr)
    ZSTD_freeDCtx(dctx);
}

} // namespace

namespace weaver
{
SAMZstdReader::SAMZstdReader(int file_index) :
  fp(nullptr),
  buffer_size(ZSTD_DStreamInSize()),
  buffer_output_size(ZSTD_DStreamOutSize()),
  order_line(-1),
  dctx(ZSTD_createDCtx(), destroy_dctx)

{
  assert(buffer_size > 0);
  assert(buffer_output_size > 0);
  order_line.file_index = file_index;

  buffer = malloc(buffer_size);
  buffer_output = malloc(buffer_output_size);

  if (buffer == nullptr || buffer_output == nullptr)
  {
    print_error(_HERE_, " Could not allocated memory for buffers.");
    std::exit(1);
  }

  output = {buffer_output, buffer_output_size, 0};
  // dctx = ZSTD_createDCtx();

  if (dctx == nullptr)
  {
    print_error(_HERE_, " Could not create DCtx");
    std::exit(1);
  }

  assert(buffer_size > 0);
  input = {/*src=*/buffer, /*size=*/0, /*pos=*/0};
}

void SAMZstdReader::close()
{
  output.pos = 0;
  assert(output.size == buffer_output_size);
  input.pos = 0;
  input.size = 0;
  buffer_output_current_pos = 0;

  // close the zstd sam file
  if (fp != nullptr)
  {
    assert(feof(fp));
    fclose(fp);
    fp = nullptr;
  }
}

void SAMZstdReader::open(filesystem::path const & sam_zstd_path)
{
  print_debug(_HERE_, " opening ", sam_zstd_path);

  assert(fp == nullptr); // should have been closed
  assert(filesystem::exists(sam_zstd_path));

  fp = std::fopen(sam_zstd_path.c_str(), "rb");
  assert(fp != nullptr);
  std::size_t const bytes_read = fread(buffer, 1, buffer_size, fp);
  assert(bytes_read <= buffer_size);

  if (bytes_read == 0)
  {
    print_warning(_HERE_, " Tried to read empty file: ", sam_zstd_path);
  }
  else
  {
    print_debug(_HERE_, " bytes_read=", bytes_read, " is_eof=", feof(fp));

    // read a new buffer
    input.size = bytes_read;
    input.pos = 0;
    ZSTD_decompressStream(dctx.get(), &output, &input);
  }
}

SAMOrderLine const * SAMZstdReader::read_line()
{
  // if (fp == nullptr)
  //  return nullptr;
  char const * begin = reinterpret_cast<char *>(output.dst) + buffer_output_current_pos;
  char const * end;

  // check if there is some data in output
  if (buffer_output_current_pos < output.pos)
  {
    end = reinterpret_cast<char *>(output.dst) + output.pos;
    char const * newline = std::find(begin, end, '\n');

    // Handle the case when a newline is found in the buffer
    if (newline != end)
    {
      assert(begin != newline);
      get_sam_order_line(order_line, begin, newline);
      buffer_output_current_pos += (newline - begin) + 1u; // +1 because the \n
      return &order_line;
    }
  }
  else
  {
    end = begin;
  }

  // No newline was found, means we need to decompress more from input buffer
  std::string prev_str(begin, end); // first store the data from the previous buffer
  buffer_output_current_pos += (end - begin);
  assert(input.pos <= input.size);

  // Check if there is any more data to decompress from input buffer
  if (input.pos == input.size && not feof(fp))
  {
    // We need to read into the input buffer from the file
    print_debug(_HERE_, " I need a brand new buffer, please.");
    std::size_t bytes_read = fread(buffer, 1, buffer_size, fp);
    print_debug(_HERE_, " bytes_read=", bytes_read, " is_eof=", feof(fp));

    // create a new view of the buffer
    input = {buffer /*src*/, bytes_read /*size*/, 0 /*pos*/};
  }

  if (input.pos < input.size)
  {
    // we have some data yet to decompress in a old buffer
    output.pos = 0;
    buffer_output_current_pos = 0;

    ZSTD_decompressStream(dctx.get(), &output, &input);
    begin = reinterpret_cast<char *>(output.dst);
    end = begin + output.pos;
    char const * newline = std::find(begin, end, '\n');

    assert(newline != end);

    if (prev_str.empty())
    {
      get_sam_order_line(order_line, begin, newline);
    }
    else
    {
      prev_str += std::string(begin, newline);
      get_sam_order_line(order_line, prev_str.c_str(), prev_str.c_str() + prev_str.size());
    }

    buffer_output_current_pos += (newline - begin) + 1u; // +1 because the \n

    return &order_line;
  }

  assert(feof(fp));
  close(); // close the file

  // reached end of file
  print_debug(_HERE_, " reached end of file with file_index = ", order_line.file_index);
  return nullptr;
}

int SAMZstdReader::get_file_index() const
{
  return order_line.file_index;
}

SAMZstdReader::~SAMZstdReader()
{
  // free dctx
  // ZSTD_freeDCtx(dctx);

  // free allocated memory for buffers
  free(buffer);
  free(buffer_output);
}

ParallelSAMZstdReader::ParallelSAMZstdReader(int num_chunks)
{
  // reserve the size we need
  heap.reserve(num_chunks);
  readers.reserve(num_chunks);

  for (int i{0}; i < num_chunks; ++i)
    readers.push_back(std::make_unique<SAMZstdReader>(i));
}

void ParallelSAMZstdReader::open(filesystem::path const & sam_zstd_path, int const file_index)
{
  // open a new reader and give it an index
  assert(file_index < static_cast<int>(readers.size()));

  if (filesystem::exists(sam_zstd_path))
  {
    print_debug(_HERE_, " path exists=", sam_zstd_path);
    assert(readers[file_index] != nullptr);
    SAMZstdReader & reader = *readers[file_index];
    assert(reader.get_file_index() == file_index);
    reader.open(sam_zstd_path);

    // read the first line and add it to the heap
    SAMOrderLine const * sam_line = reader.read_line();

    if (sam_line == nullptr)
    {
      // no lines in file
      print_warning(_HERE_, " I did not read anything from ", sam_zstd_path);
    }
    else
    {
      // push the line to the heap
      heap.push_back(sam_line);
      std::push_heap(heap.begin(), heap.end(), cmp_sam_order_line_gt_ptrs);
      assert(std::is_heap(heap.begin(), heap.end(), cmp_sam_order_line_gt_ptrs));
    }
  }
}

bool ParallelSAMZstdReader::read_next_line(std::string & sam_line)
{
  if (heap.empty())
    return false;

  // make a copy of the minimum sam line
  sam_line = heap[0]->sam_line;
  auto const file_index = heap[0]->file_index;

  assert(file_index < static_cast<int>(readers.size()));
  assert(readers[file_index] != nullptr);

  std::pop_heap(heap.begin(), heap.end(), cmp_sam_order_line_gt_ptrs);
  SAMOrderLine const * new_line = readers[file_index]->read_line();

  if (new_line == nullptr)
  {
    // no more lines in file, remove it from the heap. We do not close the file here, this is done later
    heap.pop_back();
  }
  else
  {
    assert(new_line->file_index == heap.back()->file_index);
    std::push_heap(heap.begin(), heap.end(), cmp_sam_order_line_gt_ptrs);
  }

  // read a new line
  return true;
}

bool ParallelSAMZstdReader::read_next_sam_order_line(SAMOrderLine & sam_order_line)
{
  if (heap.empty())
    return false;

  // make a copy of the minimum sam line
  assert(not heap.empty());
  assert(heap[0] != nullptr);
  sam_order_line = (*heap[0]);
  auto const file_index = heap[0]->file_index;

  assert(file_index < static_cast<int>(readers.size()));
  assert(readers[file_index] != nullptr);

  std::pop_heap(heap.begin(), heap.end(), cmp_sam_order_line_gt_ptrs);
  SAMOrderLine const * new_line = readers[file_index]->read_line();

  if (new_line == nullptr)
  {
    // no more lines in file, remove it from the heap. We do not close the file here, this is done later
    heap.pop_back();
  }
  else
  {
    assert(new_line->file_index == heap.back()->file_index);
    std::push_heap(heap.begin(), heap.end(), cmp_sam_order_line_gt_ptrs);
  }

  return true;
}

std::vector<filesystem::path> get_sam_zstd_paths(std::vector<std::string> const & sam_lines_fn,
                                                 std::vector<uint32_t> const & chunk_counter)
{
  std::vector<filesystem::path> sam_zstd_paths;
  int const threads{static_cast<int>(chunk_counter.size())};

  for (int t{0}; t < threads; ++t)
  {
    int const num_thread_chunks = chunk_counter[t];
    print_debug(_HERE_, " chunk_counter[", t, "] = ", chunk_counter.at(t));
    std::string const & sam_line_fn = sam_lines_fn[t];

    for (int c{0}; c < num_thread_chunks; ++c)
      sam_zstd_paths.push_back(sam_line_fn + std::to_string(c));
  }

  return sam_zstd_paths;
}

} // namespace weaver
