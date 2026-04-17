#include "sam_zstd_writer.hpp"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <zstd.h>

#include "logging.hpp"
#include "sam_order_line.hpp"

namespace weaver
{
SAMZstdWriter::SAMZstdWriter(std::string const & fn) :
  fp(nullptr), cctx(ZSTD_createCCtx()), buffer_in_size(ZSTD_DStreamInSize()), buffer_out_size(ZSTD_DStreamOutSize())
{
  open(fn);
  assert(fp != nullptr);

  int compression_level = 1;
  ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compression_level);

  assert(buffer_out_size > 0);
  buffer_out = malloc(buffer_out_size);

  str.reserve(buffer_in_size);
}

SAMZstdWriter::~SAMZstdWriter()
{
  close();
  ZSTD_freeCCtx(cctx);
  free(buffer_out);
}

void SAMZstdWriter::close()
{
  if (fp != nullptr)
  {
    flush(/*is_last_chunk=*/true);
    fclose(fp);   // close the file
    fp = nullptr; // makes its nullptr afterwards
  }
}

void SAMZstdWriter::open(std::string const & fn)
{
  assert(fp == nullptr);
  fp = fopen(fn.c_str(), "wb");
}

void SAMZstdWriter::flush(bool const is_last_chunk)
{
  if (str.empty())
    return;

  ZSTD_EndDirective const mode = is_last_chunk ? ZSTD_e_end : ZSTD_e_continue;
  ZSTD_inBuffer input = {str.data(), str.size(), 0};

  bool finished{false};

  while (not finished)
  {
    ZSTD_outBuffer output{buffer_out, buffer_out_size, 0};
    size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    std::size_t written_bytes = fwrite(buffer_out, 1, output.pos, fp);

    if (written_bytes != output.pos)
    {
      print_error(_HERE_, " could not write zstd file");
      print_error(_HERE_, " written_bytes != expected");
      print_error(_HERE_, " ", written_bytes, " != ", output.pos);
      std::exit(1);
    }

    finished = is_last_chunk ? static_cast<bool>(remaining == 0) : static_cast<bool>(input.pos == input.size);
  }

  str.clear();
}

void SAMZstdWriter::write_sam_order_line(SAMOrderLine const & sam_order_line)
{
  // check cache
  if ((str.size() + MAX_SAM_LINE_SIZE) >= buffer_in_size)
    flush(/*is_last_chunk=*/false);

  str += std::to_string(sam_order_line.order);
  str += '\t';
  str += sam_order_line.sam_line;
  str += '\n';
}

void SAMZstdWriter::write_line(std::pair<uint64_t, std::string> const & sam_order_line)
{
  // check cache
  if ((str.size() + MAX_SAM_LINE_SIZE) >= buffer_in_size)
    flush(/*is_last_chunk=*/false);

  str += std::to_string(sam_order_line.first);
  str += '\t';
  str += sam_order_line.second;
  str += '\n';
}

} // namespace weaver
