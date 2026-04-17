#include "io.hpp"

#include <cassert>
#include <fstream>
#include <string>
#include <zstd.h>

#include "logging.hpp"

namespace weaver
{
void compress_and_write(std::string const & path, std::string const & data)
{
  print_info(_HERE_, " writing to ", path);
  FILE * const fout = fopen(path.c_str(), "wb");
  ZSTD_CCtx * const cctx = ZSTD_createCCtx();
  int compression_level = 1;
  ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compression_level);

  std::size_t const c_buffer_size = ZSTD_compressBound(data.size()); // compressed buffer size
  void * const c_buffer = malloc(c_buffer_size);                     // compressed buffer

  // Compress the data. "comp_size" is the size actually written to the compressed buffer.
  std::size_t const comp_size = ZSTD_compress(c_buffer, c_buffer_size, data.data(), data.size(), 1);
  assert(comp_size <= c_buffer_size);
  std::size_t const written_size = fwrite(c_buffer, 1, comp_size, fout); // Write compressed buffer to file.

  if (written_size != comp_size)
  {
    print_error(_HERE_, " Could not write to ", path);
    std::exit(1);
    // fprintf(stderr, "fwrite: %s : %s \n", fileName, strerror(errno));
    // exit(ERROR_fwrite);
  }

  ZSTD_freeCCtx(cctx);
  fclose(fout);
  free(c_buffer);
}

} // namespace weaver
