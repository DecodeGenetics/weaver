#pragma once

#include <memory>
#include <string>
#include <variant>

#include "logging.hpp"

#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"

namespace weaver
{
/****************************
 * Type definitions for I/O *
 ****************************/
using file_ptr = std::unique_ptr<FILE, void (*)(FILE *)>;           //!< Type definition for a smart FILE pointer.
using bgzf_ptr = std::unique_ptr<BGZF, void (*)(BGZF *)>;           //!< Type definition for a smart BGZF pointer.
using hts_file_ptr = std::unique_ptr<htsFile, void (*)(htsFile *)>; //!< Type definition for a smart htsFile pointer.
using fastq_variant = std::variant<file_ptr, bgzf_ptr>;             //!< Type definition for a smart FASTQ variant.
using tbx_t_ptr = std::unique_ptr<tbx_t, void (*)(tbx_t *)>;        //!< Type definition for a smart tbx_t pointer.
using hts_itr_t_ptr = std::unique_ptr<hts_itr_t, void (*)(hts_itr_t *)>; //!< Type definition for a hts_itr_t pointer.

//! Compresses data in std::string and writes the contents to path.
void compress_and_write(std::string const & path, std::string const & data);

/************
 * FILE I/O *
 ************/
//! Closes a VCF file stream, i.e. stdout/stdin
inline void close_file_nop(FILE *)
{
  // NOP
}

//! Closes VCF file iff its not null.
inline void close_file(FILE * f)
{
  if (f != nullptr)
  {
    fclose(f);
  }
}

/*!
 * @brief Opens a uncompressed file from filename or get a pointer to stdout/stdin file handles.
 *
 * @details
 * When the filename "-" is used then a pointer to the stdout/stdin file handles will be returned.
 * stdin is opened if filemode is "r", but otherwise stdout is returned.
 *
 * @param[in] fn Filename. Use "-" to open stdout/stdin.
 * @param[in] filemode File mode of the opened file handle.
 *
 * @returns Pointer to an openend file handle.
 */
file_ptr open_file(std::string const & fn, std::string const & filemode);

/************
 * BGZF I/O *
 ************/
//! Closes a bgzf file iff its not null.
inline void close_bgzf(BGZF * bgzf)
{
  if (bgzf != nullptr)
  {
    if (bgzf_close(bgzf) != 0)
    {
      print_error("Failed closing bgzf file.");
      std::exit(1);
    }
  }
}

//! Opens a BGZF file
inline bgzf_ptr open_bgzf(const char * fn, const char * filemode)
{
  bgzf_ptr in_bgzf(bgzf_open(fn, filemode), weaver::close_bgzf);

  if (in_bgzf == nullptr)
  {
    print_error("Failed opening bgzf file ", fn);
    std::exit(1);
  }

  return in_bgzf;
}

/***************
 * htsFile I/O *
 ***************/
inline void close_hts_file(htsFile * f)
{
  if (f != nullptr)
  {
    if (hts_close(f) != 0)
    {
      print_error("Failed closing hts file.");
      std::exit(1);
    }
  }
}

inline hts_file_ptr open_hts_file(const char * fn, const char * fm)
{
  hts_file_ptr ptr(hts_open(fn, fm), weaver::close_hts_file);

  if (ptr == nullptr)
  {
    print_error("Could not open file ", fn);
    std::exit(1);
  }

  return ptr;
}

/*************
 * Tabix I/O *
 *************/
//! Closes a tabix index iff its not null.
inline void close_tbx_t(tbx_t * f)
{
  if (f != nullptr)
    tbx_destroy(f);
}

//! Opens a tabix index
inline tbx_t_ptr open_tbx_t(const char * filename)
{
  tbx_t_ptr ptr(tbx_index_load(filename), weaver::close_tbx_t);

  if (ptr == nullptr)
  {
    print_error("Could not open file ", filename, ".tbi");
    std::exit(1);
  }

  return ptr;
}

/********************
 * hts iterator I/O *
 ********************/
inline void close_hts_itr_t(hts_itr_t * f)
{
  if (f != nullptr)
    tbx_itr_destroy(f);
}

inline hts_itr_t_ptr open_hts_itr_t(tbx_t * tbx, const char * region)
{
  hts_itr_t_ptr ptr(tbx_itr_querys(tbx, region), weaver::close_hts_itr_t);

#ifndef NDEBUG
  if (ptr == nullptr)
    print_debug(_HERE_, " No records found in region ", region);
#endif // NDEBUG

  return ptr;
}

inline hts_itr_t_ptr open_hts_itr_t(tbx_t * tbx, const char * chrom, int begin, int end)
{
  if (tbx == nullptr)
  {
    print_warning("Could not query region, no index to read from.");
    return hts_itr_t_ptr(nullptr, weaver::close_hts_itr_t);
  }

  int tid = tbx_name2id(tbx, chrom);
  hts_itr_t_ptr ptr(tbx_itr_queryi(tbx, tid, begin, end), weaver::close_hts_itr_t);

#ifndef NDEBUG
  if (ptr == nullptr)
    print_debug(_HERE_, " No records found in region ", chrom, " (tid=", tid, "): ", begin, "-", end);
#endif // NDEBUG

  return ptr;
}

/*************
 * FASTQ I/O *
 *************/
//! Open a FASTQ file, which can either be gzip compressed or in raw/uncompressed format.
inline fastq_variant open_fastq(std::string const & fn, std::string const & filemode)
{
  std::size_t const n = fn.size();

  // Assume file is compressed if it has ".gz" or "bgz" ending.
  if (n > 3 && (fn[n - 3] == '.' || fn[n - 3] == 'b') && fn[n - 2] == 'g' && fn[n - 1] == 'z')
    return open_bgzf(fn.c_str(), filemode.c_str()); // open compressed FASTQ
  else
    return open_file(fn, filemode); // open raw/uncompressed FASTQ
}

//! Read a file which is either FASTQ or FASTQ.gz
inline int read_fastq(fastq_variant & fq_var, void * buf, unsigned int len)
{
  if (std::holds_alternative<bgzf_ptr>(fq_var))
  {
    BGZF * ptr = std::get<bgzf_ptr>(fq_var).get();

    if (ptr != nullptr)
      return bgzf_read(ptr, buf, len); // read compressed FASTQ

    return 0;
  }

  FILE * ptr = std::get<file_ptr>(fq_var).get();

  if (ptr != nullptr)
    return fread(buf, 1, len, ptr); // read raw/uncompressed FASTQ

  return 0;
}

//! Close the FASTQ
inline void close_fastq(fastq_variant & fq_var)
{
  if (std::holds_alternative<bgzf_ptr>(fq_var))
    close_bgzf(std::get<bgzf_ptr>(fq_var).get());
  else
    close_file(std::get<file_ptr>(fq_var).get());
}

} // namespace weaver
