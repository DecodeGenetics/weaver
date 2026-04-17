#pragma once
/*!
 * @file in.constants.hpp
 * @brief Global constants, macros and configurations set by CMake.
 */

#include <cstdint>     // int64_t
#include <string.h>    // strrchr
#include <string_view> // std::string_view

// clang-format off
// CMake variables
#define weaver_VERSION_MAJOR @weaver_VERSION_MAJOR@
#define weaver_VERSION_MINOR @weaver_VERSION_MINOR@
#define weaver_VERSION_PATCH @weaver_VERSION_PATCH@
#define weaver_SOURCE_DIRECTORY "@PROJECT_SOURCE_DIR@"
#define weaver_BINARY_DIRECTORY "@PROJECT_BINARY_DIR@"
#define GIT_BRANCH "@GIT_BRANCH@"
#define GIT_COMMIT_SHORT_HASH "@GIT_COMMIT_SHORT_HASH@"
#define GIT_COMMIT_LONG_HASH "@GIT_COMMIT_LONG_HASH@"
#define GIT_NUM_DIRTY_LINES "@GIT_NUM_DIRTY_LINES@"
// clang-format on

namespace weaver
{
std::string_view constexpr weaver_source_dir(weaver_SOURCE_DIRECTORY);
std::string_view constexpr weaver_binary_dir(weaver_BINARY_DIRECTORY);

/*!
 * @brief minimizer/k-mer size
 *
 * @details
 * Memory layout of keys for KMIN=23 :
 *
 *    [ c,current minimizer hash (46 bits) | NP key, previous min-hash XOR next min-hash (18 bits) ]
 */
int constexpr K_DEFAULT{23};

/*!
 * @brief Number of k-mers in each window.
 *
 * @detauls
 * The k-mer with the smallest hash value in any window is a minimizer.
 */
int constexpr W_DEFAULT{11};

//! Maximum value of \a w, modify if larger \a w is needed.
std::size_t constexpr MAX_W{32};

//! How many hit interval should be extracted per read.
uint64_t constexpr NUM_HIT_INTERVALS{24};

namespace SAMFlags
{
uint16_t constexpr IS_PAIRED = 1;             //!< Set when paired.
uint16_t constexpr IS_PROPER_PAIR = 2;        //!< Set when paired and proper.
uint16_t constexpr IS_UNMAPPED = 4;           //!< Set when not mapped.
uint16_t constexpr IS_MATE_UNMAPPED = 8;      //!< Set when mate exists and is not mapped.
uint16_t constexpr IS_SEQ_REVERSED = 16;      //!< Set when sequence has been reversed.
uint16_t constexpr IS_MATE_SEQ_REVERSED = 32; //!< Set when mate sequence exists and has been reversed.
uint16_t constexpr IS_FIRST_IN_PAIR = 64;     //!< Set when read is paired and first in pair.
uint16_t constexpr IS_SECOND_IN_PAIR = 128;   //!< Set when read is paired and second in pair.
uint16_t constexpr IS_SECONDARY = 256;        //!< Set when read is secondary.
uint16_t constexpr IS_QC_FAIL = 512;          //!< Set when read fails on some QC metric.
uint16_t constexpr IS_DUPLICATION = 1024;     //!< Set when read is marked as a duplication.
uint16_t constexpr IS_SUPPLEMENTARY = 2048;   //!< Set when read is marked as supplementary read.

} // namespace SAMFlags

/*!
 * @}
 */
} // namespace weaver

// Macros
#define S1_weaver_internal__(x) #x
#define S2_weaver_internal__(x) S1_weaver_internal__(x)
#define _HERE_ (strrchr("/" __FILE__ ":" S2_weaver_internal__(__LINE__), '/') + 1)
