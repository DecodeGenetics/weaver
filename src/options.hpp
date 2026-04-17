#pragma once
/*!
 * @file options.hpp
 * @brief Defines the Options class.
 */

#include <limits>
#include <string>

#include <weaver/constants.hpp>

namespace weaver
{
/*!
 * @brief Singleton object containing user-defined options.
 *
 * @headerfile options.hpp "weaver/options.hpp"
 */
class Options
{
public:
  /*!
   * @name Basic options
   * @{
   */
  bool verbose{false};                    //!< Set to get verbose log messages.
  bool vverbose{false};                   //!< Set to get very verbose log messages.
  std::string log{};                      //!< path to log file. If empty or set to "-" then standard error is used.
  int threads{1};                         //!< The maximum number of threads to use.
  std::string output{"-"};                //!< Output file name. If "-", then output is printed to standard output.
  std::string read_group_header_line{};   //!< The complete read group line for the output header.
  std::string extra_header_lines{};       //!< Insert extra header lines from this string.
  std::string debug_read_name{};          //!< In debug mode, only read with this name will be processed.
  int fastq_data_buffer_size{256 * 1024}; //!< The amount of reads sent to other threads.
  bool no_cleanup{false};                 //!< Set to skip cleanup of temporary files.
  bool is_no_adapter_removal{false};      //!< Set to skip removing adapters (not recommended).
  int max_open_files{8 * 1024};           //!< The maximum number of files allowed to be open at the same time.
  std::string command_line{};             //!< String containing the command line.
  bool no_PG{false};                      //!< Set to skip writing a @PG header line.

  /*!
   * @}
   *
   * @name Indexing options
   * @{
   */
  int max_icu_distance{1500}; //!< Maximum distance between two vertices in the graph for them to be stored in the ICU

  /*!
   * @brief
   * Minimum effect an edit needs to have such that it is stored in the index.
   *
   * @details
   * If the graph contains many (hundreds or more) haplotypes, very rare edits will have practially no effect on
   * the adjusted alignment score and thus can be discard without affecting the results. If an edit would by itself
   * affect the alignment score less then this value it won't be used/stored in the index and discarding it will save
   * both time and memory.
   *
   * Default scoring matrix is assumed when deciding the effect of the edit.
   *
   * Note that even though an edit is discarded this way, the edit will still be used when minimizers are extracted so
   * these rare variants should not be discarded from the input VCF.
   */
  double min_edit_as_effect{0.02};

  /*!
   * @}
   *
   * @name Alignment options
   * @{
   */
  int match{1};                    //!< Score for an alignment match. Expected to be greater than zero.
  int mismatch{4};                 //!< Penalty for an alignment mismatch. Expected to be greater than zero.
  int gap_open{7};                 //!< Penalty for an alignment gap open. Expected to be greater than "gap_extend"
  int gap_extend{1};               //!< Penalty for an alignment gap extend. Expected to be greater than zero.
  int clip{6};                     //!< Penalty for clipping an alignment. Expected to be greater than zero.
  double gc_content{0.41};         //!< Approximate GC content.
  double log_base{-1.0};           //!< Log base value, if log_base<=0 it is calculated score values and gc content.
  double identity_penalty{0.0};    //!< How much low alignment score should penalize mapping quality. 0 to disable.
  int mininum_score_to_output{30}; //!< Alignment with a lower alignment score will be changed to unmapped alignments.
  int max_read_buffers{std::numeric_limits<int>::max()}; //! Maximum number of buffers to read from input FASTQ.
  bool rta3_quals{false}; //!< Set to bin quality values similar to RTA3 to reduce storage.

  /*!
   * @}
   *
   * @name Region options
   * @{
   */
  //! Target region to output reads
  std::string region{};

  //! Target region sfa_idx, -1 if there is no target region.
  int region_sfa_idx{-1};

  //! Target region rank
  int region_rank{0};

  //! Inclusive lower bound of the order to be inside target region
  int region_lower_pos{0};

  //! Inclusive upper bound of the order to be inside target region
  int region_upper_pos{std::numeric_limits<int>::max()};

  /*!
   * @}
   *
   * @name Methods
   * @{
   */
  void check_read_group_header_line() const; //!< Check if the provided read group header line is ok

  static Options * instance();             //!< Gets the only instance (singleton pattern)
  static const Options * const_instance(); //!< The read-only version of Options::instance()

private:
  Options();                                     //!< Prevent construction of new instances
  Options(Options const &) = delete;             //!< Prevent copy-construction
  Options & operator=(Options const &) = delete; //!< Prevent copy-assignment

  static Options * _instance; //!< Singleton instance
};

} // namespace weaver
