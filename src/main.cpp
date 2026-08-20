/*!
 * @file main.cpp
 * @brief Entry point for weaver. Handles argument parsing.
 */
#include <iostream>
#include <memory>

#include <paw/parser.hpp>

#include "check_options.hpp"
#include "filesystem.hpp"
#include "gfa.hpp"
#include "hashmap.hpp"
#include "index.hpp"    // make_index()
#include "index_io.hpp" // serialize_index()
#include "log_base.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "pesr.hpp"

namespace weaver
{
//! Sets up logger based on verbosity flags.
void setup_logger()
{
  Options & opts = *(Options::instance());
  log_severity severity{};

  if (opts.vverbose)
    severity = log_severity::debug;
  else if (opts.verbose)
    severity = log_severity::info;
  else
    severity = log_severity::warning;

  if (opts.log.size() == 0 || opts.log == "-")
    log_singleton = std::make_unique<Logging>(severity, std::clog);
  else
    log_singleton = std::make_unique<Logging>(severity, opts.log);
}

//! Parses the "weaver map" subcommand
int subcmd_map(paw::Parser & parser)
{
  Options & opts = *(Options::instance());

  {
    // Create a command line string
    std::vector<std::string> const & raw_args = parser.get_raw_args_reference();

    if (raw_args.size() > 0)
    {
      opts.command_line += raw_args[0];

      for (int r{1}; r < static_cast<int>(raw_args.size()); ++r)
      {
        opts.command_line += ' ';
        opts.command_line += raw_args[r];
      }
    }
  }

  std::string graph_fn;
  std::string vcf_fn;
  std::string fastq1;
  std::string fastq2;
  int k{-1};
  int w{-1};
  bool see_advanced_options{false};

  parser.parse_positional_argument(graph_fn, "gfa", "Map to this GFA graph.");
  parser.parse_positional_argument(fastq1, "read.fastq[.gz]", "Map reads in this FASTQ file.");

  parser.parse_option(fastq2,
                      '2',
                      "fq2",
                      "If there are separate files for read1 and read2, this specifies the read2 file.");

  // Parse options
  parser.parse_option(opts.read_group_header_line,
                      'R',
                      "read-group-header-line",
                      "The full read group header line"
                      " such as '@RG\\tID:foo\\tSM:bar'.");

  parser.parse_option(opts.extra_header_lines,
                      'H',
                      "extra-header-lines",
                      "Insert extra header lines from STR. "
                      "If STR does not start with '@' then lines from a file with path STR will be inserted instead.");

  parser.parse_option(opts.max_open_files, 'm', "max-open-files", "Max. number of open files allowed when sorting.");

  parser.parse_option(opts.threads,
                      '@',
                      "threads",
                      "Max. number of threads to use. Default value is determined based on current system.");

  parser.parse_option(see_advanced_options,
                      'a',
                      "advanced",
                      "Set to enable advanced options. "
                      "See a list of all options (including advanced) with 'weaver map --advanced --help'");

  if (see_advanced_options)
    parser.see_advanced_options(true);

  parser.parse_advanced_option(opts.no_PG, 'z', "no-PG", "Set to skip writing a weaver PG-line in the header.");

  parser.parse_advanced_option(opts.no_cleanup, ' ', "no-cleanup", "Set to skip cleanup of temporary files.");

  parser.parse_advanced_option(opts.debug_read_name,
                               ' ',
                               "debug-read-name",
                               "(debug builds only) Read name which should only be parsed.");

  parser.parse_advanced_option(k,
                               'k',
                               "k",
                               "Note: Typically this option should be when calling the the 'index' subcommand. "
                               "Size of the minimizer k-mers. "
                               "If <= 0, then use what the index has. If no index, use " +
                                 std::to_string(K_DEFAULT));

  parser.parse_advanced_option(w,
                               'w',
                               "w",
                               "Note: Typically this option should be when calling the the 'index' subcommand. "
                               "Window size to use when minimizer are found. Smaller values make the index largest."
                               " If <= 0, then use what the index has. If no index, use " +
                                 std::to_string(W_DEFAULT));

  parser.parse_advanced_option(vcf_fn,
                               'v',
                               "vcf",
                               "Path to a VCF.gz file containing small variants (only used if no index is found).");

  parser.parse_advanced_option(opts.match, 'A', "match", "Alignment score for a match.");
  parser.parse_advanced_option(opts.mismatch, 'B', "mismatch", "Alignment penalty for a mismatch.");
  parser.parse_advanced_option(opts.gap_open, 'O', "gap-open", "Alignment penalty for a openning a gap.");
  parser.parse_advanced_option(opts.gap_extend, 'E', "gap-extend", "Alignment penalty for a extending a gap.");

  parser.parse_advanced_option(opts.mininum_score_to_output,
                               'T',
                               "minimum-score-to-output",
                               "Minimum alignment score in output. Alignments with a lower score will be unmapped.");

  parser.parse_advanced_option(opts.fastq_data_buffer_size,
                               ' ',
                               "fastq-data-buffer-size",
                               "The number of reads each chunk contains "
                               "in multi-threading mode. The chunks are "
                               "given to different threads.");

  parser.parse_advanced_option(opts.is_no_adapter_removal,
                               ' ',
                               "no-adapter-removal",
                               "Set to disable adapter removal.");

  parser.parse_advanced_option(opts.gc_content, ' ', "gc_content", "GC content of target genome.");
  parser.parse_advanced_option(opts.identity_penalty,
                               ' ',
                               "identity_penalty",
                               "How much low alignment score should penalize mapping quality. Set as 0 to disable.");
  parser.parse_advanced_option(opts.max_read_buffers,
                               ' ',
                               "max-read-buffers",
                               "Stop reading input once this many read buffers have been read (for debugging).");

  parser.parse_advanced_option(opts.region,
                               'r',
                               "region",
                               "Target region. Only output alignments that overlap the region or have a mate "
                               "alignments that overlaps the region.");

  parser.parse_advanced_option(opts.rta3_quals,
                               ' ',
                               "rta3-quals",
                               "Set to bin quality values similar as introduced in the RTA3 software (NovaSeq) "
                               "to reduce storage footprint.");

  parser.finalize();
  setup_logger();

  // check options
  bool are_options_valid{true};
  is_k_option_valid(are_options_valid, k);
  is_w_option_valid(are_options_valid, w);
  is_read_group_header_line_valid(are_options_valid, opts.read_group_header_line);
  are_fastq_options_valid(are_options_valid, fastq1, fastq2);

  if (not are_options_valid)
  {
    print_error("Invalid options given. Check the above message(s) for problems.");
    return 1;
  }

  print_info("Reads will be tagged with the following read group line: '", opts.read_group_header_line, "'");

  if (fastq2.empty())
  {
    print_info("Mapping in interleaved mode.");
    fastq2 = fastq1; // set same as fastq1 to indicate an interleaved file
  }

  // By default, the log base will be calculated (this is recommended)
  if (opts.log_base <= 0.0)
  {
    opts.log_base = calculate_log_base(); // calculate log_base from other alignment options
    print_debug(_HERE_, " Set log_base = ", opts.log_base);
  }

  pesr(graph_fn, vcf_fn, fastq1, fastq2, k, w);
  return 0;
}

//! Parses the "weaver idxstats" subcommand
int subcmd_idxstats(paw::Parser & parser)
{
  std::string index_fn;
  parser.parse_positional_argument(index_fn, "INDEX", "Path to a weaver index.");
  parser.finalize();
  setup_logger();
  return print_index_stats(index_fn);
}

//! Parses the "weaver index" subcommand
int subcmd_index(paw::Parser & parser)
{
  Options & opts = *(Options::instance());

  std::string graph_fn;
  std::string vcf_fn;
  int k{K_DEFAULT};
  int w{W_DEFAULT};

  parser.parse_positional_argument(graph_fn, "GRAPH", "Path to graph.");
  parser.parse_option(k, 'k', "k", "Size of the minimizer k-mers.");

  parser.parse_option(w,
                      'w',
                      "w",
                      "Window size to use when minimizer are found. Smaller values make the index larger. Max is: " +
                        std::to_string(MAX_W));

  parser.parse_option(opts.threads,
                      '@',
                      "threads",
                      "Max. number of threads to use. Default value is determined based on current system.");

  parser.parse_option(opts.min_edit_as_effect,
                      ' ',
                      "min-edit-as-effect",
                      "Minimum effect an edit needs to have such that it is stored in the index.");

  parser.parse_option(vcf_fn, 'v', "vcf", "Path to a VCF.gz file containing small variants.");
  parser.finalize();
  setup_logger();

  GFA gfa(graph_fn);
  MMI mmi;   // minimizer index
  T_icu icu; // icu index

  if (vcf_fn.size() > 0 && not weaver::filesystem::exists(vcf_fn))
  {
    print_error(" Input VCF file does not exist: '", vcf_fn, "'");
    return 1;
  }

  // check options
  bool are_options_valid{true};
  is_k_option_valid(are_options_valid, k);
  is_w_option_valid(are_options_valid, w);

  if (not are_options_valid)
  {
    print_error(" Invalid options given. Check the above message(s) for problems.");
    return 1;
  }

  make_index(gfa, mmi, icu, vcf_fn, k, w);
  serialize_index(graph_fn, mmi, icu, k, w);
  return 0;
}

} // namespace weaver

int main(int argc, char ** argv)
{
#ifdef NDEBUG
  std::ios::sync_with_stdio(false);
#endif

  int ret{0};
  paw::Parser parser(argc, argv);
  parser.set_name("weaver");
  parser.set_version(weaver_VERSION_MAJOR, weaver_VERSION_MINOR, weaver_VERSION_PATCH);

  // uncomment to print command line arguments for any debugging
  // std::copy(argv + 1, argv + argc, std::ostream_iterator<const char *>(std::cerr, "\n"));

  try
  {
    std::string subcmd{};

    parser.add_subcommand("index", "Create and store an index.");
    parser.add_subcommand("map", "Paired-end short read mapping algorithm.");
    parser.add_subcommand("idxstats", "Print the statistics that are stored in a weaver index.");
    parser.add_subcommand("version", "Print version number only.");

    parser.parse_subcommand(subcmd);

    weaver::Options & opts = *(weaver::Options::instance());

    parser.parse_option(opts.log,
                        'l',
                        "log",
                        "Set path to log file. All logger messages will be written to this file. "
                        "Default behaviour is to write logger messages to stderr.");

    parser.parse_option(opts.verbose,
                        'v',
                        "verbose",
                        "Set to output INFO logger messages in addition to the default WARNING and ERROR messages.");

    parser.parse_option(opts.vverbose,
                        ' ',
                        "vverbose",
                        "Set to output DEBUG logger messages in addition to all other types.");

    if (argc == 1)
      parser.throw_help();

    if (subcmd == "idxstats")
    {
      ret = weaver::subcmd_idxstats(parser);
    }
    else if (subcmd == "index")
    {
      ret = weaver::subcmd_index(parser);
    }
    else if (subcmd == "map")
    {
      ret = weaver::subcmd_map(parser);
    }
    else if (subcmd == "version")
    {
      std::cout << weaver_VERSION_MAJOR << "." //
                << weaver_VERSION_MINOR << "." //
                << weaver_VERSION_PATCH;

      std::string_view constexpr git_commit_short_hash(GIT_COMMIT_SHORT_HASH);
      std::string_view constexpr git_commit_long_hash(GIT_COMMIT_LONG_HASH);

      if (git_commit_short_hash.size() > 0)
        std::cout << '-' << git_commit_short_hash << '\n';
      else
        std::cout << '\n';

      if (git_commit_long_hash.size() > 0)
        std::cout << git_commit_long_hash << '\n';
    }
    else if (subcmd.size() == 0)
    {
      parser.finalize();
    }
    else
    {
      ret = 1;
      parser.finalize();
    }
  }
  catch (paw::exception::help const & e)
  {
    std::cout << e.what();
    return 0;
  }
  catch (std::exception const & e)
  {
    std::cerr << e.what();
    return 1;
  }

  return ret;
}
