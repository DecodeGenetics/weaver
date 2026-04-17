#include <cstdint>
#include <sstream>
#include <string>

#include <parallel_hashmap/phmap.h>

#include <weaver/alignment_utils.hpp>
#include <weaver/fastq_data.hpp>
#include <weaver/haplotype_stats.hpp>
#include <weaver/log_base.hpp>
#include <weaver/logging.hpp>
#include <weaver/sequence_utils.hpp>

#include "test.hpp"

#include <weaver.hpp>

namespace
{
weaver::FastqData make_dummy_fastq_data(std::string const & seq)
{
  weaver::FastqData fastq_data;
  fastq_data.name = "test";
  fastq_data.seq = seq;
  fastq_data.qual = std::string(seq.size(), static_cast<char>(33 + 30));
  return fastq_data;
}

} // namespace

namespace weaver::test
{
//! Tests where reads should align to the graph without any softclips/multi-mappings.
static void test_align_simple()
{
  print_info(_HERE_, " test_align_simple..");
  GFA gfa = read_graph("test_small.gfa");
  uint32_t const num_segments = gfa.get_num_segments();
  REQUIRE(num_segments == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  weaver::Options & opts = *(weaver::Options::instance());
  opts.match = 1;
  opts.mismatch = 4;
  opts.gap_open = 6;
  opts.gap_extend = 1;
  opts.clip = 5;
  opts.log_base = weaver::calculate_log_base();

  int k{17};
  int w{3};
  MMI mmi;
  T_icu icu;
  make_index(gfa, mmi, icu, "", k, w);

  HaplotypeStats hap_stats(&mmi.haplotypes);

  // placeholders
  SRSeed empty_seed;
  SRAlignment empty_alignment;

  std::string const segment1 = "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC";

  // Trivial case, both reads are identical to segment1
  {
    FastqData seq1 = make_dummy_fastq_data(segment1);
    FastqData seq2 = make_dummy_fastq_data(get_reverse_complement(segment1));
    REQUIRE(seq1.seq.size() == 50);
    REQUIRE(seq2.seq.size() == 50);

    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    print_info(_HERE_, " ", sam_data_lines[0]);
    print_info(_HERE_, " ", sam_data_lines[1]);

    auto sam0_fields = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam0_fields.size() == 17);
    REQUIRE(sam0_fields[0] == "test");  // qname
    REQUIRE(sam0_fields[1] == "99");    // flag
    REQUIRE(sam0_fields[2] == "chr20"); // rname
    REQUIRE(sam0_fields[3] == "1");     // pos
    REQUIRE(sam0_fields[4] == "60");    // mapq
    REQUIRE(sam0_fields[5] == "50M");   // cigar

    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t1\t60\t50M\tchr20\t1\t50\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC\t"
            "??????????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:50\tWS:i:50\tMQ:i:60\tMC:Z:50M\tms:i:1500");

    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t1\t60\t50M\tchr20\t1\t-50\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC\t"
            "??????????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:50\tWS:i:50\tMQ:i:60\tMC:Z:50M\tms:i:1500");
  }

  // s1 has 1bp deletion
  {
    std::string s1 = segment1.substr(0, 39) + segment1.substr(40);
    std::string const & s2 = get_reverse_complement(segment1);

    FastqData seq1 = make_dummy_fastq_data(s1);
    FastqData seq2 = make_dummy_fastq_data(s2);

    REQUIRE(seq1.seq.size() == 49);
    REQUIRE(seq2.seq.size() == 50);
    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);
    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t1\t60\t39M1D10M\tchr20\t1\t50\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAGAGCAGATC\t"
            "?????????????????????????????????????????????????\t"
            "NM:i:1\tAS:i:43\tWS:i:43\tMQ:i:60\tMC:Z:50M\tms:i:1500");
    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t1\t60\t50M\tchr20\t1\t-50\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC\t"
            "??????????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:50\tWS:i:50\tMQ:i:60\tMC:Z:39M1D10M\tms:i:1470");
  }

  // s1 has 1bp insertion
  // s2 has 3bp deletion
  {
    std::string s1 = segment1.substr(0, 39) + segment1.substr(38);
    std::string const s2 = get_reverse_complement(segment1.substr(0, 39) + segment1.substr(42));

    FastqData seq1 = make_dummy_fastq_data(s1);
    FastqData seq2 = make_dummy_fastq_data(s2);

    REQUIRE(seq1.seq.size() == 51);
    REQUIRE(seq2.seq.size() == 47);
    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);
    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t1\t60\t38M1I12M\tchr20\t1\t50\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGGAAGAGCAGATC\t"
            "???????????????????????????????????????????????????\t"
            "NM:i:1\tAS:i:44\tWS:i:44\tMQ:i:60\tMC:Z:38M3D9M\tms:i:1410");
    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t1\t60\t38M3D9M\tchr20\t1\t-50\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAGCAGATC\t"
            "???????????????????????????????????????????????\t"
            "NM:i:3\tAS:i:39\tWS:i:39\tMQ:i:60\tMC:Z:38M1I12M\tms:i:1530");
  }

  print_info(_HERE_, " test_align_simple DONE");
}

//! Tests for the alignment should have some soft clipping
static void test_align_softclip()
{
  print_info(_HERE_, " test_align_softclip..");
  GFA gfa = read_graph("test_small.gfa");
  uint32_t const num_segments = gfa.get_num_segments();
  REQUIRE(num_segments == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  weaver::Options & opts = *(weaver::Options::instance());
  opts.match = 1;
  opts.mismatch = 4;
  opts.gap_open = 6;
  opts.gap_extend = 1;
  opts.clip = 5;
  opts.mininum_score_to_output = 0;

  int k{17};
  int w{3};
  MMI mmi;
  T_icu icu;
  make_index(gfa, mmi, icu, "", k, w);

  HaplotypeStats hap_stats(&mmi.haplotypes);

  // placeholders
  SRSeed empty_seed;
  SRAlignment empty_alignment;

  std::string const s1 = "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC";

  // forward read1 is soft clipped at back
  {
    // AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCA|CATGAAGAGCAGATC" clip should be where the "|" is
    FastqData seq1 = make_dummy_fastq_data(s1.substr(0, 35) + std::string(10, 'A'));
    FastqData seq2 = make_dummy_fastq_data(get_reverse_complement(seq1.seq));
    // REQUIRE(seq1.seq.size() == 45);
    // REQUIRE(seq2.seq.size() == 45);

    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);
      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    auto sam_fields0 = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam_fields0[5] == "35M10S");
    auto sam_fields1 = split_string(sam_data_lines[1], '\t');
    REQUIRE(sam_fields1[5] == "35M10S");

    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t1\t60\t35M10S\tchr20\t1\t35\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "?????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:30\tWS:i:30\tMQ:i:60\tMC:Z:35M10S\tms:i:1350");
    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t1\t60\t35M10S\tchr20\t1\t-35\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "?????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:30\tWS:i:30\tMQ:i:60\tMC:Z:35M10S\tms:i:1350");
  }

  // forward read1 is soft clipped at front
  {
    // AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCA|CATGAAGAGCAGATC" clip should be where the "|" is
    FastqData seq1 = make_dummy_fastq_data(std::string(10, 'A') + s1.substr(12, 30));
    FastqData seq2 = make_dummy_fastq_data(get_reverse_complement(seq1.seq));
    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);
      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    auto sam_fields0 = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam_fields0[5] == "10S30M");
    auto sam_fields1 = split_string(sam_data_lines[1], '\t');
    REQUIRE(sam_fields1[5] == "10S30M");

    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t13\t60\t10S30M\tchr20\t13\t30\t"
            "AAAAAAAAAAAAGATAGCTAGAGGATGGGAGCACATGAAG\t"
            "????????????????????????????????????????\t"
            "NM:i:0\tAS:i:25\tWS:i:25\tMQ:i:60\tMC:Z:10S30M\tms:i:1200");

    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t13\t60\t10S30M\tchr20\t13\t-30\t"
            "AAAAAAAAAAAAGATAGCTAGAGGATGGGAGCACATGAAG\t"
            "????????????????????????????????????????\t"
            "NM:i:0\tAS:i:25\tWS:i:25\tMQ:i:60\tMC:Z:10S30M\tms:i:1200");
  }

  // forward read1 is soft clipped at back, and doesn't start at position 1
  {
    // AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCA|CATGAAGAGCAGATC" clip should be where the "|" is
    FastqData seq1 = make_dummy_fastq_data(s1.substr(2, 33) + std::string(10, 'A'));
    FastqData seq2 = make_dummy_fastq_data(get_reverse_complement(seq1.seq));
    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    auto sam_fields0 = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam_fields0[5] == "33M10S");
    auto sam_fields1 = split_string(sam_data_lines[1], '\t');
    REQUIRE(sam_fields1[5] == "33M10S");

    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t3\t60\t33M10S\tchr20\t3\t33\t"
            "CAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "???????????????????????????????????????????\t"
            "NM:i:0\tAS:i:28\tWS:i:28\tMQ:i:60\tMC:Z:33M10S\tms:i:1290");
    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t3\t60\t33M10S\tchr20\t3\t-33\t"
            "CAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "???????????????????????????????????????????\t"
            "NM:i:0\tAS:i:28\tWS:i:28\tMQ:i:60\tMC:Z:33M10S\tms:i:1290");
  }

  // reverse read1 is soft clipped at back, and doesn't start at position 1
  {
    // AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCA|CATGAAGAGCAGATC" clip should be where the "|" is
    std::string const seq = s1.substr(2, 33) + std::string(10, 'A');
    FastqData seq1 = make_dummy_fastq_data(get_reverse_complement(seq));
    FastqData seq2 = make_dummy_fastq_data(seq);
    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    auto sam_fields0 = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam_fields0[5] == "33M10S");
    auto sam_fields1 = split_string(sam_data_lines[1], '\t');
    REQUIRE(sam_fields1[5] == "33M10S");

    REQUIRE(sam_data_lines[0] ==
            "test\t83\tchr20\t3\t60\t33M10S\tchr20\t3\t33\t"
            "CAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "???????????????????????????????????????????\t"
            "NM:i:0\tAS:i:28\tWS:i:28\tMQ:i:60\tMC:Z:33M10S\tms:i:1290");

    REQUIRE(sam_data_lines[1] ==
            "test\t163\tchr20\t3\t60\t33M10S\tchr20\t3\t-33\t"
            "CAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "???????????????????????????????????????????\t"
            "NM:i:0\tAS:i:28\tWS:i:28\tMQ:i:60\tMC:Z:33M10S\tms:i:1290");
  }

  // test adapter removal
  {
    print_info(_HERE_, " adapter removal test");
    // AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCA|CATGAAGAGCAGATC" clip should be where the "|" is
    FastqData seq1 = make_dummy_fastq_data(s1.substr(12, 33));
    FastqData seq2 = make_dummy_fastq_data(get_reverse_complement(s1.substr(9, 33)));

    // 3 bp should be removed as adapters
    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);

      print_info(_HERE_, " pos=", record1.pos, " cig=", cigar2string(record1.cig.begin(), record1.cig.end()));
      print_info(_HERE_, " pos=", record2.pos, " cig=", cigar2string(record2.cig.begin(), record2.cig.end()));
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    auto sam_fields0 = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam_fields0[2] == "chr20"); // rname
    REQUIRE(sam_fields0[3] == "13");    // pos
    REQUIRE(sam_fields0[5] == "30M3S"); // cigar

    auto sam_fields1 = split_string(sam_data_lines[1], '\t');
    REQUIRE(sam_fields1[2] == "chr20"); // rname
    REQUIRE(sam_fields1[3] == "13");    // pos
    REQUIRE(sam_fields1[5] == "3S30M"); // cigar
  }

  print_info(_HERE_, " test_align_softclip DONE");
}

//! Tests for the alignment on the full graph
static void test_align_full_graph()
{
  print_info(_HERE_, " test_align_full_graph..");
  GFA gfa = read_graph("test_human_10k.gfa.gz");
  uint32_t const num_segments = gfa.get_num_segments();
  REQUIRE(num_segments == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  weaver::Options & opts = *(weaver::Options::instance());
  opts.match = 1;
  opts.mismatch = 4;
  opts.gap_open = 6;
  opts.gap_extend = 1;
  opts.clip = 5;

  int k{17};
  int w{3};
  MMI mmi;
  T_icu icu;
  make_index(gfa, mmi, icu, "", k, w);

  HaplotypeStats hap_stats(&mmi.haplotypes);

  // placeholders
  SRSeed empty_seed;
  SRAlignment empty_alignment;

  std::string const s1 =
    "CCCAGCTACTTGGGAAGCTGAGGCAGGAGAATCGCTTGAATCCAGGAGGCAGAGGTTGCCATGAGCCGAGATTGTGCCAC"
    "NGCATTCCAGCCTGGGCGATAGAGTGAGACTCCGCCTCAAAATAAATAAATAAATAGATAAATAAATAAAT";

  std::string const s2 =
    "GCATGTCCAAAAATAGCACCAAATCAAATAGTTTTAATCCTCAATTTTGTTAAATAGCAGGATATATTGGATATTTTTCC"
    "AAAATAAATCTAAATAAATATTGCAAAACTAATAAATATACCAACATGGTATACTTTAGAACTATGAGTGA";

  {
    FastqData seq1 = make_dummy_fastq_data(s1);
    FastqData seq2 = make_dummy_fastq_data(s2);

    std::ostringstream ss;

    // scope around the sam writer
    {
      SAMWriter writer(ss);
      SAMRecord record1(seq1);
      SAMRecord record2(seq2);

      sr_align_pair(/*is_unique_mapping=*/false,
                    gfa,
                    mmi,
                    icu,
                    hap_stats,
                    record1,
                    record2,
                    k,
                    w,
                    seq1.seq,
                    seq2.seq,
                    /*is_debug=*/true);

      sam_write(writer, record1, record2);
    }

    std::string sam_data = ss.str();
    auto sam_data_lines = split_string(sam_data, '\n');

    REQUIRE(sam_data_lines.size() == 2);

    /*
    auto sam_fields0 = split_string(sam_data_lines[0], '\t');
    REQUIRE(sam_fields0[5] == "35M10S");
    auto sam_fields1 = split_string(sam_data_lines[1], '\t');
    REQUIRE(sam_fields1[5] == "35M10S");

    REQUIRE(sam_data_lines[0] ==
            "test\t99\tchr20\t1\t60\t35M10S\tchr20\t1\t35\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "?????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:35\tMQ:i:60\tMC:Z:35M10S\tms:i:1350");

    REQUIRE(sam_data_lines[1] ==
            "test\t147\tchr20\t1\t60\t35M10S\tchr20\t1\t-35\t"
            "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCAAAAAAAAAAA\t"
            "?????????????????????????????????????????????\t"
            "NM:i:0\tAS:i:35\tMQ:i:60\tMC:Z:35M10S\tms:i:1350");
    */
  }

  print_info(_HERE_, " test_align_full_graph DONE");
}

//! \cond TESTS
TEST_CASE("Tests for sr_align.cpp.", "[sam]")
{
  test_align_simple();
  test_align_softclip();
  test_align_full_graph();
}
//! \endcond
} // namespace weaver::test
