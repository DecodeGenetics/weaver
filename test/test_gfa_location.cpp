#include <string> // std::string

#include <weaver/gfa.hpp>            // weaver::GFA
#include <weaver/gfa_location.hpp>   // weaver::GFALocation
#include <weaver/sequence_utils.hpp> // weaver::get_reverse_complement

#include "test.hpp" // read_graph

namespace weaver::test
{
//! Tests for GFALocation::GFALocation()
static void test_gfa_location_constructor()
{
  GFA gfa = read_graph("test_small.gfa");

  {
    GFALocation zero(0);
    REQUIRE(zero.is_valid());
    REQUIRE(zero.rid == 0);
    REQUIRE(zero.pos == 0);
    REQUIRE(zero.strand == 0);

    GFALocation one(1);
    REQUIRE(one.is_valid());
    REQUIRE(one.rid == 0);
    REQUIRE(one.pos == 0);
    REQUIRE(one.strand == 1);

    GFALocation two(2);
    REQUIRE(two.is_valid());
    REQUIRE(two.rid == 0);
    REQUIRE(two.pos == 1);
    REQUIRE(two.strand == 0);
  }
}

//! Tests for the GFALocation::get_base()
static void test_gfa_location_get_base()
{
  GFA gfa = read_graph("test_small.gfa");

  {
    GFALocation test(9ull << 1 | 1ull);
    std::string seq;
    test.get_base(seq);
    REQUIRE(seq == "G");
  }

  // first base
  {
    GFALocation test(0ull << 1 | 0ull);
    std::string seq;
    test.get_base(seq);
    REQUIRE(seq == "A");

    test.flip_strand();
    test.get_base(seq);
    REQUIRE(seq == "AT");
  }

  // second base
  {
    GFALocation test(1ull << 1 | 0ull);
    std::string seq;
    test.get_base(seq);
    REQUIRE(seq == "A");
    test.flip_strand();
    test.get_base(seq);
    REQUIRE(seq == "AT");
  }

  // third base
  {
    GFALocation test(2ull << 1 | 0ull);
    std::string seq;
    test.get_base(seq);
    REQUIRE(seq == "C");
    test.flip_strand();
    test.get_base(seq);
    REQUIRE(seq == "CG");
  }

  // ninth base
  {
    GFALocation test(9ull << 1 | 0ull);
    std::string seq;
    test.get_base(seq);
    REQUIRE(seq == "C");
    test.flip_strand();
    test.get_base(seq);
    REQUIRE(seq == "CG");
  }

  // tenth base
  {
    GFALocation test(10ull << 1 | 0ull);
    std::string seq;
    test.get_base(seq);
    REQUIRE(seq == "A");
    test.flip_strand();
    test.get_base(seq);
    REQUIRE(seq == "AT");
  }
}

//! Tests for GFALocation::is_same_location()
static void test_gfa_location_is_same_location()
{
  GFALocation loc1(5 << 1ull);
  GFALocation loc2(5 << 1ull);

  REQUIRE(loc1.is_same_location(loc2));
  REQUIRE(loc2.is_same_location(loc1));

  // flipped strand is a different location
  loc1.flip_strand();
  REQUIRE(!loc1.is_same_location(loc2));
  REQUIRE(!loc2.is_same_location(loc1));

  // copied location is the same
  GFALocation loc3(loc1.get_value());
  REQUIRE(loc3.is_same_location(loc1));
  REQUIRE(loc1.is_same_location(loc3));
  REQUIRE(!loc3.is_same_location(loc2));
  REQUIRE(!loc2.is_same_location(loc3));

  // advanced location is different
  loc3.pos++;
  REQUIRE(!loc3.is_same_location(loc1));
  REQUIRE(!loc1.is_same_location(loc3));
  REQUIRE(!loc3.is_same_location(loc2));
  REQUIRE(!loc2.is_same_location(loc3));
}

//! Tests for GFALocation::advance_arc()
static void test_gfa_location_advance_arc()
{
  GFA gfa = read_graph("test_small.gfa");

  // advance arc forward
  {
    GFALocation loc(0);
    std::vector<gfa_arc_t const *> arcs;
    loc.advance_when_same_contig(55, arcs); // get an arc
    REQUIRE(arcs.size() == 1);
    REQUIRE(loc.is_valid());

    gfa_arc_t const * arc = arcs[0];
    GFALocation loc2(0);
    loc2.advance_arc(arc);

    REQUIRE(loc2.rid == 1);
    REQUIRE(loc2.pos == 0);
    REQUIRE(loc2.strand == 0);
  }

  // advance arc in reverse
  {
    GFALocation loc(1ull << 32ull | 5ull << 1ull | 1ull);
    std::vector<gfa_arc_t const *> arcs;
    loc.advance_when_same_contig(50, arcs); // get an arc
    REQUIRE(arcs.size() == 1);
    REQUIRE(loc.is_valid());

    gfa_arc_t const * arc = arcs[0];
    GFALocation loc2(1ull << 32ull | 5ull << 1ull | 1ull);
    loc2.advance_arc(arc);

    REQUIRE(loc2.rid == 0);
    REQUIRE(loc2.pos == 49);
    REQUIRE(loc2.strand == 1);
  }
}

//! Tests for GFALocation::flip_strand()
static void test_gfa_location_flip_strand()
{
  GFA gfa = read_graph("test_small.gfa");

  // Check that strand can be flipped
  {
    GFALocation zero(0);
    GFALocation one(1);
    REQUIRE(!zero.is_same_location(one));

    zero.flip_strand();
    REQUIRE(zero.is_same_location(one));
    REQUIRE(zero.strand == 1);

    one.flip_strand();
    REQUIRE(!zero.is_same_location(one));
    REQUIRE(one.strand == 0);
  }
}

//! Tests for GFALocation::advance_on_segment() and GFALocation::advance_on_segment_and_get_sequence()
static void test_gfa_location_advance_on_segment()
{
  GFA gfa = read_graph("test_small.gfa");

  // No way to move
  {
    GFALocation test(1); // rid=0, pos=0, strand=reverse
    test.advance_on_segment(0);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 1);

    test.advance_on_segment(10);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 1);
  }

  // Simple test
  {
    GFALocation test(0); // rid=0, pos=0, strand=forward
    test.advance_on_segment(0);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 0);

    test.advance_on_segment(10);
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 10);
    REQUIRE(test.strand == 0);

    test.advance_on_segment(100);
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 49);
    REQUIRE(test.strand == 0);

    test.advance_on_segment(30);
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 49);
    REQUIRE(test.strand == 0);
  }

  // A more complex in both directions
  {
    GFALocation test(1ull << 32 | 22ull << 1); // rid=0, pos=0, strand=forward
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 22);
    REQUIRE(test.strand == 0);

    REQUIRE(test.advance_on_segment(0) == 0);
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 22);
    REQUIRE(test.strand == 0);

    REQUIRE(test.advance_on_segment(5) == 0);
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 27);
    REQUIRE(test.strand == 0);

    REQUIRE(test.advance_on_segment(31) == 0);
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 58);
    REQUIRE(test.strand == 0);

    test.flip_strand();
    REQUIRE(test.advance_on_segment(60) == 2);
    test.flip_strand();
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 0);

    REQUIRE(test.advance_on_segment(100) == 41);
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 59);
    REQUIRE(test.strand == 0);

    REQUIRE(test.advance_on_segment(30) == 30);
    REQUIRE(test.rid == 1);
    REQUIRE(test.pos == 59);
    REQUIRE(test.strand == 0);
  }

  {
    GFALocation test(1); // rid=0, pos=0, strand=reverse
    std::string seq;
    test.advance_on_segment_and_get_sequence(0, seq);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 0);

    REQUIRE(test.advance_on_segment_and_get_sequence(10, seq) == 10);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 0);
  }

  {
    GFALocation test(0); // rid=0, pos=0, strand=forward
    std::string seq;
    test.advance_on_segment_and_get_sequence(0, seq);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 0);
    REQUIRE(seq.size() == 0);

    REQUIRE(test.advance_on_segment_and_get_sequence(10, seq) == 0);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 10);
    REQUIRE(test.strand == 0);
    REQUIRE(seq.size() == 10);
    REQUIRE(seq == "AACAAGTTCC");

    test.get_base(seq);
    REQUIRE(seq.size() == 11);
    REQUIRE(seq == "AACAAGTTCCA");
  }

  {
    GFALocation test(0); // rid=0, pos=0, strand=forward
    std::string seq;
    long by = test.advance_on_segment_and_get_sequence(50, seq);

    // No changes
    REQUIRE(by == 1);
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 49);
    REQUIRE(test.strand == 0);
    REQUIRE(seq.size() == 49);
    REQUIRE(seq == "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGAT");
  }

  {
    // A|AC|AAGTTCCAGAAGATAGCTA...
    GFALocation test(2ull << 1 | 1ull); // rid=0, pos=0, strand=forward
    std::string seq;
    test.advance_on_segment_and_get_sequence(0, seq);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 2);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 0);

    REQUIRE(test.advance_on_segment_and_get_sequence(1, seq) == 0);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 1);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 1);
    REQUIRE(seq == "G");

    test.get_base(seq);
    REQUIRE(seq.size() == 2);
    REQUIRE(seq == "GT");
  }

  {
    // AACAAGTTCC|AGAAGATAGCTA...
    GFALocation test(9ull << 1 | 1ull); // rid=0, pos=0, strand=forward
    std::string seq;
    test.advance_on_segment_and_get_sequence(0, seq);

    // No changes
    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 9);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 0);

    REQUIRE(test.advance_on_segment_and_get_sequence(9, seq) == 0);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 9);
    REQUIRE(seq == "GGAACTTGT");

    test.get_base(seq);
    REQUIRE(seq.size() == 10);
    REQUIRE(seq == "GGAACTTGTT");
  }
}

//! Tests for GFALocation::advance_until_and_get_sequence()
static void test_gfa_location_advance_until()
{
  GFA gfa = read_graph("test_small.gfa");
  REQUIRE(gfa.get_num_segments() == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  GFALocation::gfa = &gfa;

  {
    // A|AC|AAGTTCCAGAAGATAGCTA...
    GFALocation test(1ull << 1); // rid=0, pos=1, strand=forward
    GFALocation end(2ull << 1);

    std::string seq;
    std::vector<gfa_arc_t const *> arcs; // no arcs
    test.advance_until_and_get_sequence(end, arcs, seq);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 2);
    REQUIRE(test.strand == 0);
    REQUIRE(seq.size() == 2);
    REQUIRE(seq == "AC");
  }

  {
    // A|AC|AAGTTCCAGAAGATAGCTA...
    GFALocation test(2ull << 1 | 1ull); // rid=0, pos=1, strand=forward
    GFALocation end(1ull << 1 | 1ull);

    std::string seq;
    std::vector<gfa_arc_t const *> arcs; // no arcs
    test.advance_until_and_get_sequence(end, arcs, seq);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 1);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 2);
    REQUIRE(seq == "GT");
  }

  {
    // AACAAGTTCC|AGAAGATAGCTA...
    GFALocation test(0); // rid=0, pos=0, strand=forward
    GFALocation end(9ull << 1);

    std::string seq;
    std::vector<gfa_arc_t const *> arcs; // no arcs
    test.advance_until_and_get_sequence(end, arcs, seq);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 9);
    REQUIRE(test.strand == 0);
    REQUIRE(seq.size() == 10);
    REQUIRE(seq == "AACAAGTTCC");
  }

  {
    // .. TAGCTATCTTCT|GGAACTTGTT
    GFALocation test(9ull << 1 | 1ull);
    GFALocation end(0ull | 1ull);

    std::string seq;
    std::vector<gfa_arc_t const *> arcs; // no arcs
    test.advance_until_and_get_sequence(end, arcs, seq);

    REQUIRE(test.rid == 0);
    REQUIRE(test.pos == 0);
    REQUIRE(test.strand == 1);
    REQUIRE(seq.size() == 10);
    REQUIRE(seq == "GGAACTTGTT");
  }
}

//! Tests for GFALocation::advance_when_same_contig() and GFALocation::advance_when_same_contig_and_get_sequence()
static void test_gfa_location_advance_when_same_contig()
{
  GFA gfa = read_graph("test_small.gfa");
  REQUIRE(gfa.get_num_segments() == 6);
  REQUIRE(gfa.get_num_arcs() == 14);

  // zero advance
  {
    std::vector<gfa_arc_t const *> arcs;
    GFALocation loc(0);
    long by = loc.advance_when_same_contig(0, arcs);

    REQUIRE(by == 0);
    REQUIRE(loc.rid == 0);
    REQUIRE(loc.pos == 0);
    REQUIRE(loc.strand == 0);
  }

  // advance inside the segment
  {
    std::vector<gfa_arc_t const *> arcs;
    GFALocation loc(0);
    loc.advance_when_same_contig(49, arcs);

    REQUIRE(loc.rid == 0);
    REQUIRE(loc.pos == 49);
    REQUIRE(loc.strand == 0);
  }

  // advance to end of segment
  {
    std::vector<gfa_arc_t const *> arcs;
    GFALocation loc(0);
    loc.advance_when_same_contig(50, arcs);

    REQUIRE(loc.rid == 1);
    REQUIRE(loc.pos == 0);
    REQUIRE(loc.strand == 0);
    REQUIRE(arcs.size() == 1);

    // flip and go back
    loc.flip_strand();
    loc.advance_when_same_contig(50, arcs);

    REQUIRE(loc.rid == 0);
    REQUIRE(loc.pos == 0);
    REQUIRE(loc.strand == 1);
    REQUIRE(arcs.size() == 2);
  }

  // advance to end of segment
  {
    std::vector<gfa_arc_t const *> arcs;
    GFALocation loc(0);
    std::string seq;
    loc.advance_when_same_contig_and_get_sequence(50, arcs, seq);

    REQUIRE(loc.rid == 1);
    REQUIRE(loc.pos == 0);
    REQUIRE(loc.strand == 0);
    REQUIRE(arcs.size() == 1);
    REQUIRE(seq.size() == 50);
    REQUIRE(seq == "AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATC");

    // flip and go back
    std::string seq2;
    loc.flip_strand();
    loc.advance_when_same_contig(1, arcs);
    loc.advance_when_same_contig_and_get_sequence(50, arcs, seq2);

    REQUIRE(loc.rid == 0);
    REQUIRE(loc.pos == 0);
    REQUIRE(loc.strand == 1);
    REQUIRE(arcs.size() == 2);
    REQUIRE(seq2.size() == 49);

    loc.get_base(seq2);

    REQUIRE(seq2.size() == 50);
    REQUIRE(seq == get_reverse_complement(seq2));
  }
}

//! Tests for GFALocation::get_value()
static void test_gfa_location_get_value()
{
  {
    uint64_t in = 2ull << 1 | 1ull;
    GFALocation in_loc(in);
    REQUIRE(in_loc.get_value() == in);
  }

  {
    uint64_t in = 1ull << 32 | 2ull << 1 | 1ull;
    GFALocation in_loc(in);
    REQUIRE(in_loc.get_value() == in);
  }

  {
    uint64_t in = 3ull << 32 | 5ull << 1;
    GFALocation in_loc(in);
    REQUIRE(in_loc.get_value() == in);
    in_loc.flip_strand();
    REQUIRE(in_loc.get_value() != in);
    REQUIRE(in_loc.get_value() == (in | 1ull));
    in_loc.flip_strand();
    in_loc.pos += 3;
    REQUIRE(in_loc.get_value() == (in + (3ull << 1)));
  }
}

//! \cond TESTS
TEST_CASE("Tests for gfa_location.cpp", "[gfa]")
{
  test_gfa_location_constructor();
  test_gfa_location_get_base();
  test_gfa_location_is_same_location();
  test_gfa_location_advance_arc();
  test_gfa_location_flip_strand();
  test_gfa_location_is_same_location();
  test_gfa_location_advance_on_segment();
  test_gfa_location_advance_until();
  test_gfa_location_advance_when_same_contig();
  test_gfa_location_get_value();
}
//! \endcond
} // namespace weaver::test
