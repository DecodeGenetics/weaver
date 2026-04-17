/*!
 * \file test_sketch.cpp
 * \brief Unit tests for sketching.
 */
#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string> // std::string
#include <vector> // std::vector

#include <weaver/logging.hpp>
#include <weaver/sequence_utils.hpp>
#include <weaver/sketch.hpp>
#include <weaver/sketch_cache.hpp>     // weaver::SketchCache
#include <weaver/sketch_to_string.hpp> // weaver::sketch_to_string
#include <weaver/sketch_value.hpp>     // weaver::sketch_value_rid

#include "test.hpp"

//! \cond TESTS
TEST_CASE("Sketch simple sequences", "[sketch]")
{
  using namespace weaver;

  using u64_pair = std::pair<uint64_t, uint64_t>;
  std::string str;
  int w{5};
  int k{23};
  bool is_reading_forward{true};

  SECTION("poly A")
  {
    SketchCache<u64_pair> cache;
    int rid{1};
    std::vector<u64_pair> sketches;
    str = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    sketch(sketches, cache, str, /*w=*/5, k, rid, is_reading_forward);
    push_last_sketch(sketches, cache, -1, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif // NDEBUG

    REQUIRE(sketches.size() == 7);

    for (auto const & sketch : sketches)
    {
      std::string sketch_key = sketch_key_to_string(sketch.first);
      REQUIRE(sketch_key == "1000011110111011110000000110111100101010010000");

      REQUIRE(sketch_value_rid(sketch.second) == 1);
      REQUIRE(sketch_value_strand(sketch.second) == 1);
    }
  }

  SECTION("Seq 1 with w=10")
  {
    SketchCache<u64_pair> cache;
    int rid{2};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;

    //    0        10 |      20        30
    str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT";
    //     TGCCCTTGCTTAGGATGTATTTA                     11 at A
    //      GCCCTTGCTTAGGATGTATTTAA                    12 at G
    //       CCCTTGCTTAGGATGTATTTAAA                   13 at G <= minimizer
    //        CCTTGCTTAGGATGTATTTAAAA                  14 at A
    //         CTTGCTTAGGATGTATTTAAAAG                 15 at T
    //          TTGCTTAGGATGTATTTAAAAGT                16 <= minimizer in window [14,24[, but if w=11 then there is no
    //           TGCTTAGGATGTATTTAAAAGTT               17    window with this minimizer
    //            GCTTAGGATGTATTTAAAAGTTA              18
    //             CTTAGGATGTATTTAAAAGTTAA             19
    //              TTAGGATGTATTTAAAAGTTAAC            20
    //               TAGGATGTATTTAAAAGTTAACG           21
    //                AGGATGTATTTAAAAGTTAACGT          22
    //                 GGATGTATTTAAAAGTTAACGTA         23
    //                  GATGTATTTAAAAGTTAACGTAC        24 <= minimizer
    //                   ATGTATTTAAAAGTTAACGTACA       25
    //                    TGTATTTAAAAGTTAACGTACAT      26

    sketch(sketches, cache, str, /*w=*/10, k, rid, is_reading_forward);
    push_last_sketch(sketches, cache, -1, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 3);

    {
      auto const & min = sketches.at(0);
      REQUIRE(sketch_key_to_string(min.first) == "0000011100000011001001111111011110101000000011");
      REQUIRE(sketch_value_rid(min.second) == 2);
      REQUIRE(sketch_value_pos(min.second) == 13);
      REQUIRE(sketch_value_strand(min.second) == 1);
    }

    {
      auto const & min = sketches.at(1);
      REQUIRE(sketch_key_to_string(min.first) == "0000100000111000110110001111100010000100001001");
      REQUIRE(sketch_value_rid(min.second) == 2);
      REQUIRE(sketch_value_pos(min.second) == 16);
      REQUIRE(sketch_value_strand(min.second) == 0);
    }

    {
      auto const & min = sketches.at(2);
      REQUIRE(sketch_key_to_string(min.first) == "0000000011010110100111110000000101101001100100");
      REQUIRE(sketch_value_rid(min.second) == 2);
      REQUIRE(sketch_value_pos(min.second) == 24);
      REQUIRE(sketch_value_strand(min.second) == 1);
    }
  }

  SECTION("Seq 1 with w=11")
  {
    // This is the exact same test as the previous one, except has w=11 which makes the minimizer at pos=16 go away
    SketchCache<u64_pair> cache;
    int rid{2};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;

    str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT"; // pos
    //     TGCCCTTGCTTAGGATGTATTTA                     11
    //      GCCCTTGCTTAGGATGTATTTAA                    12
    //       CCCTTGCTTAGGATGTATTTAAA                   13 <= minimizer
    //        CCTTGCTTAGGATGTATTTAAAA                  14
    //         CTTGCTTAGGATGTATTTAAAAG                 15
    //          TTGCTTAGGATGTATTTAAAAGT                16 <= minimizer in window [14,24[, but if w=11 then there is no
    //           TGCTTAGGATGTATTTAAAAGTT               17    window with this minimizer
    //            GCTTAGGATGTATTTAAAAGTTA              18
    //             CTTAGGATGTATTTAAAAGTTAA             19
    //              TTAGGATGTATTTAAAAGTTAAC            20
    //               TAGGATGTATTTAAAAGTTAACG           21
    //                AGGATGTATTTAAAAGTTAACGT          22
    //                 GGATGTATTTAAAAGTTAACGTA         23
    //                  GATGTATTTAAAAGTTAACGTAC        24 <= minimizer
    //                   ATGTATTTAAAAGTTAACGTACA       25
    //                    TGTATTTAAAAGTTAACGTACAT      26

    sketch(sketches, cache, str, /*w=*/11, k, rid, is_reading_forward);

    // push last sketch
    if (cache.min.first != std::numeric_limits<uint64_t>::max())
      sketches.push_back(cache.min);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 2);

    {
      auto const & min = sketches.at(0);
      REQUIRE(sketch_key_to_string(min.first) == "0000011100000011001001111111011110101000000011");
      REQUIRE(sketch_value_rid(min.second) == 2);
      REQUIRE(sketch_value_pos(min.second) == 13);
      REQUIRE(sketch_value_strand(min.second) == 1);
    }

    {
      auto const & min = sketches.at(1);
      REQUIRE(sketch_key_to_string(min.first) == "0000000011010110100111110000000101101001100100");
      REQUIRE(sketch_value_rid(min.second) == 2);
      REQUIRE(sketch_value_pos(min.second) == 24);
      REQUIRE(sketch_value_strand(min.second) == 1);
    }
  }

  SECTION("Arc on same strand")
  {
    SketchCache<u64_pair> cache;
    int rid{3};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;
    str = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 7 sketches in windows in str with rid=3
    std::string str2 = "AAAAAAAAAAAAAA";   // 11 sketches in arc with rid=3
    // in total 18 sketches with rid=3

    sketch(sketches,
           cache,
           str,
           w,
           k,
           rid,
           true, // is_reading_forward,
           -1,   // max_read_size
           0);   // read_init

    int rid2{4};
    sketch(sketches, cache, str2, w, k, rid2,
           is_reading_forward); // str2 and rid2
    push_last_sketch(sketches, cache, -1, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 21);

    for (int s{0}; s < static_cast<int>(sketches.size()); ++s)
    {
      auto const & sketch = sketches[s];
      std::string sketch_key = sketch_key_to_string(sketch.first);
      REQUIRE(sketch_key == "1000011110111011110000000110111100101010010000");

      if (s < 18) // first 18 sketches with rid=3
        REQUIRE(sketch_value_rid(sketch.second) == 3);
      else
        REQUIRE(sketch_value_rid(sketch.second) == 4);

      REQUIRE(sketch_value_strand(sketch.second) == 1);
    }
  }

  SECTION("Seq 1 with arc")
  {
    for (int _i{0}; _i < 2; ++_i)
    {
      std::string str2;
      int local_w{11};

      if (_i == 0)
        str2 = "CAT"; // forward
      else
        str2 = "ATG"; // reverse complement

      int rid{5};
      std::vector<std::pair<uint64_t, uint64_t>> sketches;
      SketchCache<u64_pair> cache;
      str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTA";

      sketch(sketches,
             cache,
             str,
             /*w=*/local_w,
             k,
             rid,
             is_reading_forward,
             -1, // max_read_size
             0); // read_init

      int rid2{6};
      sketch(sketches, cache, str2, /*w=*/local_w, k, rid2, is_reading_forward ^ static_cast<bool>(_i == 1));
      push_last_sketch(sketches, cache, -1, k);

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 2);

      {
        auto const & min = sketches.at(0);
        REQUIRE(sketch_key_to_string(min.first) == "0000011100000011001001111111011110101000000011");
        REQUIRE(sketch_value_rid(min.second) == 5);
        REQUIRE(sketch_value_pos(min.second) == 13);
        REQUIRE(sketch_value_strand(min.second) == 1);
      }

      {
        auto const & min = sketches.at(1);
        REQUIRE(sketch_key_to_string(min.first) == "0000000011010110100111110000000101101001100100");
        REQUIRE(sketch_value_rid(min.second) == 5);
        REQUIRE(sketch_value_pos(min.second) == 24);
        REQUIRE(sketch_value_strand(min.second) == 1);
      }
    }
  }

  SECTION("Seq 1 with arc 2")
  {
    // seq1 TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT

    for (int _i{0}; _i < 2; ++_i)
    {
      int local_w{11};
      SketchCache<u64_pair> cache;
      str = "TGCCCTTGCTTAGGATGTATTTAAAAGTT";
      std::string seq = "AACGTACATGACCATAGGGA";
      std::string str2;

      if (_i == 0)
        str2 = seq;
      else
        str2 = get_reverse_complement(seq);

      int rid = 5;
      std::vector<std::pair<uint64_t, uint64_t>> sketches;

      sketch(sketches,
             cache,
             str,
             /*w=*/local_w,
             k,
             rid,
             true, // is_reading_forward,
             -1,   // max_read_size
             0);   // read_init

      {
        int rid2 = 6;
        std::vector<std::pair<uint64_t, uint64_t>> new_sketches;

        sketch(new_sketches,
               cache,
               str2,
               /*w=*/local_w,
               k,
               rid2,
               is_reading_forward ^ static_cast<bool>(_i == 1),
               -1, // max_read_size
               0); // read_init

        if (cache.min.first != std::numeric_limits<uint64_t>::max())
          new_sketches.push_back(cache.min);

        // push_last_sketch(new_sketches, cache, -1, k);
        std::copy(new_sketches.begin(), new_sketches.end(), std::back_inserter(sketches));
      }

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 3);

      REQUIRE(sketch_key_to_string(sketches.at(0).first) == "0000011100000011001001111111011110101000000011");
      REQUIRE(sketch_value_rid(sketches.at(0).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(0).second) == 13);
      REQUIRE(sketch_value_strand(sketches.at(0).second) == 1);

      REQUIRE(sketch_key_to_string(sketches.at(1).first) == "0000000011010110100111110000000101101001100100");
      REQUIRE(sketch_value_rid(sketches.at(1).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(1).second) == 24);
      REQUIRE(sketch_value_strand(sketches.at(1).second) == 1);

      REQUIRE(sketch_value_rid(sketches.at(2).second) == 6);

      if (_i == 0)
        REQUIRE(sketch_value_pos(sketches.at(2).second) == 1);
      else
        REQUIRE(sketch_value_pos(sketches.at(2).second) == static_cast<int>(str2.size()) - 1 - 1);

      REQUIRE(sketch_value_strand(sketches.at(2).second) == static_cast<unsigned>(_i));
    }
  }

  SECTION("Same as above except moved 1 base between sequences - only changes the positions of the index")
  {
    // seq1 TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT
    for (int _i{0}; _i < 2; ++_i)
    {
      int local_w{11};
      SketchCache<u64_pair> cache;
      str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTA";
      std::string seq = "ACGTACATGACCATAGGGA";
      std::string str2;

      if (_i == 0)
        str2 = seq;
      else
        str2 = get_reverse_complement(seq);

      int rid = 5;
      std::vector<std::pair<uint64_t, uint64_t>> sketches;

      sketch(sketches,
             cache,
             str,
             /*w=*/local_w,
             k,
             rid,
             true, // is_reading_forward,
             -1,   // max_read_size
             0);   // read_init

      {
        int rid2 = 6;
        std::vector<std::pair<uint64_t, uint64_t>> new_sketches;

        sketch(new_sketches,
               cache,
               str2,
               /*w=*/local_w,
               k,
               rid2,
               is_reading_forward ^ static_cast<bool>(_i == 1),
               -1, // max_read_size
               0); // read_init

        push_last_sketch(new_sketches, cache, -1, k);
        std::copy(new_sketches.begin(), new_sketches.end(), std::back_inserter(sketches));
      }

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 3);

      REQUIRE(sketch_key_to_string(sketches.at(0).first) == "0000011100000011001001111111011110101000000011");
      REQUIRE(sketch_value_rid(sketches.at(0).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(0).second) == 13);
      REQUIRE(sketch_value_strand(sketches.at(0).second) == 1);

      REQUIRE(sketch_value_rid(sketches.at(1).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(1).second) == 24);
      REQUIRE(sketch_value_strand(sketches.at(1).second) == 1);

      REQUIRE(sketch_value_rid(sketches.at(2).second) == 6);

      if (_i == 0)
        REQUIRE(sketch_value_pos(sketches.at(2).second) == 0);
      else
        REQUIRE(sketch_value_pos(sketches.at(2).second) == static_cast<int>(str2.size()) - 1);

      REQUIRE(sketch_value_strand(sketches.at(2).second) == static_cast<unsigned>(_i));
    }
  }

  SECTION("(again) Same as above except moved 1 base between sequences - now it has more changes")
  {
    for (int _i{0}; _i < 2; ++_i)
    {
      int local_w{11};
      // seq1 TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT
      SketchCache<u64_pair> cache;
      str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAA";
      REQUIRE(str.size() == 31);
      std::string seq = "CGTACATGACCATAGGGATAGGA";
      std::string str2;

      if (_i == 0)
        str2 = seq;
      else
        str2 = get_reverse_complement(seq);

      int rid = 5;
      std::vector<std::pair<uint64_t, uint64_t>> sketches;

      sketch(sketches,
             cache,
             str,
             local_w,
             k,
             rid,
             true, // is_reading_forward,
             -1,   // max_read_size
             0);   // read_init

      {
        int rid2 = 6;
        std::vector<std::pair<uint64_t, uint64_t>> new_sketches;
        sketch(new_sketches,
               cache,
               str2,
               local_w,
               k,
               rid2,
               is_reading_forward ^ static_cast<bool>(_i == 1),
               -1, // max_read_size
               0); // read_init

        push_last_sketch(new_sketches, cache, -1, k);
        std::copy(new_sketches.begin(), new_sketches.end(), std::back_inserter(sketches));
      }

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 4);

      REQUIRE(sketch_key_to_string(sketches.at(0).first) == "0000011100000011001001111111011110101000000011");
      REQUIRE(sketch_value_rid(sketches.at(0).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(0).second) == 13);
      REQUIRE(sketch_value_strand(sketches.at(0).second) == 1);

      REQUIRE(sketch_key_to_string(sketches.at(1).first) == "0000000011010110100111110000000101101001100100");
      REQUIRE(sketch_value_rid(sketches.at(1).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(1).second) == 24);
      REQUIRE(sketch_value_strand(sketches.at(1).second) == 1);

      REQUIRE(sketch_value_rid(sketches.at(2).second) == 5);
      REQUIRE(sketch_value_pos(sketches.at(2).second) == 30);
      REQUIRE(sketch_value_strand(sketches.at(2).second) == 0);

      REQUIRE(sketch_value_rid(sketches.at(3).second) == 6);

      if (_i == 0)
        REQUIRE(sketch_value_pos(sketches.at(3).second) == 6);
      else
        REQUIRE(sketch_value_pos(sketches.at(3).second) == static_cast<int>(str2.size()) - 6 - 1);

      REQUIRE(sketch_value_strand(sketches.at(3).second) == !_i);
    }
  }
}
//! \endcond

//! \cond TESTS
TEST_CASE("Sketch arcs with max_read_size set", "[sketch]")
{
  using namespace weaver;
  using u64_pair = std::pair<uint64_t, uint64_t>;

  std::string str{};
  int w{5};
  int k{23};
  bool is_reading_forward{true};

  SketchCache<u64_pair> cache;
  int rid{5};
  std::vector<std::pair<uint64_t, uint64_t>> sketches;
  str = "TGCCCTTGCTTAGGATGTATTTAAAAGTT"; // size 29
  sketch(sketches, cache, str, w, k, rid, is_reading_forward);
  int rid2{6};
  std::string str2 = "AACGTACATCATGGAT";

  SECTION("Max read size is 12")
  {
    int max_read_size{12};
    sketch(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 5);
  }

  SECTION("Max read size is 13")
  {
    int max_read_size{13};
    sketch(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 6);
  }

  SECTION("Max read size is 14")
  {
    int max_read_size{14};
    sketch(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 7);
  }

  SECTION("Max read size is -1")
  {
    int max_read_size{-1};
    sketch(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 7);
  }
}
//! \endcond

//! \cond TESTS
TEST_CASE("Sketch same segment twice with arc", "[sketch]")
{
  using namespace weaver;
  using u64_pair = std::pair<uint64_t, uint64_t>;

  SECTION("forward case")
  {
    SketchCache<u64_pair> cache;
    int rid = 8;
    int w = 11;
    int k = 23;
    bool is_reading_forward{true};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;
    std::string str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT";
    sketch(sketches, cache, str, w, k, rid, is_reading_forward);

    REQUIRE(sketches.size() == 1);

    int rid2 = 9;
    std::string str2 = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT";
    int max_read_size{-1};

    sketch(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 6);
  }

  SECTION("reversed case")
  {
    SketchCache<u64_pair> cache;
    int rid = 8;
    int w = 11;
    int k = 23;
    bool is_reading_forward{false};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;
    std::string str = "ATGTACGTTAACTTTTAAATACATCCTAAGCAAGGGCA"; // len 48
    sketch(sketches, cache, str, w, k, rid, is_reading_forward);

    REQUIRE(sketches.size() == 1);

    int rid2 = 9;
    std::string str2 = "ATGTACGTTAACTTTTAAATACATCCTAAGCAAGGGCA";
    int max_read_size{-1};

    sketch(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 6);
  }
}
//! \endcond

//! \cond TESTS
TEST_CASE("Sketch same reversed segment with max_read_size", "[sketch]")
{
  using namespace weaver;
  using u64_pair = std::pair<uint64_t, uint64_t>;

  {
    SketchCache<u64_pair> cache;
    SketchCache<u64_pair> cache_rev;
    int rid = 8;
    int w = 11;
    int k = 23;
    std::vector<std::pair<uint64_t, uint64_t>> sketches;
    std::vector<std::pair<uint64_t, uint64_t>> sketches_rev;
    std::string str =
      "ATGTACGTTAACTTTTAAATACATCCTAAGCAAGGGCAATGTACGTTAACTTTTAA"
      "ATACATCCTAAGCAAGGGCA"; // 76
    std::string str_rev =
      "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACATTGCCCTTGCTTAGG"
      "ATGTATTTAAAAGTTAACGTACAT";

    SECTION("Read everything")
    {
      int max_read_size{-1};

      sketch(sketches, cache, str, w, k, rid, true, max_read_size);
      push_last_sketch(sketches, cache, max_read_size, k);

      sketch(sketches_rev, cache_rev, str_rev, w, k, rid, false, max_read_size);
      push_last_sketch(sketches_rev, cache_rev, max_read_size, k);

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " forward=", sketch_to_string(sketch));

      for (auto const & sketch : sketches_rev)
        print_debug(_HERE_, " reverse=", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 6);
      REQUIRE(sketches.size() == sketches_rev.size());
    }

    SECTION("Max read size k+w")
    {
      int max_read_size{k + w};

      sketch(sketches, cache, str, w, k, rid, true, max_read_size);
      push_last_sketch(sketches, cache, max_read_size, k);

      sketch(sketches_rev, cache_rev, str_rev, w, k, rid, false, max_read_size);
      push_last_sketch(sketches_rev, cache_rev, max_read_size, k);

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " forward=", sketch_to_string(sketch));

      for (auto const & sketch : sketches_rev)
        print_debug(_HERE_, " reverse=", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 1);
      REQUIRE(sketches.size() == sketches_rev.size());
    }

    SECTION("Max read size k+1")
    {
      int max_read_size{k + 1};

      sketch(sketches, cache, str, w, k, rid, true, max_read_size);
      push_last_sketch(sketches, cache, max_read_size, k);

      sketch(sketches_rev, cache_rev, str_rev, w, k, rid, false, max_read_size);
      push_last_sketch(sketches_rev, cache_rev, max_read_size, k);

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " forward=", sketch_to_string(sketch));

      for (auto const & sketch : sketches_rev)
        print_debug(_HERE_, " reverse=", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 1);
      REQUIRE(sketches.size() == sketches_rev.size());
      REQUIRE(sketches[0].first == sketches_rev[0].first);
      REQUIRE(sketch_value_rid(sketches[0].second) == sketch_value_rid(sketches_rev[0].second));
      REQUIRE(sketch_value_pos(sketches[0].second) != sketch_value_pos(sketches_rev[0].second));
      REQUIRE(sketch_value_strand(sketches[0].second) != sketch_value_strand(sketches_rev[0].second));
    }

    SECTION("Max read size k")
    {
      int max_read_size{k};

      sketch(sketches, cache, str, w, k, rid, true, max_read_size);
      push_last_sketch(sketches, cache, max_read_size, k);

      sketch(sketches_rev, cache_rev, str_rev, w, k, rid, false, max_read_size);
      push_last_sketch(sketches_rev, cache_rev, max_read_size, k);

#ifndef NDEBUG
      for (auto const & sketch : sketches)
        print_debug(_HERE_, " forward=", sketch_to_string(sketch));

      for (auto const & sketch : sketches_rev)
        print_debug(_HERE_, " reverse=", sketch_to_string(sketch));
#endif

      REQUIRE(sketches.size() == 0);
      REQUIRE(sketches.size() == sketches_rev.size());
    }
  }
}

TEST_CASE("Sketching boundaries", "[sketch]")
{
  using namespace weaver;
  using u64_pair = std::pair<uint64_t, uint64_t>;

  SECTION("small string is same as normal sketching")
  {
    SketchCache<u64_pair> cache;
    int rid = 8;
    int w = 11;
    int k = 23;
    bool is_reading_forward{true};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;
    std::string str = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT";
    sketch_boundary(sketches, cache, str, w, k, rid, is_reading_forward);

    REQUIRE(sketches.size() == 1);

    int rid2 = 9;
    std::string str2 = "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT";
    int max_read_size{-1};

    sketch_boundary(sketches, cache, str2, w, k, rid2, is_reading_forward, max_read_size);
    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_debug(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 6);
  }

  SECTION("small string is same as normal sketching")
  {
    SketchCache<u64_pair> cache;
    int rid = 8;
    int w = 11;
    int k = 23;
    bool is_reading_forward{true};
    std::vector<std::pair<uint64_t, uint64_t>> sketches;

    std::string str =
      "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACATTGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT"
      "TGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACATTGCCCTTGCTTAGGATGTATTTAAAAGTTAACGTACAT";

    sketch_boundary(sketches, cache, str, w, k, rid, is_reading_forward);

    REQUIRE(sketches.size() == 1);

    {
      std::vector<std::pair<uint64_t, uint64_t>> extra_sketches;
      SketchCache<u64_pair> cache_cp(cache);
      push_last_sketch(extra_sketches, cache_cp, /*max_read_size=*/-1, k);
      REQUIRE(extra_sketches.size() == 1); // for second boundary
    }

    int rid2 = 9;
    std::string str2 = str;
    int max_read_size{-1};

    sketch_boundary(sketches,
                    cache,
                    str2,
                    w,
                    k,
                    rid2,
                    is_reading_forward,
                    max_read_size,
                    /*read_init=*/0,
                    /*add_pos=*/0);

    push_last_sketch(sketches, cache, max_read_size, k);

#ifndef NDEBUG
    for (auto const & sketch : sketches)
      print_info(_HERE_, " ", sketch_to_string(sketch));
#endif

    REQUIRE(sketches.size() == 6);
  }
}
//! \endcond
