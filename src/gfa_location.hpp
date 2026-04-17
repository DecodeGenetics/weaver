#pragma once
/*!
 * @file gfa_location.hpp
 * @brief Defines the GFALocation class
 */

#include <cstdint> // uint64_t
#include <gfa.h>   // gfa_arc_t
#include <string>  // std::string

#include "gfa.hpp" // GFA

namespace weaver
{
/*!
 * @brief Stores locations on a GFA graph and provides functions to walk in the graph.
 *
 * @headerfile gfa_location.hpp "weaver/gfa_location.hpp"
 *
 * @details
 * GFALocation objects store a static pointer to a GFA graph, which is set on GFA construction.
 * Therefore, if you have a single GFA graph you can assume this pointer will be set to the correct graph.
 * But if you have multiple graphs these objects will assume you are working the most recently opened graph.
 */
class GFALocation
{
public:
  /*!
   * @name Public instance variables
   * @{
   */

  int rid{};     //!< Reference id where the path starts at.
  int pos{};     //!< Position the path starts at.
  bool strand{}; //!< Strand the path starts at. false=forward, true=reverse.

  /*!
   * @}
   *
   * @name Static variables
   * @{
   */

  /*!
   * @brief Static pointer to the graph which the segment is on.
   *
   * @details
   * Pointer is set on GFA construction, initialized as nullptr in gfa_location.cpp
   */
  static GFA const * gfa;

  /*!
   * @}
   *
   * @name Constructors
   * @{
   */

  GFALocation() = delete;                     //!< Empty constructor is deleted.
  explicit GFALocation(uint64_t begin_value); //!< Construct a GFA location.

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  /*!
   * @brief Returns the value represented the location compacted in a single 8-byte (64bit) value.
   *
   * @details
   * This value is also sometimes called "sketch value". The 64 bits have the memory layout:
   *
   *      32 bits |  31 bits  |   1 bit
   *       rid    |    pos    |  strand
   *
   * Creating a GFALocation from the returned value should be the same location.
   *
   * @see SketchCache uses values with the same memory layout.
   * @see test_gfa_location_get_value() tests the functionality of this method.
   */
  uint64_t get_value() const;

  /*!
   * @brief Returns true iff the two locations are the same.
   *
   * @details
   * No checks are made whether the two location are on the same graph. The method compares rid, pos, and strand.
   *
   * @see test_gfa_location_is_same_location() tests the functionality.
   */
  bool is_same_location(GFALocation const & o) const;

  /*!
   * @brief Get the base at location.
   *
   * @details
   * Note that this function only gets the base the gfa location but makes no update to the location. Therefore,
   * calling this function and then another function that gets the sequence will get this base twice.
   *
   * @param[in, out] seq Sequence to add to.
   *
   * @see advance_on_segment_and_get_sequence(int by, std::string & seq)
   * @see advance_when_same_contig_and_get_sequence()
   */
  void get_base(std::string & seq) const;

  //! Returns true iff strand is forward
  bool is_strand_forward() const;

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  //! Call to flip strands, forward strand becomes reverse strand and vice versa.
  void flip_strand();

  /*!
   * @}
   *
   * @name Debugging methods
   * @{
   */

  //! Get a string representation of the gfa location.
  std::string to_string() const;

  /*!
   * @brief Checks if the location is valid, i.e. that it is on the graph.
   *
   * @details
   * For debugging purposes. In debug mode, warning messages will be printed if the location is not valid.
   */
  bool is_valid() const;

  /*!
   * @}
   *
   * @name Methods for advancing in the graph
   * @{
   */

  /*!
   * @brief Advances a given arc.
   *
   * Assertion error will be triggered in debug mode if the arc is invalid or it cannot be followed.
   *
   * @param[in] arc Pointer to the arc to follow/advance.
   *
   * @see test_gfa_location_advance_arc()
   */
  void advance_arc(gfa_arc_t const * arc);

  /*!
   * @brief Advances the location value by \a a many bases but will not change segments.
   *
   * @param[in] by how many base should be advanced.
   *
   * @returns How bases it was unable to advance.
   *
   * @see advance_on_segment_and_get_sequence for also getting the advanced sequence.
   */
  int advance_on_segment(int by);

  /*!
   * @brief Advances the location value by \a by amount and get the underlying sequence stored in the graph,
   *        at the same time.
   *
   * @param[in] by The maximum amount to advance.
   * @param[in] seq The sequence advanced.
   *
   * @returns How many bases it was unable to advance.
   */
  int advance_on_segment_and_get_sequence(int by, std::string & seq);

  /*!
   * @brief Advances whenever we can continue on the same contig.
   *
   * @returns How many bases it was unable to advance.
   *
   * @see advance_when_same_contig_and_get_sequence()
   */
  int advance_when_same_contig(int by, std::vector<gfa_arc_t const *> & new_arcs);

  /*!
   * @brief Advanced whenever we can continue on the same contig and also get the sequence in the graph at the same
   *        time.
   *
   * @see advance_when_same_contig()
   *
   * @returns How many bases it was unable to advance.
   */
  int advance_when_same_contig_and_get_sequence(int by, std::vector<gfa_arc_t const *> & new_arcs, std::string & seq);

  /*!
   * @brief Advanced until a certain position has been seen.
   *
   * @note You need to be sure that the location \a end is reachable from current location.
   *
   * @param[in] end Location to reach.
   * @param[in] arcs which gfa arcs to follow.
   * @param[in,out] seq extracted graph sequence will be appended to this string.
   *
   * @returns How many bases it was unable to advance.
   */
  int advance_until_and_get_sequence(GFALocation const & end,
                                     std::vector<gfa_arc_t const *> const & arcs,
                                     std::string & seq);

  /*!
   * @}
   */
};

} // namespace weaver
