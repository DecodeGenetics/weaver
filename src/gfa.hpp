#pragma once
/*!
 * @file gfa.hpp
 * @brief Defines the GFA class.
 */

#include <gfa.h>       // part of gfatools
#include <string>      // std::string
#include <string_view> // std::string_view
#include <vector>      // std::vector

namespace weaver
{
/*!
 * @brief Read and manipulate a GFA/rGFA file.
 *
 * @headerfile gfa.hpp "weaver/gfa.hpp"
 *
 * @details
 * GFA objects are created from filenames to a GFA or an rGFA file. Default constructor is deleted.
 *
 * The underlying graph object is retrieved using the GFA::get_graph() or GFA::get_graph_ptr() functions.
 * It's type is gfa_t. Memory alignment of gfa_t (from gfatools' gfa.h):
\code{.cpp} typedef struct {
  // segments
  uint32_t m_seg, n_seg, max_rank;
  gfa_seg_t *seg;
  void *h_names;

  // persistent names
  uint32_t m_sseq, n_sseq;
  gfa_sseq_t *sseq;
  void *h_snames;

  // links
  uint64_t m_arc, n_arc;
  gfa_arc_t *arc;
  gfa_aux_t *link_aux;
  uint64_t *idx;
} gfa_t;
\endcode
 *
 * @see gfa_arc.hpp defines further methods for the links/arcs of the GFA.
 * @see gfa_location.hpp defines further methods for
 */
class GFA
{
private:
  gfa_t * g{nullptr}; //!< Underlying graph object from gfatools.

public:
  std::string fn{}; //!< GFA file name.

  /*!
   * @name Constructors and destructors
   * @{
   */

  /*!
   * @brief Default constructor deleted, it is only possible to construct from filename.
   *
   * @see GFA(std::string const & fn) for constructing from filename.
   */
  GFA() = delete;

  /*!
   * @brief Construct from filename.
   *
   * @details
   * The GFA/rGFA graph can be gzipped or in uncompressed GFA format.
   *
   * @see test_gfa_constructor() tests funcionality.
   */
  explicit GFA(std::string const & fn);

  /*!
   * @brief Destructor for GFA.
   *
   * @details
   * Destructing the GFA will destroy the underlying gfa_t object, invalidating all gfa_arc_t pointers,
   * segment sequence views, etc.
   */
  ~GFA();

  /*! @brief Move constructing is allowed.
   *
   * @details
   * This is handled the same as changing the ownership of the underlying graph to the newly constructed object.
   */
  GFA(GFA &&) noexcept;

  //! Copy constructor deleted.
  GFA(GFA const &) = delete;

  /*!
   * @}
   *
   * @name Operators
   * @{
   */

  GFA & operator=(GFA const &) = delete; //!< Copy assignment deleted.
  GFA & operator=(GFA &&) = delete;      //!< Move assignment deleted.

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  /*!
   * @brief Check a single arc if it has the same rank as the segments it links together.
   *
   * @param[in] arc Pointer to an arc to check. Must not be a nullptr.
   *
   * @returns True iff the arc has the same rank as the segments it links together.
   *
   * @see are_arcs_with_same_rank_as_segments() for checking multiple arcs.
   */
  bool is_arc_with_same_rank_as_segments(gfa_arc_t const * arc);

  /*!
   * @brief Check if all arcs are have the same rank as the segments they link together
   *
   * @param[in] arcs Arcs to check.
   *
   * @returns True iff the check was successfull.
   */
  bool are_arcs_with_same_rank_as_segments(std::vector<gfa_arc_t const *> const & arcs) const;

  /*!
   * @brief Get a read-only pointer to the underlying graph object of type gfa_t (from gfa.h).
   *
   * @see get_graph() for getting a graph reference instead.
   */
  gfa_t const * get_graph_ptr() const;

  /*!
   * @brief Get a read-only reference to the underlying graph object of type gfa_t (from gfa.h).
   *
   * @see get_graph_ptr() for getting a graph pointer instead.
   */
  gfa_t const & get_graph() const;

  int get_num_segments() const;                            //!< Get the number of segments in the graph.
  int get_num_stable_segments() const;                     //!< Get the number of stable segments in the graph
  int get_num_arcs() const;                                //!< Get the number of arcs in the graph.
  int get_num_arcs_from_vertex(uint32_t const v) const;    //!< Get the number of arcs from a given vertex ID, v.
  gfa_seg_t const & get_segment(uint32_t const rid) const; //!< Get a reference to a segment of type gfa_seg_t.
  std::string_view view_segment(uint32_t const rid) const; //!< Get a view of a segment sequence.

  /*!
   * @brief Get a reference to the stable sequence.
   *
   * @details
   * The stable sequence is of type gfa_sseq_t, which is a type from gfa.h in gfatools
   *
   * @param[in] segment Segment to get the stable sequence of.
   *
   * @returns Reference to stable sequence struct.
   */
  gfa_sseq_t const & get_stable_sequence(gfa_seg_t const & segment) const;

  /*!
   * @brief Find a arc that continues a stable sequence from vertex that has \a vertex_id .
   *
   * @param[in] vertex_id Vertex ID to search from.
   *
   * @returns Pointer to the arc that continues a stable sequence.
   * @retval nullptr if no such arc exists.
   */
  gfa_arc_t const * get_arc_to_next_stable_sequence(uint32_t vertex_id) const;

  // TODO document
  int64_t get_approximate_stable_position(uint64_t const value) const;

  //! Get the lowest stable position. If rid segment has rank>0 then ignore pos to get the lowest coordinate.
  uint64_t get_lowest_possible_stable_position(uint64_t const value) const;

  //! Order the mmi index values.
  bool mmi_value_order(uint64_t const a, uint64_t const b) const;

  /*!
   * @}
   **/
};

} // namespace weaver
