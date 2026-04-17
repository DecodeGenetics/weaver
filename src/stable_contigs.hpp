#pragma once

#include <cstdint> // uint64_t
#include <string>  // std::string
#include <vector>  // std::vector

#include "hashmap.hpp" // phmap::flat_hash_map

namespace weaver
{
class GFA;

/*!
 * @brief Contig information. Used only internally in StableContigs.
 *
 * @headerfile stable_contigs.hpp "weaver/stable_contigs.hpp"
 */
class Contig
{
public:
  std::string name{}; //!< Contig name, i.e. chr20
  int min{};          //!< Small position on contig, 0-based
  int max{};          //!< Largest position on contig
  int rank{0};        //!< Rank of the contig

  //! Accumulated number of buckets in this contig and all prior contigs.
  int accumulated_num_buckets{0};

  Contig() = default;
  Contig(std::string const & _name, int _min, int _max, int _rank, int _accumulated_num_buckets) :
    name(_name), min(_min), max(_max), rank(_rank), accumulated_num_buckets(_accumulated_num_buckets)
  {
    assert(accumulated_num_buckets > 0);
  }

  //! Get the length of the contig.
  inline int get_length() const
  {
    return max - min;
  }
};

/*!
 * @brief Contains an ordered list of stable contigs.
 *
 * @details
 * The name of the stable contig matches the name of the contig from a FASTA
 * transformed file from the rGFA. This can be done with
 *
 *  gfatools gfa2fa -s <rGFA>
 *
 * Segments are identified with \a rid .
 * The segments that come from the same assembly have the same \a rank .
 * The segments that from the same assembled contig they share \a snid .
 * A formed continuous path on the same \a snid is the stable contig.
 *
 * @headerfile stable_contigs.hpp "weaver/stable_contigs.hpp"
 */
class StableContigs
{
public:
  /*!
   * @name Constructors and destructor
   * @{
   */
  StableContigs() = default;                                 //!< Explicit default construction.
  explicit StableContigs(int const _bucket_size);            //!< Constructor with custom bucket size.
  StableContigs(StableContigs const &) = delete;             //!< Deleted, singleton pattern.
  StableContigs(StableContigs &&) = delete;                  //!< Deleted, singleton pattern.
  StableContigs & operator=(StableContigs const &) = delete; //!< Deleted, singleton pattern.
  StableContigs & operator=(StableContigs &&) = delete;      //!< Deleted, singleton pattern.
  ~StableContigs() = default;                                //!< Explicit default destruction.

  /*!
   * @}
   *
   * @name Public instance variables
   * @{
   */

  //! List of all stable contigs.
  std::vector<Contig> contigs;

  //! Maps contig names to their stable fasta index in \a contigs .
  phmap::flat_hash_map<std::string, int> name2sfa_idx;

  //! Maps snid to their stable fasta index in \a contigs .
  phmap::flat_hash_map<int, std::vector<int>> snid2sfa_idx;

  //! Size of a bucket. Used in sorting.
  int const bucket_size{10'000'000};

  /*!
   * @}
   *
   * @name Modifying methods
   * @{
   */

  /*!
   * @brief Add a contig to the list of stable contigs
   *
   * @param[in] name Contig name of the new contig.
   * @param[in] snid Sequence contig ID.
   * @param[in] min Smallest position on contig.
   * @param[in] max Largest position on contig.
   * @param[in] rank Rank of the contig. Base contig has rank 0.
   *
   * @see StableContigs for exact description of snid.
   */
  void add_contig(std::string const & name, int const snid, int const min, int const max, int const rank);

  //! Clears all instance members of the StableContig class
  void clear();

  /*!
   * @}
   *
   * @name Read-only methods
   * @{
   */

  /*!
   * @brief Get the bucket index that SAM order \a order has.
   *
   * @param[in] order Order to get the bucket index for.
   * @returns The bucket index.
   * @retval get_num_buckets() if RID is missing, i.e. the read is unmapped.
   *
   * @see SAMRecord::get_sam_order() for definition of the SAM order value.
   */
  int get_bucket_index(uint64_t order) const;

  //! Get a read-only reference to a contig with a given name
  Contig const & get_contig(std::string const & contig_name) const;

  //! Get how many buckets there are.
  int get_num_buckets() const;

  //! Get how many stable contigs are stored.
  int get_num_contigs() const;

  //! Retrieve the index to stable fasta contig.
  int get_sfa_idx(int const snid, int const pos) const;

  /*!
   * @brief Retrieve a read-only reference to the contig with a given name.
   *
   * @details
   * If it is not certain the snid exists it must be checked first to avoid undefined behaviour.
   *
   * @param[in] contig_name Contig name
   *
   * @returns reference to the contig
   */
  Contig const & get_contig_reference(std::string const & contig_name) const;

  /*!
   * @}
   */
};

//! Global instance of the stable contigs
extern StableContigs stable_contigs;

/*! @brief Gather all stable contigs from graph and add them to the global instance of StableContigs .
 *
 * @param[in] graph Get stable contigs from this graph.
 */
void set_stable_contigs_from_graph(GFA const & graph);

} // namespace weaver
