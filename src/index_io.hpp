#pragma once
/*!
 * @file index_io.hpp
 * @brief Declares the function for the weaver index I/O operations.
 */

#include <string> // std::string

#include "filesystem.hpp" // filesystem::path
#include "gfa.hpp"        // GFA
#include "icu.hpp"        // T_icu

namespace weaver
{
class MMI;

/*!
 * @brief Deserialize minimizer and ICU index from disk.
 *
 * The serialized index object is expected at location \a graph_path + ".wmi". The stored object also stores the size of
 * the graph when it was built. If the current size of the graph is not the same as the stored size an error will be
 * thrown.
 *
 * @param[in]  graph_path Filaname path to the GFA graph.
 * @param[out] mmi The minimizer index.
 * @param[out] icu The "I see you" index.
 * @param[in]  k kmer size wanted or -1 to use the index k. Deserialization will abort early if the index has a
 *               different kmer size.
 * @param[in]  w minimizer window size or -1 to use the index w.
 *
 * @returns true iff deserialization finished successfully.
 *
 * @see deserialize_or_make_index() for creating an index if it doesn't exist.
 * @see serialize_index() for serializing the index.
 * @see MMI the minimizer index type description.
 * @see T_icu the ICU index type description.
 */
bool deserialize_index(filesystem::path const & graph_path, MMI & mmi, T_icu & icu, int & k, int & w);

/*!
 * @brief Deserilizes if a serialized object exists, otherwise make the indexes instead.
 *
 * @param[in]  graph_path Filaname path to the GFA graph.
 * @param[in]  gfa GFA graph.
 * @param[out] mmi The minimizer index.
 * @param[out] icu The "I see you" index.
 * @param[in]  vcf_fn Filename path to a VCF file. Empty if no VCF.
 * @param[in]  k minimizer index kmer size wanted. Use -1 to use the stored k. The index will be recreated if it has a
 *               different kmer size.
 * @param[in]  w minimizer index window size. Use -1 to use the stored w.
 *
 * @see deserialize_index() for deserializing an index without making one if it doesn't exist.
 * @see make_index() for making the index without checking if another exists already.
 * @see serialize_index() for serializing the index.
 * @see MMI the minimizer index type description.
 * @see T_icu the ICU index type description.
 */
void deserialize_or_make_index(filesystem::path const & graph_path, //
                               GFA const & gfa,
                               MMI & mmi,
                               T_icu & icu,
                               std::string const & vcf_fn,
                               int & k,
                               int & w);

/*!
 * @brief Serialize a minimizer and ICU index to disk.
 *
 * @param[in] graph_fn Filaname path for writing the graph to.
 * @param[in] mmi The minimizer index.
 * @param[in] icu The "I see you" index.
 * @param[in] k minimizer index kmer size.
 * @param[in] w minimizer index window size.
 *
 * @see make_index() for making the index before it will be serialized.
 */
void serialize_index(filesystem::path const & graph_fn, MMI const & mmi, T_icu const & icu, int k, int w);

//! Print the statistics that are stored in a weaver index at \a index_path
int print_index_stats(filesystem::path const & index_path);
} // namespace weaver
