/*!
 * \file test.cpp
 * \brief Compiles the Catch unit testing library.
 */

//! \cond TESTS
#define CATCH_CONFIG_MAIN
//! \endcond

#include "test.hpp"

#include <weaver/constants.hpp>
#include <weaver/filesystem.hpp>
#include <weaver/gfa.hpp>

#include <catch2/catch.hpp>

namespace weaver::test
{
GFA read_graph(std::string const & graph_name)
{
  std::string graph_path = std::string(weaver_SOURCE_DIRECTORY) + std::string("/test/data/") + graph_name;

  if (!weaver::filesystem::exists(graph_path))
    print_warning(_HERE_, " No such graph file ", graph_path);

  GFA graph(graph_path);
  return graph;
}
} // namespace weaver::test
