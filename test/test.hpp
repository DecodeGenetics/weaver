#pragma once
/*!
 * \file test.hpp
 * \brief Sets up a logger for the unit tests.
 */

#include <iostream>
#include <memory>
#include <string>

#include "../src/gfa.hpp"
#include "../src/logging.hpp"
#include <catch2/catch.hpp>

namespace weaver::test
{
inline void setup_logging()
{
  if (!weaver::log_singleton)
    weaver::log_singleton = std::make_unique<weaver::Logging>(weaver::log_severity::info, std::clog);
}

/*! \brief This is a global variable during whose initialisation the log_singleton is also initialised */
inline int IGNOREME = []()
{
  setup_logging();
  return 0;
}();

GFA read_graph(std::string const & graph_name);

} // namespace weaver::test
