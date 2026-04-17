/*!
 * @file logging.cpp
 * @brief Implements the Logging class.
 */
#include "logging.hpp"

#include <fstream>
#include <memory>
#include <ostream>
#include <string>

namespace weaver
{
std::unique_ptr<Logging> log_singleton = nullptr;

Logging::Logging(log_severity _severity, std::ostream & _sink) noexcept :
  severity{_severity}, sink{&_sink, stream_deleter_noop}
{
}

Logging::Logging(log_severity _severity, std::string const & filename) :
  severity{_severity}, sink{new std::ofstream{filename, std::ios::binary}, stream_deleter_default}
{
}

void Logging::flush_stream()
{
  assert(sink);

  std::lock_guard guard{mutex};
  sink->flush();
}

} // namespace weaver
