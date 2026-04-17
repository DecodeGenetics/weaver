/*!
 * @file options.cpp
 * @brief Implements the Options class.
 */

#include "options.hpp"

#include <string>
#include <thread>

#include "logging.hpp"

namespace
{
void set_default_threads(weaver::Options & opts)
{
  opts.threads = std::thread::hardware_concurrency();
}

} // namespace

namespace weaver
{
void Options::check_read_group_header_line() const
{
  if (read_group_header_line.size() == 0)
    return;

  if (read_group_header_line.size() <= 8)
  {
    print_error(_HERE_, " Bad read_group_header_line=", read_group_header_line, " is too small.");
    print_error(_HERE_, " Should on the format: @RG\\tID:foo\tSM:bar");
    std::exit(1);
  }

  std::string const check_str = "@RG\\tID:";
  int cmp_val = read_group_header_line.compare(0, 8, check_str, 0, 8);

  if (cmp_val != 0)
  {
    print_error(_HERE_, " Bad read_group_header_line='", read_group_header_line, "' != '", check_str, "'");
    print_error(_HERE_, " returned comparison value=", cmp_val);
    std::exit(1);
  }

  if (read_group_header_line.find("\\tSM:") == std::string::npos)
  {
    print_error(_HERE_, " Bad read_group_header_line='", read_group_header_line, " that has no SM: tag");
    std::exit(1);
  }
}

Options * Options::instance()
{
  return _instance;
}

const Options * Options::const_instance()
{
  return _instance;
}

Options::Options()
{
  ::set_default_threads(*this);
}

Options * Options::_instance = new Options;

} // namespace weaver
