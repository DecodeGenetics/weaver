#pragma once

#include <string>

#include "filesystem.hpp"

namespace weaver
{
//! Create a \a directory dir with \a mode set.
void create_dir(std::string const & dir, unsigned mode);
std::string create_temp_dir();
void remove_file(filesystem::path const & path);

} // namespace weaver
