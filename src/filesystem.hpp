#pragma once
/*!
 * @file filesystem.hpp
 * @brief Includes either the \<filesystem\> or \<experimental/filesystem\> header.
 *
 * Instead of including \<filesystem\>, all other files should include this file and use the weaver::filesystem
 * namespace alias.
 *
 * Source:
 * https://askubuntu.com/questions/1256440/how-to-get-libstdc-with-c17-filesystem-headers-on-ubuntu-18-bionic
 */
#if __has_include(<filesystem>)
#  include <filesystem>
namespace weaver
{
namespace filesystem = std::filesystem;
}
#elif __has_include(<experimental/filesystem>)
#  include <experimental/filesystem>
namespace weaver
{
namespace filesystem = std::experimental::filesystem;
}
#else
#  error "Missing the <filesystem> header."
#endif
