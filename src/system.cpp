#include "system.hpp"

#include <algorithm>
#include <cassert>
#include <dirent.h> // DIR
#include <random>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

#include "filesystem.hpp"
#include "logging.hpp"

namespace
{
std::string get_env_var(std::string const & key, std::string const & _default = "")
{
  char * val = getenv(key.c_str());
  return val ? val : _default;
}

std::string get_random_string(int length)
{
  // inspired from https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
  std::vector<char> char_vec{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', //
                             'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', //
                             'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', //
                             'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', //
                             'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', //
                             'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', //
                             'y', 'z'};

  std::default_random_engine rng(std::random_device{}());
  std::uniform_int_distribution<> dist(0, char_vec.size() - 1);

  std::string str(length, 0);
  std::generate_n(str.begin(), length, [&]() { return char_vec[dist(rng)]; });

  assert(std::find(str.begin(), str.end(), '\0') == str.end()); // No NULLs in final string
  return str;
}

std::string current_sec()
{
  time_t now = time(0);
  struct tm time_structure;
  char buf[32];
  time_structure = *localtime(&now);
  strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &time_structure);
  std::string sec(buf);
  return sec;
}

} // namespace

namespace weaver
{
void create_dir(std::string const & dir, unsigned mode)
{
  mkdir(dir.c_str(), mode);
}

std::string create_temp_dir()
{
  std::ostringstream ss;
  ss << get_env_var("TMPDIR", "/tmp") << "/weaver_" << current_sec() << "." << get_random_string(6);

  std::string tmp = ss.str();
  create_dir(tmp, 0700);

  std::error_code ec;
  filesystem::space_info const si = filesystem::space(tmp, ec);

  print_info("Temporary directory created: ", tmp);
  print_info("Available space on the directory: ", (si.available / 1024 / 1024 / 1024), " GB");

  return tmp;
}

void remove_file(filesystem::path const & path)
{
#ifndef NDEBUG
  bool const ret = filesystem::remove(path);

  if (ret)
    print_debug(_HERE_, " removed a file: ", path);
  else
    print_warning(_HERE_, " Could not remove file: ", path);
#else
  filesystem::remove(path);
#endif // NDEBUG
}

} // namespace weaver
