#pragma once

#include <cstdint>
#include <limits>
#include <string>

namespace weaver
{
//! Class containing a SAM record/line, the file index which the record is from, and it's order for sorting purposes.
class SAMOrderLine
{
public:
  //! The order of the SAM record
  uint64_t order{std::numeric_limits<uint64_t>::max()};

  //! The raw SAM line
  std::string sam_line{};

  //! File index the line is from
  int file_index{-1};

  //! Default constructor
  SAMOrderLine() = default;
  SAMOrderLine(SAMOrderLine const &) = default;
  SAMOrderLine(SAMOrderLine &&) = default;

  //! Assignment copy constructor is deleted.
  SAMOrderLine & operator=(SAMOrderLine const &) = default;

  //! Assignment move constructor is deleted.
  SAMOrderLine & operator=(SAMOrderLine &&) = default;

  //! Destructor defaulted.
  ~SAMOrderLine() = default;

  //! Constructs an instance with the selected file_index
  explicit SAMOrderLine(int _file_index) :
    order(std::numeric_limits<uint64_t>::max()), sam_line(), file_index(_file_index)
  {
  }

  //! Clears all fields except file_index
  inline void clear()
  {
    order = std::numeric_limits<uint64_t>::max();
    sam_line.clear();
  }
};

} // namespace weaver
