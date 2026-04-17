#pragma once

#include <string_view>

namespace weaver
{
//! Views a sequence in the graph
class SeqView
{
public:
  std::string_view ref{};
  std::string_view call{};
  int rid{};
  int rid_pos{};

  //! Constructor for different reference and call.
  SeqView(std::string_view _ref, std::string_view _call, int _rid, int _rid_pos) :
    ref(_ref), call(_call), rid(_rid), rid_pos(_rid_pos)
  {
  }

  //! Constructor for when ref and call are the same.
  SeqView(std::string_view _seq, int _rid, int _rid_pos) : ref(_seq), call(_seq), rid(_rid), rid_pos(_rid_pos)
  {
  }

  inline bool is_ref_equal_call()
  {
    return ref == call;
  }
};

} // namespace weaver
