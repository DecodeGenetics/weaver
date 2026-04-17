#include <algorithm>
#include <vector>

#include <weaver/edit.hpp>

#include "test.hpp"

namespace weaver::test
{
static void test_edit_index_sort()
{
  {
    std::vector<int> vec = {-1, 3};
    REQUIRE(std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
  }

  {
    std::vector<int> vec = {1, -3};
    REQUIRE(std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
  }

  {
    std::vector<int> vec = {1, -2};
    REQUIRE(std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
  }

  {
    std::vector<int> vec = {-2, 1};
    REQUIRE(std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
  }

  {
    std::vector<int> vec = {-1, 4, -10, 20};
    REQUIRE(std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
  }

  {
    std::vector<int> vec = {-1, -10, 4, 20};
    REQUIRE(not std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
    std::sort(vec.begin(), vec.end(), is_less_edit_index);
    REQUIRE(std::is_sorted(vec.begin(), vec.end(), is_less_edit_index));
  }
}

//! \cond TESTS
TEST_CASE("Tests for edit.cpp.", "[edit]")
{
  test_edit_index_sort();
}
//! \endcond

} // namespace weaver::test
