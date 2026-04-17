#include "index.hpp"

#include <string>

#include "gfa.hpp"
#include "icu.hpp"
#include "logging.hpp"
#include "make_mmi.hpp"
#include "mmi.hpp"
#include "options.hpp"

namespace weaver
{
void make_index(GFA const & gfa, MMI & mmi, T_icu & icu, std::string const & vcf_fn, int k, int w)
{
  Options const & copts = (*Options::const_instance());

  print_info("Making minimizer index (mmi)...");
  mmi = make_mmi_index(gfa, vcf_fn, k, w);
  print_info("Done making MMI. Number of keys in the MMI=", mmi.num_keys());

  print_info("Making ICU index...");
  icu = make_icu_index(gfa, copts.max_icu_distance);
  print_info("ICU index ready. Number of keys in the ICU=", icu.size());
}

} // namespace weaver
