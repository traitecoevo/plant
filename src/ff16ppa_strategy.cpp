// Built from  src/ff16_strategy.cpp on Thu Oct 29 11:14:41 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16ppa_strategy.h>

namespace plant {

FF16ppa_Strategy::FF16ppa_Strategy() {
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16ppa";
}

FF16ppa_Strategy::ptr make_strategy_ptr(FF16ppa_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16ppa_Strategy>(s);
}
}
