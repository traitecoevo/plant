// Built from  src/ff16_strategy.cpp on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16r_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
FF16r_Strategy::FF16r_Strategy() {
  // height above hmat at which allocation to reproduction is half its max value

  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16r";
}


FF16r_Strategy::ptr make_strategy_ptr(FF16r_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16r_Strategy>(s);
}
}
