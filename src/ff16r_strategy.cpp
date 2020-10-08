// Built from  src/ff16_strategy.cpp on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16r_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
FF16r_Strategy::FF16r_Strategy() {
  // height above hmat at which allocation to reproduction is half its max value
  a_f2   = 2; // [dimensionless]

  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16r";
}


// [eqn 16] Fraction of production allocated to reproduction
double FF16r_Strategy::fraction_allocation_reproduction(double height) const {
  if(height <= hmat)
    return 0.0;
  else
    return a_f1 * (height - hmat) / (a_f2  + (height - hmat));
}


FF16r_Strategy::ptr make_strategy_ptr(FF16r_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16r_Strategy>(s);
}
}
