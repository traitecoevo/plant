// Inherit from FF16, like FF16r
#include <plant/models/water_strategy.h>

namespace plant {
Water_Strategy::Water_Strategy() {
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "Water";
}

Water_Strategy::ptr make_strategy_ptr(Water_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<Water_Strategy>(s);
}
}
