#include <plant/models/tf24f_strategy.h>

namespace plant {

// Phase A: pure inheritance skeleton. The base TF24_Strategy constructor sets
// collect_all_auxiliary, builds the state/aux index maps (refresh_indices) and
// the name; here we only override the name so the strategy is reported as TF24f.
TF24f_Strategy::TF24f_Strategy() {
  name = "TF24f";
}

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Strategy>(s);
}

}
