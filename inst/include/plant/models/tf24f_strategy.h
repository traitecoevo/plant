// -*-c++-*-
#ifndef PLANT_PLANT_TF24F_STRATEGY_H_
#define PLANT_PLANT_TF24F_STRATEGY_H_

#include <plant/models/tf24_strategy.h>

namespace plant {

// TF24f ("f" for fast / forecasting): a variant of TF24 that will let the
// optimal leaf hydraulic state chase its optimum via an extra ODE state
// (gradient-ascent), instead of re-solving the nested leaf optimisation /
// root-collar root-find from scratch at every step (issue #525). It inherits
// TF24_Strategy and reuses TF24_Environment + TF24_Pars; only the leaf-solve
// hook and the extra state are overridden.
//
// Phase A: pure inheritance skeleton — overrides nothing behaviourally, so a
// TF24f run reproduces TF24 exactly. The acclimation biology is added in later
// phases.
class TF24f_Strategy : public TF24_Strategy {
public:
  typedef std::shared_ptr<TF24f_Strategy> ptr;
  TF24f_Strategy();
};

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s);

}

#endif
