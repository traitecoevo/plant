// -*-c++-*-
#ifndef TREE_EBT_MUTANT_RUNNER_H_
#define TREE_EBT_MUTANT_RUNNER_H_

#include "ebt.h"
#include "fake_light_environment.h"

namespace model {

// Implementing the bare minimum for now.  The other option would be
// to inherit or compose EBT and then get in and muck with the time
// setting.

class EBTMutantRunner : public EBT {
public:
  EBTMutantRunner(Parameters  p, interpolator::FakeLightEnvironment e);
  EBTMutantRunner(Parameters *p, interpolator::FakeLightEnvironment e);

  ode::iterator_const set_ode_values(double time_, ode::iterator_const it);

private:
  interpolator::FakeLightEnvironment light_environment;
};

}

#endif
