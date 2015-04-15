#include "ebt_mutant_runner.h"

namespace model {

// TODO: check that all strategies are marked as mutants?  Probably
// does not matter here because that only affects how the light
// environment is calculated and we side-step that.
//
// TODO: Set control.environment_light_skip to FALSE, or enforce that
// it is set that way.

// Bypass set_ode_values mucking with the light environment; it does
// this:
//   1. set values
//   2. update light environment
//   3. update physiological variables
// We set our own light environment and then force the EBT to skip the
// check.
EBTMutantRunner::EBTMutantRunner(Parameters p,
				 interpolator::FakeLightEnvironment e)
  : EBT(p),
    light_environment(e) {
  if (!p.get_control().environment_light_skip)
    Rcpp::stop("Turn environment_light_skip to true");
}
EBTMutantRunner::EBTMutantRunner(Parameters *p,
				 interpolator::FakeLightEnvironment e)
  : EBT(p),
    light_environment(e) {
  if (!p->get_control().environment_light_skip)
    Rcpp::stop("Turn environment_light_skip to true");
}

ode::iterator_const EBTMutantRunner::set_ode_values(double time,
						    ode::iterator_const it) {
  patch.set_light_environment(light_environment(time));
  return EBT::set_ode_values(time, it);
}

}
