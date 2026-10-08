// -*-c++-*-
#ifndef PLANT_PLANT_PARAMETERS_H_
#define PLANT_PLANT_PARAMETERS_H_

#include <vector>
#include <RcppCommon.h> // SEXP

#include <plant/control.h>
#include <plant/models/ff16_strategy.h>
#include <plant/node_schedule.h>
#include <plant/scm_utils.h> // Unfortunately needed for setup_node_schedule

#include <memory>
#include <plant/disturbance_regime.h>
#include <plant/disturbances/no_disturbance.h>
#include <plant/disturbances/weibull_disturbance.h>

// TODO(#483): I will possibly move out the "Patch" parameters out into
// their own simple list class at some point, to make this a bit more
// coherent.
//
// TODO(#483): Will require some free functions on the R side:
//   * add_strategy (with flag for mutant/non mutant)

namespace plant {

template <typename T, typename E>
struct Parameters {
  typedef T strategy_type;
  typedef E environment_type;

  Parameters() :
    patch_area(1.0),
    n_patches(1),
    patch_type("meta-population"),
    max_patch_lifetime(105.32), // designed to agree with Daniel's implementation
    initial_time(0.0)
  {
    validate();
  }

  // Data -- public for now (see github issue #17).
  double patch_area; // Size of the patch (m^2)
  size_t n_patches;  // Number of patches in the metacommunity
  std::string patch_type;
  double max_patch_lifetime; // Disturbance interval (years)
  std::vector<strategy_type> strategies;
  // Every active value the strategies here carry. They are copies the patch does
  // not run, and a recording does not read them -- but they hold slots, and a
  // slot left registered is one the next clear reissues.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, strategies, strategy_default);
  }

  // Owned, and shared with every copy of these parameters and with the patch
  // that reads it. A raw pointer here was allocated on every validate() and
  // freed nowhere, and the leak was load-bearing: the pointer is copied into
  // each Parameters copy and aliased by the patch, so freeing it once would have
  // left the others dangling.
  std::shared_ptr<Disturbance_Regime> disturbance;

  // Default strategy.
  strategy_type strategy_default;

  // Node information.
  std::vector<double> node_schedule_times_default;
  std::vector<std::vector<double> > node_schedule_times;
  std::vector<double> ode_times;
  // The size of the step that reached each of ode_times, NaN first. Recorded
  // beside the times because a schedule replayed by times alone takes different
  // steps: a size differenced back out of two recorded times is not the size
  // that was taken. A run carrying both replays itself exactly.
  std::vector<double> ode_step_sizes;

  // Initial patch state. When initial_state is non-empty the patch is seeded
  // with these nodes at reset() instead of starting empty -- used to resume an
  // exported patch run or to seed an arbitrary initial size distribution.
  // Carried here so a run stays fully self-describing/serialisable (see
  // agents.md), and so the seeding survives the reset() at the start of every
  // SCM::run()/refine_schedule().
  //   initial_state             flat ODE state (all nodes, then environment),
  //                             in the order Patch::set_ode_state expects.
  //   n_initial_cohorts         number of nodes per species.
  //   initial_node_times        per-node introduction time (flat across species).
  //   initial_patch_density     per-node patch-age density at birth (flat).
  //   initial_pr_patch_survival per-node pr_patch_survival at birth (flat).
  //   initial_time              patch age to resume at (0 for a fresh seed).
  // Length consistency is checked in Patch::set_initial_state(), where the
  // node/ode sizes are known.
  std::vector<double> initial_state;
  std::vector<size_t> n_initial_cohorts;
  std::vector<double> initial_node_times;
  std::vector<double> initial_patch_density;
  std::vector<double> initial_pr_patch_survival;
  double initial_time;

  // Some little query functions for use on the C side:
  size_t size() const;
  void validate();

private:
  void setup_node_schedule();
};

template <typename T, typename E>
size_t Parameters<T,E>::size() const {
  return strategies.size();
}

// NOTE: this will be called *every time* that the object is passed in
// from R -> C++.  That's unlikely to be that often, but it does incur
// a penalty.  So don't put anything too stupidly heavy in here.
template <typename T, typename E>
void Parameters<T,E>::validate() {
  const size_t n_spp = size();

  setup_node_schedule();
  if (node_schedule_times.size() != n_spp) {
    util::stop("Incorrect length node_schedule_times");
  }

  // Disturbances used to describe evolution of a metapopulation of patches
  // when calculating fitness, otherwise defaults to fixed-duration run without
  // disturbance
  if(patch_type == "meta-population") {
    disturbance = std::make_shared<Weibull_Disturbance_Regime>(max_patch_lifetime);
  }
  else {
    disturbance = std::make_shared<No_Disturbance>();
  }
}

// Separating this out just because it's a bit crap:
// TODO(#483): Consider adding this to scm_utils.h perhaps?
template <typename T, typename E>
void Parameters<T,E>::setup_node_schedule() {
  node_schedule_times_default =
      plant::node_schedule_times_default(max_patch_lifetime);

  if ((node_schedule_times.empty() && size() > 0)) {
    node_schedule_times.clear();
    for (size_t i = 0; i < size(); ++i) {
      node_schedule_times.push_back(node_schedule_times_default);
    }
  }
}
}

#endif
