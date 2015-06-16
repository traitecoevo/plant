// -*-c++-*-
#ifndef PLANT_PLANT_PARAMETERS_H_
#define PLANT_PLANT_PARAMETERS_H_

#include <vector>
#include <RcppCommon.h> // SEXP

#include <plant/control.h>
#include <plant/ffw16_strategy.h>
#include <plant/cohort_schedule.h>
#include <plant/ebt_utils.h> // Unfortunately needed for setup_cohort_schedule

// TODO: I will possibly move out the "Patch" parameters out into
// their own simple list class at some point, to make this a bit more
// coherent.
//
// TODO: Will require some free functions on the R side:
//   * add_strategy (with flag for mutant/non mutant)

namespace plant {

template <typename T>
struct Parameters {
  typedef T strategy_type;
  Parameters();
  // Data -- public for now (see github issue #17).
  double c_ext;      // Light extinction coefficient
  double patch_area; // Size of the patch (m^2)
  size_t n_patches;  // Number of patches in the metacommunity
  double disturbance_mean_interval; // Disturbance interval (years)

  std::vector<strategy_type> strategies;
  std::vector<double> seed_rain;
  std::vector<bool>   is_resident;

  // Disturbance regime.

  // Algorithm control.
  Control control;

  // Default strategy.
  strategy_type strategy_default;

  // Cohort information.
  double cohort_schedule_max_time;
  std::vector<double> cohort_schedule_times_default;
  std::vector<std::vector<double> > cohort_schedule_times;
  std::vector<double> cohort_schedule_ode_times;

  // An R function that will be used to hyperparametrise the model.
  SEXP hyperpar;

  // Some little query functions for use on the C side:
  size_t size() const;
  size_t n_residents() const;
  size_t n_mutants() const;
  void validate();
private:
  void setup_cohort_schedule();
};

template <typename T>
Parameters<T>::Parameters()
  : c_ext(0.5),
    patch_area(1.0),
    n_patches(1),
    disturbance_mean_interval(30),
    cohort_schedule_max_time(NA_REAL),
    hyperpar(util::get_from_package("FFW16_hyperpar")) {
}

template <typename T>
size_t Parameters<T>::size() const {
  return strategies.size();
}

template <typename T>
size_t Parameters<T>::n_residents() const {
  return static_cast<size_t>
    (std::count(is_resident.begin(), is_resident.end(), true));
}

template <typename T>
size_t Parameters<T>::n_mutants() const {
  return size() - n_residents();
}

// NOTE: this will be called *every time* that the object is passed in
// from R -> C++.  That's unlikely to be that often, but it does incur
// a penalty.  So don't put anything too stupidly heavy in here.
template <typename T>
void Parameters<T>::validate() {
  const size_t n_spp = size();

  // Set some defaults and check lengths.  Number of strategies is
  // taken as the "true" size.
  if (seed_rain.empty()) {
    seed_rain = std::vector<double>(n_spp, 1.0);
  } else if (seed_rain.size() != n_spp) {
    util::stop("Incorrect length seed_rain");
  }
  if (is_resident.empty()) {
    is_resident = std::vector<bool>(n_spp, true);
  } else if (is_resident.size() != n_spp) {
    util::stop("Incorrect length is_resident");
  }

  setup_cohort_schedule();
  if (cohort_schedule_times.size() != n_spp) {
    util::stop("Incorrect length cohort_schedule_times");
  }

  // This is not a lot of checking, but should be enough.  There's no
  // way of telling if the function is a good idea without running it,
  // anyway.
  if (hyperpar != R_NilValue && !util::is_function(hyperpar)) {
    util::stop("hyperpar must be NULL or a function");
  }

  // Overwrite all strategy control objects so that they take the
  // Parameters' control object.
  for (auto& s : strategies) {
    s.control = control;
  }
}

// Separating this out just because it's a bit crap:
// TODO: Consider adding this to ebt_utils.h perhaps?
template <typename T>
void Parameters<T>::setup_cohort_schedule() {
  const double max_time = cohort_schedule_max_time_default(*this);
  const bool update =
    !(util::is_finite(cohort_schedule_max_time) &&
      util::identical(cohort_schedule_max_time, max_time));

  if (update || !util::is_finite(cohort_schedule_max_time)) {
    cohort_schedule_max_time = max_time;
  }
  if (update || cohort_schedule_times_default.empty()) {
    cohort_schedule_times_default =
      plant::cohort_schedule_times_default(cohort_schedule_max_time);
  }

  if (update || (cohort_schedule_times.empty() && size() > 0)) {
    cohort_schedule_times.clear();
    for (size_t i = 0; i < size(); ++i) {
      cohort_schedule_times.push_back(cohort_schedule_times_default);
    }
  }
}

}

#endif
