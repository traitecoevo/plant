// -*-c++-*-
#ifndef PLANT_PLANT_PARAMETERS_H_
#define PLANT_PLANT_PARAMETERS_H_

#include <RcppCommon.h> // SEXP
#include <vector>

#include <plant/cohort_schedule.h>
#include <plant/scm_utils.h> // Unfortunately needed for setup_cohort_schedule

// TODO: Will require some free functions on the R side:
//   * add_strategy (with flag for mutant/non mutant)

namespace plant {

template <typename T> struct SpeciesParameters {
  typedef T strategy_type;

  SpeciesParameters() { validate(); }

  // Species info
  std::vector<strategy_type> species;
  std::vector<double> birth_rate;
  std::vector<bool> is_resident;

  // Default strategy.
  strategy_type strategy_default;

  // Cohort information.
  std::vector<double> cohort_schedule_times_default;
  std::vector<std::vector<double>> cohort_schedule_times;
  std::vector<double> cohort_schedule_ode_times;

  // Some little query functions for use on the C side:
  size_t size() const;
  size_t n_residents() const;
  size_t n_mutants() const;
  void validate();

private:
  void setup_cohort_schedule();
};

template <typename T> size_t SpeciesParameters<T>::size() const {
  return species.size();
}

template <typename T> size_t SpeciesParameters<T>::n_residents() const {
  return static_cast<size_t>(
      std::count(is_resident.begin(), is_resident.end(), true));
}

template <typename T> size_t SpeciesParameters<T>::n_mutants() const {
  return size() - n_residents();
}

// NOTE: this will be called *every time* that the object is passed in
// from R -> C++.  That's unlikely to be that often, but it does incur
// a penalty.  So don't put anything too stupidly heavy in here.
template <typename T> void SpeciesParameters<T>::validate() {
  const size_t n_spp = size();

  // Set some defaults and check lengths
  if (birth_rate.empty()) {
    birth_rate = std::vector<double>(n_spp, 1.0);
  } else if (birth_rate.size() != n_spp) {
    util::stop("Incorrect length birth_rate");
  }
  if (is_resident.empty()) {
    is_resident = std::vector<bool>(n_spp, true);
  } else if (is_resident.size() != n_spp) {
    util::stop("Incorrect length is_resident");
  }

  // Set up cohort introductions
  setup_cohort_schedule();
  if (cohort_schedule_times.size() != n_spp) {
    util::stop("Incorrect length cohort_schedule_times");
  }
}

// Separating this out just because it's a bit crap:
// TODO: Consider adding this to scm_utils.h perhaps?
template <typename T> void SpeciesParameters<T>::setup_cohort_schedule() {
  if ((cohort_schedule_times.empty() && size() > 0)) {
    cohort_schedule_times.clear();
    for (size_t i = 0; i < size(); ++i) {
      cohort_schedule_times.push_back(cohort_schedule_times_default);
    }
  }
}

} // namespace plant

#endif
