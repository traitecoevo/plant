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

  SpeciesParameters() {}

  // Species info
  std::vector<strategy_type> species;
  std::vector<double> birth_rate;
  std::vector<bool> is_resident;

  // Cohort information.
  std::vector<double> cohort_schedule_times_default;
  std::vector<std::vector<double>> cohort_schedule_times;
  std::vector<double> cohort_schedule_ode_times;

  // Some little query functions for use on the C side:
  size_t size() const;
  size_t n_residents() const;
  size_t n_mutants() const;
  void validate();

  // Expose cohort schedule fns to r
  void r_setup_cohort_schedule(const Control &c) { setup_cohort_schedule(c); }

private:
  void setup_cohort_schedule(const Control &c);
  std::vector<double> generate_default_cohort_schedule_times(double max_time);
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

  if (cohort_schedule_times.size() != n_spp) {
    util::stop("Incorrect length cohort_schedule_times, did you run "
               "`setup_cohort_schedule()`?");
  }
}

// Separating this out just because it's a bit crap:
// TODO: Consider adding this to scm_utils.h perhaps?
template <typename T>
void SpeciesParameters<T>::setup_cohort_schedule(const Control &c) {

  if (cohort_schedule_times_default.empty()) {
    double max_time = c.patch_max_lifetime;
    cohort_schedule_times_default =
      generate_default_cohort_schedule_times(max_time);
  }

  if ((cohort_schedule_times.empty() && size() > 0)) {
    cohort_schedule_times.clear();
    for (size_t i = 0; i < size(); ++i) {
      cohort_schedule_times.push_back(cohort_schedule_times_default);
    }
  }
}

//' Generate a suitable set of default cohort introduction times,
//' biased so that introductions are more closely packed at the
//' beginning of time, become increasingly spread out.
//'
//' The reason for the stepped distribution is to keep step sizes as
//' series of doublings.  Doing this limits the range of possible
//' introduction times from an infinite set of possible values to a
//' very limited subset of values (based on combinations of 1, 0.5,
//' 0.25, 0.125 etc).  The reason for doing this is to minimise the
//' number of unique introduction times across all species. The ODE
//' stepper needs to stop at each point where a cohort is introduced.
//' If each species was selecting a bunch of points that was
//' essentially unique (compared to those selected for all other
//' species), the number of unique cohort introductions times could
//' get very large, requiring more ODE steps.
//'
//' @title Generate Default Cohort Introduction Times
//' @param max_time Time to generate introduction times up to (the
//' last introduction time will be at least \code{max_time}).
//' @return Vector of introduction times.
//' @export
//' @author Rich FitzJohn, adapted from original C++ code by Daniel
//' S. Falster.
template <typename T>
std::vector<double> SpeciesParameters<T>::generate_default_cohort_schedule_times(double max_time) {
  const double multiplier = 0.2, min_step_size = 1e-5, max_step_size = 2.0;
  if (min_step_size <= 0) {
    util::stop("The minimum step size must be greater than zero");
  }
  double dt = 0.0, time = 0.0;
  std::vector<double> times;
  times.push_back(time);
  while (time <= max_time) {
    dt = std::exp2(std::floor(std::log2(time * multiplier)));
    time += util::clamp(dt, min_step_size, max_step_size);
    times.push_back(time);
  }
  // Drop the last time; that's not going to be needed:
  times.resize(times.size() - 1);
  return times;
}

} // namespace plant

#endif
