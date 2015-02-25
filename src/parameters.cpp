#include <tree2/parameters.h>
#include <tree2/tree_utils.h>
#include <tree2/util.h>
#include <Rcpp.h>

namespace tree2 {

Parameters::Parameters()
  : c_ext(0.5),
    patch_area(1.0),
    Pi_0(0.25),
    n_patches(1),
    disturbance_mean_interval(30),
    cohort_schedule_max_time(NA_REAL),
    hyperpar(R_NilValue) {
}

size_t Parameters::size() const {
  return strategies.size();
}

size_t Parameters::n_residents() const {
  return static_cast<size_t>
    (std::count(is_resident.begin(), is_resident.end(), true));
}

size_t Parameters::n_mutants() const {
  return size() - n_residents();
}

// NOTE: this will be called *every time* that the object is passed in
// from R -> C++.  That's unlikely to be that often, but it does incur
// a penalty.  So don't put anything too stupidly heavy in here.
void Parameters::validate() {
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
  if (hyperpar != R_NilValue && !Rcpp::is<Rcpp::Function>(hyperpar)) {
    util::stop("hyper must be NULL or a function");
  }

  // Overwrite all strategy control objects so that they take the
  // Parameters' control object.
  for (auto& s : strategies) {
    s.control = control;
  }
}

// Separating this out just because it's a bit crap:
void Parameters::setup_cohort_schedule() {
  if (!util::is_finite(cohort_schedule_max_time)) {
    cohort_schedule_max_time = cohort_schedule_max_time_default(*this);
  }
  if (cohort_schedule_times_default.empty()) {
    cohort_schedule_times_default =
      tree2::cohort_schedule_times_default(cohort_schedule_max_time);
  }
  if (cohort_schedule_times.empty() && size() > 0) {
    for (size_t i = 0; i < size(); ++i) {
      cohort_schedule_times.push_back(cohort_schedule_times_default);
    }
  }
}

}
