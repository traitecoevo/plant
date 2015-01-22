#include <tree2/ebt_utils.h>
#include <tree2.h>

namespace tree2 {

std::vector<double> cohort_schedule_times_default(double max_time) {
  const double multiplier=0.2, min_step_size=1e-5, max_step_size=2.0;
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

double cohort_schedule_max_time_default(const Parameters& p) {
  Disturbance d(p.disturbance_mean_interval);
  return d.cdf(p.control.schedule_patch_survival);
}

CohortSchedule cohort_schedule_default(const Parameters& p) {
  const double max_time = cohort_schedule_max_time_default(p);
  CohortSchedule schedule(0);
  schedule.r_set_max_time(max_time);
  std::vector<double> times = cohort_schedule_times_default(max_time);
  if (times.size() < 1) {
    util::stop("Did not generate any times, surprisingly");
  }
  return schedule.expand(p.size(), times);
}

// NOTE: can't be a const reference because of validation, which might
// set defaults within p.
CohortSchedule make_cohort_schedule(const Parameters& p) {
  CohortSchedule ret(p.size());
  ret.r_set_max_time(p.cohort_schedule_max_time);
  ret.set_times(p.cohort_schedule_times);
  ret.r_set_ode_times(p.cohort_schedule_ode_times);
  return ret;
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
//' @param multiplier The rate of increase of step size with time.
//' The greater the number the faster step size will increase.
//' @param min_step_size The smallest gap between introduction times
//' (must be greater than zero, and will be the first introduction
//' time).
//' @param max_step_size The largest gap between introduction times
//' (may be infinite).
//' @return Vector of introduction times.
//' @export
//' @author Rich FitzJohn, adapted from original C++ code by Daniel
//' S. Falster.
// [[Rcpp::export]]
std::vector<double>
cohort_schedule_times_default(double max_time) {
  return tree2::cohort_schedule_times_default(max_time);
}

// These two have the slighty odd export to prevent argument dependent
// lookup, which renders the version in tree2:: a candidate.  It might
// be better to have RcppExports refer to
// ::cohort_schedule_max_time_default, but that's not how it's
// implemented.
// [[Rcpp::export(cohort_schedule_max_time_default)]]
double r_cohort_schedule_max_time_default(const tree2::Parameters& p) {
  return tree2::cohort_schedule_max_time_default(p);
}

// [[Rcpp::export(cohort_schedule_default)]]
tree2::CohortSchedule r_cohort_schedule_default(const tree2::Parameters& p) {
  return tree2::cohort_schedule_default(p);
}

// [[Rcpp::export(make_cohort_schedule)]]
tree2::CohortSchedule r_make_cohort_schedule(const tree2::Parameters& p) {
  return tree2::make_cohort_schedule(p);
}
