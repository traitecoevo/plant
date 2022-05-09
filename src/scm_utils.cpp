#include <plant/scm_utils.h>
#include <plant.h>

namespace plant {

std::vector<double> node_schedule_times_default(double max_time) {
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

}

//' Generate a suitable set of default node introduction times,
//' biased so that introductions are more closely packed at the
//' beginning of time, become increasingly spread out.
//'
//' The reason for the stepped distribution is to keep step sizes as
//' series of doublings.  Doing this limits the range of possible
//' introduction times from an infinite set of possible values to a
//' very limited subset of values (based on combinations of 1, 0.5,
//' 0.25, 0.125 etc).  The reason for doing this is to minimise the
//' number of unique introduction times across all species. The ODE
//' stepper needs to stop at each point where a node is introduced.
//' If each species was selecting a bunch of points that was
//' essentially unique (compared to those selected for all other
//' species), the number of unique node introductions times could
//' get very large, requiring more ODE steps.
//'
//' @title Generate Default Node Introduction Times
//' @param max_time Time to generate introduction times up to (the
//' last introduction time will be at least \code{max_time}).
//' @return Vector of introduction times.
//' @export
//' @author Rich FitzJohn, adapted from original C++ code by Daniel
//' S. Falster.
// [[Rcpp::export]]
std::vector<double> node_schedule_times_default(double max_time) {
  return plant::node_schedule_times_default(max_time);
}
