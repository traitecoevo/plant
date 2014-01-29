#include "fake_light_environment.h"

namespace interpolator {

FakeLightEnvironment::FakeLightEnvironment(std::vector<double> times_,
					   std::vector<Interpolator> splines_)
  : times(times_), splines(splines_), time(0) {
  initialise();
}

FakeLightEnvironment::FakeLightEnvironment(std::vector<double> times_,
					   Rcpp::List splines_)
  : times(times_), time(0) {
  for (int i = 0; i < splines_.size(); ++i)
    splines.push_back(Rcpp::as<Interpolator>(splines_[i]));
  initialise();
}

void FakeLightEnvironment::initialise() {
  util::check_length(times.size(), splines.size());
  if (times.size() < 2)
    Rcpp::stop("Need at least two times");
  if (!util::identical(times[0], 0.0))
    Rcpp::stop("First time must be zero");
  if (!util::is_sorted(times.begin(), times.end()))
    Rcpp::stop("Times must be sorted");
  time_max = times.back();
  set_time(0);
}

void FakeLightEnvironment::set_time(double t) {
  if (t < 0)
    Rcpp::stop("Can't set negative time");
  if (t > time_max)
    Rcpp::stop("Time exceeds maximum");
  time = t;

  std::vector<double>::const_iterator it1 =
    std::lower_bound(times.begin(), times.end(), t);
  std::vector<double>::const_iterator
    it0 = it1 == times.begin() ? it1 : it1 - 1;

  current = merge(time,
		  *it0, splines[static_cast<size_t>(it0 - times.begin())],
		  *it1, splines[static_cast<size_t>(it1 - times.begin())]);
}

Interpolator
FakeLightEnvironment::merge(double t,
			    double t0, const Interpolator& s0,
			    double t1, const Interpolator& s1) const {
  if (util::identical(t, t1))
    return s1;

  std::vector<double> h0 = s0.get_x(), h1 = s1.get_x();
  std::vector<double> h(h0.size() + h1.size());
  std::vector<double>::iterator it =
    std::set_union(h0.begin(), h0.end(), h1.begin(), h1.end(), h.begin());
  h.resize(static_cast<size_t>(it - h.begin()));
  std::vector<double> y;
  y.reserve(h.size());

  const double w = (t - t0) / (t1 - t0);
  for (std::vector<double>::const_iterator
	 hi = h.begin(); hi != h.end(); ++hi) {
    const double
      y0 = *hi > s0.max() ? 1.0 : s0.eval(*hi),
      y1 = *hi > s1.max() ? 1.0 : s1.eval(*hi);
    y.push_back(y0 + (y1 - y0) * w);
  }

  Interpolator ret;
  ret.init(h, y);
  return ret;
}

}
