#include "integration.h"
#include "util.h"

namespace integration {

QAG::QAG(size_t rule, size_t max_iterations, double atol, double rtol)
  : q(rule),
    limit(max_iterations),
    epsabs(atol),
    epsrel(rtol),
    area(NA_REAL),
    error(NA_REAL),
    iteration(0),
    roundoff_type1(0),
    roundoff_type2(0) {
}

double QAG::integrate(util::DFunctor *f, double a, double b) {
  bool success = initialise(f, a, b);
  while (!success && iteration < limit) {
    success = refine(f);
    iteration++;
  }

  if (!success) {
    if (iteration == limit)
      ::Rf_error("Maximum number of subdivisions reached");
    else
      ::Rf_error("Could not integrate function");
  }

  area = w.total_area();

  return area;
}

double QAG::integrate_with_intervals(util::DFunctor *f,
				     intervals_type intervals) {
  std::vector<double>::const_iterator
    a     = intervals[0].begin(),
    b     = intervals[1].begin(),
    a_end = intervals[0].end();
  w.clear();
  while (a != a_end) {
    w.push_back(do_integrate(f, *a, *b));
    ++a;
    ++b;
  }
  area  = w.total_area();
  error = w.total_error();
  return area;
}

QAG::intervals_type QAG::get_last_intervals() const {
  return w.get_intervals();
}

Rcpp::List QAG::r_get_last_intervals() const {
  intervals_type tmp = get_last_intervals();
  return Rcpp::List::create(Rcpp::_["a"] = tmp[0],
			    Rcpp::_["b"] = tmp[1]);
}

// * R interface
double QAG::r_integrate(util::RFunctionWrapper fun, double a, double b) {
  return integrate(&fun, a, b);
}

double QAG::r_integrate_with_intervals(util::RFunctionWrapper fun,
				       Rcpp::List intervals) {
  util::check_length(static_cast<size_t>(intervals.size()), 2);
  intervals_type tmp;
  tmp.push_back(Rcpp::as< std::vector<double> >(intervals[0]));
  tmp.push_back(Rcpp::as< std::vector<double> >(intervals[1]));

  if (tmp[0].size() != tmp[1].size())
    ::Rf_error("Intervals must have the same length");

  return integrate_with_intervals(&fun, tmp);
}

internal::workspace::point QAG::do_integrate(util::DFunctor *f,
					     double a, double b) {
  const double tmp = q.integrate(f, a, b);
  internal::workspace::point p(a, b, tmp, q.get_last_error());
  return p;
}

bool QAG::initialise(util::DFunctor *f, double a, double b) {
  internal::workspace::point p = do_integrate(f, a, b);
  area  = p.area;
  error = p.error;
  // More information from the last integration:
  const double resabs = q.get_last_area_abs();
  const double resasc = q.get_last_area_asc();
  // Limits for checking against:
  const double tolerance = std::max(epsabs, epsrel * std::abs(area));
  const double round_off =
    50 * std::numeric_limits<double>::epsilon() * resabs;

  w.clear();
  w.push_back(p);

  roundoff_type1 = 0;
  roundoff_type2 = 0;
  iteration = 1;
  bool success = false;

  // Check that we can procede:
  if (error < round_off && error > tolerance) {
    ::Rf_error("cannot reach tolerance because of roundoff error");
  } else if ((error <= tolerance && !util::identical(error, resasc)) ||
	     util::identical(error, 0.0)) {
    success = true;
  } else if (limit == 1) {
    ::Rf_error("Need more than one iteration to achive requested accuracy");
  }

  return success;
}

bool QAG::refine(util::DFunctor *f) {
  const internal::workspace::point p0 = w.worst_point();
  const double a1 = p0.a, b2 = p0.b;
  const double a2 = (a1 + b2) * 0.5;
  const double b1 = a2;
  // Integrate lhs:
  const internal::workspace::point p1 = do_integrate(f, a1, b1);
  const double area1 = p1.area, error1 = p1.error;
  const double resasc1 = q.get_last_area_asc();
  // Integrate rhs:
  const internal::workspace::point p2 = do_integrate(f, a2, b2);
  const double area2 = p2.area, error2 = p2.error;
  const double resasc2 = q.get_last_area_asc();

  const double area12 = area1 + area2;
  const double error12 = error1 + error2;

  area  += area12  - p0.area;
  error += error12 - p0.error;
  const double tolerance = std::max(epsabs, epsrel * std::abs(area));

  // Update roundoff calculations
  if (!util::identical(resasc1, error1) &&
      !util::identical(resasc2, error2)) {
    const double delta = p0.area - area12;
    if (std::abs(delta) <= 1.0e-5 * std::abs(area12) &&
	error12 >= 0.99 * p0.error)
      roundoff_type1++;
    if (iteration >= 10 && error12 > p0.error)
      roundoff_type2++;
  }

  // Check if we have a problem
  if (error > tolerance) {
    if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
      ::Rf_error("roundoff error prevents tolerance from being achieved");
    if (subinterval_too_small(a1, a2, b2))
      ::Rf_error("bad integrand behavior found in the integration interval");
  }

  w.update(p1, p2);

  return error <= tolerance;
}

bool QAG::subinterval_too_small(double a1, double mid,
				double b2) const {
  const double e = std::numeric_limits<double>::epsilon();
  const double u = std::numeric_limits<double>::min();
  const double tmp = (1 + 100 * e) * (std::abs(mid) + 1000 * u);
  return std::abs(a1) <= tmp && std::abs(b2) <= tmp;
}

}
