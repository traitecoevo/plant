// -*-c++-*-
#ifndef TREE_QAG_H_
#define TREE_QAG_H_

#include <tree2/qk.h>
#include <tree2/qag_internals.h>
#include <tree2/util.h> // util::stop
#include <RcppCommon.h> // SEXP

namespace tree2 {
namespace quadrature {

// This is the "QAG" algorithm from QUADPACK -- quadrature, adaptive,
// Gaussian.  Does not handle infinite intervals or singularities.
class QAG {
public:
  QAG(); // need a default constructor.
  QAG(size_t rule, size_t max_iterations, double atol, double rtol);

  template <typename Function>
  double integrate(Function f, double a, double b);
  template <typename Function>
  double integrate_with_intervals(Function f, intervals_type intervals);
  template <typename Function>
  double integrate_with_last_intervals(Function f, double a, double b);

  double get_last_area()       const {return area;}
  double get_last_error()      const {return error;}
  size_t get_last_iterations() const {return iteration;}
  intervals_type get_last_intervals() const;

  bool is_adaptive() const {return adaptive;}

  // * R interface.  These avoid referencing full Rcpp types until
  // implementation in qag.cpp.
  double r_integrate(SEXP f, double a, double b);
  double r_integrate_with_intervals(SEXP f, SEXP intervals);
  double r_integrate_with_last_intervals(SEXP f, double a, double b);

private:
  template <typename Function>
  double integrate_adaptive(Function f, double a, double b);
  template <typename Function>
  double integrate_fixed(Function f, double a, double b);
  template <typename Function>
  internal::workspace::point do_integrate(Function f, double a, double b);
  template <typename Function>
  bool initialise(Function f, double a, double b);
  template <typename Function>
  bool refine(Function f);
  static bool subinterval_too_small(double a1, double mid, double b2);

  bool adaptive;

  QK q;
  internal::workspace w;

  // Control parameters
  size_t limit;
  double epsabs;
  double epsrel;

  // Intermediates
  double area, error;
  size_t iteration;
  size_t roundoff_type1, roundoff_type2;
};

template <typename Function>
double QAG::integrate(Function f, double a, double b) {
  return adaptive ? integrate_adaptive(f, a, b) : integrate_fixed(f, a, b);
}

template <typename Function>
double QAG::integrate_adaptive(Function f, double a, double b) {
  bool success = initialise(f, a, b);
  while (!success && iteration < limit) {
    success = refine(f);
    iteration++;
  }

  if (!success) {
    if (iteration == limit) {
      util::stop("Maximum number of subdivisions reached");
    } else {
      util::stop("Could not integrate function");
    }
  }

  area = w.total_area();

  return area;
}

template <typename Function>
double QAG::integrate_fixed(Function f, double a, double b) {
  area  = q.integrate(f, a, b);
  error = q.get_last_error();
  return area;
}

template <typename Function>
double QAG::integrate_with_intervals(Function f,
                                     intervals_type intervals) {
  if (!adaptive) {
    util::stop("This really does not make any sense...");
  }

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

template <typename Function>
double QAG::integrate_with_last_intervals(Function f,
                                          double a, double b) {
  intervals_type intervals = w.get_intervals();
  if (intervals[0].empty()) {
    util::stop("No stored intervals to use");
  }
  intervals_type intervals_scaled =
    internal::rescale_intervals(intervals, a, b);
  return integrate_with_intervals(f, intervals_scaled);
}

template <typename Function>
internal::workspace::point QAG::do_integrate(Function f, double a, double b) {
  const double tmp = q.integrate(f, a, b);
  internal::workspace::point p(a, b, tmp, q.get_last_error());
  return p;
}

template <typename Function>
bool QAG::initialise(Function f, double a, double b) {
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
    util::stop("cannot reach tolerance because of roundoff error");
  } else if ((error <= tolerance && !util::identical(error, resasc)) ||
             util::identical(error, 0.0)) {
    success = true;
  } else if (limit == 1) {
    util::stop("Need more than one iteration to achive requested accuracy");
  }

  return success;
}

template <typename Function>
bool QAG::refine(Function f) {
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
  if (!util::identical(resasc1, error1) && !util::identical(resasc2, error2)) {
    const double delta = p0.area - area12;
    if (std::abs(delta) <= 1.0e-5 * std::abs(area12) &&
	error12 >= 0.99 * p0.error) {
      roundoff_type1++;
    }
    if (iteration >= 10 && error12 > p0.error) {
      roundoff_type2++;
    }
  }

  // Check if we have a problem
  if (error > tolerance) {
    if (roundoff_type1 >= 6 || roundoff_type2 >= 20) {
      util::stop("roundoff error prevents tolerance from being achieved");
    }
    if (subinterval_too_small(a1, a2, b2)) {
      util::stop("bad integrand behavior found in the integration interval");
    }
  }

  w.update(p1, p2);

  return error <= tolerance;
}

}
}

#endif
