// -*-c++-*-
#ifndef TREE_SPLINE_H_
#define TREE_SPLINE_H_

#include <vector>
#include <cstddef>

// Help for GSL -- see:
// https://www.gnu.org/software/gsl/manual/html_node/Inline-functions.html
#define HAVE_INLINE
// Could cause some platform non-independence.

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <Rcpp.h>

#include "util.h"

namespace spline {

class Spline {
public:
  typedef util::PtrWrapper<Spline> ptr;
  Spline();
  Spline(bool is_akima, bool is_linear);
  ~Spline();
  Spline(const Spline& x);
  Spline& operator=(Spline other);

  std::string type() const;

  void init(const std::vector<double>& x_,
	    const std::vector<double>& y_);
  void init_self();

  void add_point(double xi, double yi);
  void clear();

  double eval(double u) const;
  double deriv(double u) const;
  size_t size() const;

  double min() const;
  double max() const;

  std::vector<double> get_x() const;
  std::vector<double> get_y() const;

  // * R interface
  Rcpp::NumericMatrix r_get_xy() const;
  std::vector<double> r_eval(std::vector<double> u) const;
  std::vector<double> r_deriv(std::vector<double> u) const;

private:
  void gsl_free_interp();
  void gsl_free_interp_accel();
  void check_initialised() const;

  static const gsl_interp_type*
  select_interpolation_type(bool is_akima, bool is_linear);
  static const double *
  ref(const std::vector<double>& x) {return &x.front();}

  std::vector<double> x, y;
  gsl_interp_type const *interpolation_type;
  gsl_interp            *interp;
  gsl_interp_accel      *interp_accel;
};

}

RCPP_EXPOSED_CLASS(spline::Spline)

#endif
