// -*-c++-*-
#ifndef TREE_SPLINE_H_
#define TREE_SPLINE_H_

#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <Rcpp.h>

namespace spline {

class Spline {
public:
  Spline();
  ~Spline();
  Spline(const Spline& x);

  void init(std::vector<double> x_, std::vector<double> y_);
  void init_self();

  void add_point(double xi, double yi);
  void reset();

  double eval(double u) const;
  size_t size() const;

  Rcpp::NumericVector r_get_x() const;
  Rcpp::NumericVector r_get_y() const;
  Rcpp::NumericMatrix r_get_xy() const;
  std::vector<double> r_eval(std::vector<double> u) const;

protected:
  std::vector<double> x, y;

private:
  gsl_interp_accel *acc;
  gsl_spline *spline;  

  void gsl_free_spline();
  void gsl_free_acc();
};

}

RCPP_EXPOSED_CLASS(spline::Spline)

#endif
