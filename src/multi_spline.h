// -*-c++-*-
#ifndef TREE_MULTI_SPLINE_H_
#define TREE_MULTI_SPLINE_H_

#include "spline.h"

namespace spline {

class MultiSpline {
public:
  MultiSpline(int n);

  void init(std::vector<double> x, 
	    std::vector< std::vector<double> > y);
  void init_self();

  void add_point(double xi, std::vector<double> yi);
  void reset();

  double eval(double u, size_t i) const;
  std::vector<double> eval(double u) const;

  size_t size() const;
  size_t dim() const;

  void r_init(std::vector<double> x, Rcpp::NumericMatrix y);
  void r_add_point(double xi, std::vector<double> yi);
  Rcpp::NumericVector r_get_x() const;
  Rcpp::NumericMatrix r_get_y() const;
  Rcpp::NumericMatrix r_eval(std::vector<double> u) const;
  double r_eval_1(double u, size_t i) const;
  std::vector<double> r_eval_r(double u) const;

private:
  std::vector<Spline> splines;
};

}

RCPP_EXPOSED_CLASS(spline::MultiSpline)

#endif
