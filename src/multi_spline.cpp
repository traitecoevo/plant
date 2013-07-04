#include "multi_spline.h"
#include "util.h"

namespace spline {

MultiSpline::MultiSpline(size_t n) : splines(n) {
  if (n == 0)
    ::Rf_error("Must have at least one response");
}

void MultiSpline::init(std::vector<double>   x, 
	  std::vector< std::vector<double> > y) {
  for (size_t i = 0; i < dim(); ++i)
    splines[i].init(x, y[i]);
}

void MultiSpline::init_self() {
  for (std::vector<Spline>::iterator si = splines.begin();
       si != splines.end(); ++si)
    si->init_self();
}

void MultiSpline::add_point(double xi, std::vector<double> yi) {
  for (size_t i = 0; i < dim(); ++i)
    splines[i].add_point(xi, yi[i]);
}

void MultiSpline::clear() {
  for (std::vector<Spline>::iterator si = splines.begin();
       si != splines.end(); ++si)
    si->clear();
}

double MultiSpline::eval(double u, size_t i) const {
  return splines[i].eval(u);
}

std::vector<double> MultiSpline::eval(double u) const {
  std::vector<double> ret;
  ret.reserve(dim());
  for (std::vector<Spline>::const_iterator si = splines.begin();
       si != splines.end(); ++si)
    ret.push_back(si->eval(u));
  return ret;
}

size_t MultiSpline::size() const {
  return splines.begin()->size();
}
size_t MultiSpline::dim() const {
  return splines.size();
}

void MultiSpline::r_init(std::vector<double> x, Rcpp::NumericMatrix y) {
  util::check_length(static_cast<size_t>(y.ncol()), dim());
  util::check_length(static_cast<size_t>(y.nrow()), x.size());
  std::vector< std::vector<double> > yy;
  for (size_t i = 0; i < dim(); ++i) {
    Rcpp::NumericVector yi = y(Rcpp::_, static_cast<int>(i));
    std::vector<double> yiv(yi.begin(), yi.end());
    yy.push_back(yiv);
  }
  init(x, yy);
}

void MultiSpline::r_add_point(double xi, std::vector<double> yi) {
  util::check_length(yi.size(), dim());
  add_point(xi, yi);
}

std::vector<double> MultiSpline::get_x() const {
  return splines.begin()->get_x();
}

Rcpp::NumericMatrix MultiSpline::r_get_y() const {
  Rcpp::NumericMatrix ret(static_cast<int>(size()),
			  static_cast<int>(dim()));
  for (size_t i = 0; i < dim(); ++i) {
    // NOTE: This involves more copying than necessary (see
    // ode_solver.h:Solver::r_run for a better solution).
    std::vector<double> yi = splines[i].get_y();
    Rcpp::NumericVector yi_r(yi.begin(), yi.end());
    ret(Rcpp::_, static_cast<int>(i)) = yi_r;
  }
  return ret;
}

Rcpp::NumericMatrix MultiSpline::r_eval(std::vector<double> u) const {
  Rcpp::NumericMatrix ret(static_cast<int>(u.size()),
			  static_cast<int>(dim()));
  for (size_t i = 0; i < dim(); ++i)
    for (size_t j = 0; j < u.size(); ++j)
      ret(j,i) = eval(u[j], i);
  return ret;
}

double MultiSpline::r_eval_1(double u, size_t i) const {
  return eval(u, util::check_bounds_r(i, dim()));
}

std::vector<double> MultiSpline::r_eval_r(double u) const {
  return eval(u);
}

}
