// -*-c++-*-
#ifndef PLANT_PLANT_INTERPOLATOR_H_
#define PLANT_PLANT_INTERPOLATOR_H_

#include <vector>
#include <RcppCommon.h> // SEXP
#include <R_ext/Arith.h> // R_PosInf, R_NegInf
#include <tk/spline.h>

namespace plant {
namespace interpolator {

class Interpolator {
public:
  void init(const std::vector<double>& x_,
            const std::vector<double>& y_);
  void initialise();

  void add_point(double xi, double yi);
  void add_point_sorted(double xi, double yi);
  void clear();

  double eval(double u) const;
  // faster, unchecked evaluation (no check_active / bound checks): safe only
  // when the caller has already guaranteed the point lies in the domain.
  double operator()(double u) const { return tk_spline(u); }
  size_t size() const { return x.size(); }

  // These are chosen so that if a Interpolator is empty, functions looking to
  // see if they will fall outside of the covered range will always find they
  // do.  This is the same principle as R's range(numeric(0)) -> c(Inf, -Inf).
  double min() const { return size() > 0 ? x.front() : R_PosInf; }
  double max() const { return size() > 0 ? x.back() : R_NegInf; }
  void set_extrapolate(bool e);

  std::vector<double> get_x() const;
  std::vector<double> get_y() const;
  

  // * R interface
  SEXP r_get_xy() const;
	// change to const& vec?
  std::vector<double> r_eval(std::vector<double> u) const;

private:
  void check_active() const;
  std::vector<double> x, y;
  tk::spline tk_spline;
  bool active = false;
  bool extrapolate = true;
};

}
}

#endif
