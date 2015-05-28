// -*-c++-*-
#ifndef PLANT_PLANT_INTERPOLATOR_H_
#define PLANT_PLANT_INTERPOLATOR_H_

#include <vector>
#include <RcppCommon.h> // SEXP
#include <tk/spline.h>

namespace plant {
namespace interpolator {

class Interpolator {
public:
  void init(const std::vector<double>& x_,
            const std::vector<double>& y_);
  void initialise();

  void add_point(double xi, double yi);
  void clear();

  double eval(double u) const;
  size_t size() const;

  double min() const;
  double max() const;

  std::vector<double> get_x() const;
  std::vector<double> get_y() const;

  // * R interface
  SEXP r_get_xy() const;
  std::vector<double> r_eval(std::vector<double> u) const;

private:
  void check_active() const;
  std::vector<double> x, y;
  tk::spline tk_spline;
  bool active;
};

}
}

#endif
