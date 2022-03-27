#ifndef PLANT_EXTRINSIC_DRIVERS_H
#define PLANT_EXTRINSIC_DRIVERS_H

#include <plant/util.h>
#include <plant/interpolator.h>
#include <type_traits>

namespace plant {
namespace {

class Function {
public:
  Function() = default;

  Function(std::vector<double> const &x, std::vector<double> const &y) {
    variable.init(x, y);
    variable.set_extrapolate(false);
    is_variable = true;
  }

  Function(double k) {
    constant = k;
    is_variable = false;
  }

  double evaluate(double u) const {
    if (is_variable) {
      return variable.eval(u);
    } else {
      return constant;
    }
  }

  std::vector<double> evaluate_range(const std::vector<double> &u) const {
    if (is_variable) {
      return variable.r_eval(u);
    } else {
      return std::vector<double>(u.size(), constant);
    }
  }

  void set_extrapolate(bool extrapolate) {
    variable.set_extrapolate(extrapolate);
  }

private:
    interpolator::Interpolator variable;
    double constant;
    bool is_variable;
};
}

class ExtrinsicDrivers {
  using variable = interpolator::Interpolator;
public:
  // this will override any previously defined drivers with the same name
  void set_constant(std::string driver_name, double k) {
    drivers[driver_name] = Function(k);
  }

  // initialise spline of driver with x, y control points
  void set_variable(std::string driver_name, std::vector<double> const &x, std::vector<double> const &y) {
    drivers[driver_name] = Function(x, y);
  }

  void set_extrapolate(std::string driver_name, bool extrapolate) {
    drivers.at(driver_name).set_extrapolate(extrapolate);
  }

  // evaluate/query interpolated spline for driver at point u, return s(u), where s is interpolated function
  double evaluate(std::string driver_name, double u) const {
    return drivers.at(driver_name).evaluate(u);
  }

  // evaluate/query interpolated spline for driver at vector of points, return vector of values
  std::vector<double> evaluate_range(std::string driver_name, std::vector<double> u) const {
    return drivers.at(driver_name).evaluate_range(u);
  }

  // returns the name of each active driver - useful for R output
  std::vector <std::string> get_names() {
    auto ret = std::vector<std::string>();
    for (auto const &driver: drivers) {
      ret.push_back(driver.first);
    }
    return ret;
  }

private:
  std::unordered_map <std::string, Function> drivers;
};

}

#endif //PLANT_EXTRINSIC_DRIVERS_H
