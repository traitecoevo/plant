// -*-c++-*-
#ifndef PLANT_PLANT_INTERNALS_MINIMAL_H_
#define PLANT_PLANT_INTERNALS_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <plant/ode_interface.h>
#include <vector>
// #include <plant/plant_internals.h>

// TODO: extra_state bounds, upper and lower limits
namespace plant {

class Internals {
public:
  Internals(int s_size=0)
      :
      state_size(s_size),
      states(s_size, 0.0),
      rates(s_size, NA_REAL) {
    }
  int state_size;
  std::vector<double> states;
  std::vector<double> rates;

  double state(int i) const { return states[i]; }
  double  rate(int i) const { return rates[i]; }
  void set_state(int i, double v) { states[i] = v; }
  void set_rate(int i, double v) { rates[i] = v; }
};

} // namespace plant

#endif
