// -*-c++-*-
#ifndef PLANT_PLANT_INTERNALS_MINIMAL_H_
#define PLANT_PLANT_INTERNALS_MINIMAL_H_

#define HEIGHT_INDEX 0
#define MORTALITY_INDEX 1
#define FECUNDITY_INDEX 2

#include <memory> // std::shared_ptr
#include <plant/ode_interface.h>
#include <vector>
// #include <plant/plant_internals.h>

// TODO: extra_state bounds, upper and lower limits
namespace plant {

class Internals {
public:
  Internals(size_t s_size=0)
      :
      state_size(s_size),
      states(s_size, 0.0),
      rates(s_size, NA_REAL) {
    }
  size_t state_size;
  std::vector<double> states;
  std::vector<double> rates;


  double state(int i) const { return states[i]; }
  double  rate(int i) const { return rates[i]; }
  void set_state(int i, double v) { states[i] = v; }
  void set_rate(int i, double v) { rates[i] = v; }

  void resize(size_t new_size) {
    state_size = new_size;
    states.resize(new_size, 0.0);
    rates.resize(new_size, NA_REAL);
  }
};

} // namespace plant

#endif
