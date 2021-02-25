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
  Internals(size_t s_size=0, size_t a_size=0)
      :
      state_size(s_size),
      aux_size(a_size),
      states(s_size, 0.0),
      rates(s_size, NA_REAL) ,
      auxs(a_size, 0.0)
    {}
  size_t state_size;
  size_t aux_size;

  std::vector<double> states;
  std::vector<double> rates;
  std::vector<double> auxs;

  double state(int i) const { return states[i]; }
  double rate(int i) const { return rates[i]; }
  double aux(int i) const { return auxs[i]; }

  void set_state(int i, double v) { states[i] = v; }
  void set_rate(int i, double v) { rates[i] = v; }
  void set_aux(int i, double v) { auxs[i] = v; }

  void resize(size_t new_size, size_t new_aux_size) {
    state_size = new_size;
    aux_size = new_aux_size;
    states.resize(new_size, 0.0);
    rates.resize(new_size, NA_REAL);
    auxs.resize(new_aux_size, 0.0);
  }
};

} // namespace plant

#endif
