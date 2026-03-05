// -*-c++-*-
#ifndef PLANT_ODE_SOLVER_R6_H_
#define PLANT_ODE_SOLVER_R6_H_

#include <vector>
#include <odelia/ode_solver.hpp>

namespace plant {
namespace ode {

template <typename System>
class SolverR6 : public odelia::ode::Solver<System> {
public:
  using base_type = odelia::ode::Solver<System>;

  SolverR6(System sys_, odelia::ode::OdeControl control)
    : base_type(sys_, control) {}

  void advance_adaptive(double time) {
    base_type::advance_adaptive({base_type::time(), time});
  }

  void step_to(double time) {
    base_type::advance_fixed({base_type::time(), time});
  }

  void set_state_from_system() {
    base_type::set_state(base_type::state(), base_type::time());
  }

  System object() const {
    return base_type::get_system();
  }
};

}
}

#endif
