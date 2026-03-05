// -*-c++-*-
#ifndef PLANT_ODE_SOLVER_COMPAT_H_
#define PLANT_ODE_SOLVER_COMPAT_H_

#include <vector>
#include <odelia/ode_solver.hpp>

namespace odelia {
namespace ode {

template <typename System>
inline void r_advance_adaptive(Solver<System>& solver, double time) {
  solver.advance_adaptive({solver.time(), time});
}

template <typename System>
inline void r_step_to(Solver<System>& solver, double time) {
  solver.advance_fixed({solver.time(), time});
}

template <typename System>
inline void r_set_state_from_system(Solver<System>& solver) {
  solver.set_state(solver.state(), solver.time());
}

}
}

#endif
