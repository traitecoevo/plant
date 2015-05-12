// -*-c++-*-
#ifndef TREE_TREE_ODE_RUNNER_H_
#define TREE_TREE_ODE_RUNNER_H_

#include <tree/ode_solver.h>

namespace tree {
namespace ode {

// This is a little wrapper class that is meant to simplify the
// difficuly of ownership semantics around the solver and object.
template <typename T>
class Runner {
public:
  Runner(T obj_, OdeControl control) : obj(obj_), solver(obj, control) {}
  double time() const {return solver.get_time();}
  state_type state() const {return solver.get_state();}
  void set_state(state_type y, double time) {
    util::check_length(y.size(), obj.ode_size());
    ode::internal::set_ode_state(obj, y, time);
    solver.reset(obj);
    solver.set_state_from_system(obj);
  }
  void set_state_from_system() {solver.set_state_from_system(obj);}
  std::vector<double> times() const {return solver.get_times();}
  T object() const {return obj;}
  void advance(double time) {solver.advance(obj, time);}
  void advance_fixed(std::vector<double> times) {
    solver.advance_fixed(obj, times);
  }
  void step() {solver.step(obj);}
  void step_to(double time) {solver.step_to(obj, time);}

  T obj;
  Solver<T> solver;
};

}
}

#endif
