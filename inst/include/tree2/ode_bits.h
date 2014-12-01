// -*-c++-*-
#ifndef TREE_ODE_BITS_H_
#define TREE_ODE_BITS_H_

#include <vector>

namespace ode {

// Here we are preferring composition over inheritence.
//
// The basic problem: we have some object that we can set the state of
// and the compute derivatives.  We want to use that.  This is totally
// contrived for things like the Lorenz system where the derivative
// follows directly from the state, but it's important for things like
// our forest patches where we need to set the state for everything
// and then we can look up the derivatives.
//
// We'll need some set_ode_values_incremental functions.  We need this
// for all three functions.  Just assume they work and template magic
// will take care of the rest!

template <typename T>
class OdeSystem {
  // typedef std::vector<double> state_type;
  // typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type>
  // stepper_basic_type;
  // typedef boost::numeric::odeint::controlled_runge_kutta<stepper_basic_type>
  // stepper_type;

public:
  OdeSystem(T obj_) : obj(obj_) {}
  // This is the form required by odeint:
  void derivs(const std::vector<double>& y,
              std::vector<double>& dydt,
              const double /* t */) {
    obj.set_ode_values(y);
    obj.ode_rates(dydt);
  }
  // For odeint:
  void operator() (const std::vector<double>& y, 
                   std::vector<double>& dydt,
                   const double t) {
    derivs(y, dydt, t);
  }
  T obj;
};

}

#endif

