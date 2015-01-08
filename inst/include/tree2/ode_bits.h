// -*-c++-*-
#ifndef TREE_ODE_BITS_H_
#define TREE_ODE_BITS_H_

#include <vector>
#include <memory>

// #include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
// #include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
// #include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
// #include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
// #include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

#include <boost/numeric/odeint.hpp>
#include <RcppCommon.h>

namespace ode {

// The problem, and why it doesn't fit very well with most ode
// integrators.
//
// We have a system that we want to push around with the usual ode
// approaches -- it's fully defined by a set of y values and given
// these y values we can compute dy/dt to compute the future state of
// the system by numerical integration.
//
// However, there are also nontrivial auxillary variables 'z' that are
// dependent on the 'y' values.  These are important for doing things
// like changing the dimension of the system later on, so I do need to
// get them right.
//
// The other problem is that in general we can't assume that the final
// correct values will get written into the object.  That is not the
// case for some things I've tried here, so the values within an
// object are not a safe proxy for a complete system.  I think the
// deal here is that we extrapolate for the final step.  If that's the
// case then dropping full state is both safest *and* most efficient.
// So that's nice.

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

// TODO: Should this actually store a copy or should it store a reference?
// Might be best to have this be a subclass of OdeSystem and work with
// a reference, I think.
template <typename T>
class OdeRunner {
public:
  OdeRunner(T obj_) : obj(obj_) {}
  // For odeint.
  void operator()(const std::vector<double>& y,
                  std::vector<double>& dydt,
                  const double /* t */) {
    obj.set_ode_values(y);
    obj.ode_rates(dydt);
  }
  T obj;
};

template <typename T>
struct state_saver {
  typedef T state_type;
  struct observer {
    // Note this is a vector of vectors
    std::vector<state_type>& y;
    std::vector<double>& t;
    observer(std::vector<state_type>& y_, std::vector<double>& t_)
      : y(y_), t(t_) {}
    void operator()(const state_type &yi, double ti) {
      y.push_back(yi);
      t.push_back(ti);
    }
  };
  state_saver() : steps(0), obs(y, t) {}
  std::vector<state_type> y;
  std::vector<double> t;
  size_t steps;
  observer obs;
};

template <typename T>
class OdeSystem {
  typedef std::vector<double> state_type;
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type>
  stepper_basic_type;
  typedef boost::numeric::odeint::controlled_runge_kutta<stepper_basic_type>
  stepper_type;

public:
  // TODO: In my old ode system, there are additional coefficients
  // here:
  // a_y and a_dydt, plus step_size_min, step_size_max
  // Going to need to shepherd those in too I think.  Looks like we
  // have the same basic algorithm here -- see
  // http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/steppers.html (search "classical")
  OdeSystem(T obj_, double abs_tol, double rel_tol) :
    stepper_controlled(boost::numeric::odeint::make_controlled<stepper_basic_type>(abs_tol, rel_tol)),
    runner(OdeRunner<T>(obj_)),
    y(obj_.ode_values()) {
  }

  // Stepping:
  void do_step(double dt) {
    stepper_basic.do_step(runner, y, t, dt);
  }
  bool try_step(double dt) {
    return stepper_controlled.try_step(runner, y, t, dt) == 0;
  }
  // Now, what do we do with dt here?  Previously what I did was to
  // keep track of the last step size.  The initial step size was
  // arbitrarily set at 1e-6.
  void advance(double t, double dt) {
    using boost::numeric::odeint::integrate_adaptive;
    integrate_adaptive(stepper_controlled, runner, y, 0.0, t, dt);
  }
  state_saver<state_type> advance_save(double t, double dt) {
    state_saver<state_type> state;
    state.obs(y, t);
    using boost::numeric::odeint::integrate_adaptive;
    integrate_adaptive(stepper_controlled, runner, y, 0.0, t, dt, state.obs);
    return state;
  }

  // This is potentially harmful as this is not guaranteed to be
  // anywhere near up-to-date.
  // Might be better to work with a reference here perhaps?
  // Or optionally update as we go?
  T get_obj() { return runner.obj; }

  stepper_basic_type stepper_basic;
  stepper_type stepper_controlled;
  OdeRunner<T> runner;
  double t; // What is this doing here?
  state_type y;
};

}

#endif
