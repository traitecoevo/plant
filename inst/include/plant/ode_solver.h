// -*-c++-*-
#ifndef PLANT_PLANT_ODE_SOLVER_H_
#define PLANT_PLANT_ODE_SOLVER_H_

#include <plant/ode_interface.h>
#include <plant/ode_control.h>
#include <plant/ode_step.h>
#include <plant/util.h>

#include <limits>
#include <vector>
#include <cstddef>

namespace plant {
namespace ode {

template <class System>
class Solver {
public:
  Solver(const System& system, OdeControl control_);
  void reset(const System& system);
  void set_state_from_system(const System& system);

  state_type get_state() const {return y;}
  double get_time() const {return time;}
  std::vector<double> get_times() const {return prev_times;}

  void advance(System& system, double time_max_);
  void advance_fixed(System& system, const std::vector<double>& times);

  void step(System& system);
  void step_to(System& system, double time_max_);

  void set_time_max(double time_max_);

private:
  void resize(size_t size_);
  void setup_dydt_in(System& system);
  void save_dydt_out_as_in();
  void set_time(double t);

  OdeControl control;
  Step<System> stepper;

  double step_size_last; // Size of last successful step (or suggestion)

  double time;     // Current time
  double time_max; // Time we will not go past
  std::vector<double> prev_times; // Vector of previous times.

  state_type y;        // Vector of current system state
  state_type yerr;     // Vector of error estimates
  state_type dydt_in;  // Vector of dydt at beginning of step
  state_type dydt_out; // Vector of dydt during step

  // NOTE: Ideas around this may change.
  bool dydt_in_is_clean;
};

// NOTE I'm setting the initial system size to 0 here, but some
// systems are self-initialising.
template <class System>
Solver<System>::Solver(const System& system, OdeControl control_)
  : control(control_) {
  reset(system);
}

// NOTE: This resets *everything* to basically a recreated object.
template <class System>
void Solver<System>::reset(const System& system) {
  prev_times.clear();
  step_size_last = control.step_size_initial;
  time_max = std::numeric_limits<double>::infinity();
  set_state_from_system(system);
}

template <class System>
void Solver<System>::set_state_from_system(const System& system) {
  set_time(ode::ode_time(system));
  resize(system.ode_size());
  system.ode_state(y.begin());
  system.ode_rates(dydt_in.begin());
  dydt_in_is_clean = true;
}

template <class System>
void Solver<System>::advance(System& system, double time_max_) {
  set_time_max(time_max_);
  while (time < time_max) {
    std::cout << "###### Solver advancing to next integration node ####### \n";
    step(system);
    std::cout << "End of step\n\n" << std::endl;
  }
}

// NOTE: We take a vector of times {t_0, t_1, ...}.  This vector
// *must* contain a starting time, but can otherwise be empty.  We
// will step exactly to t_1, then to t_2 up to the end point.  No step
// size adjustments will be done.  This is used in the SCM.
//
// NOTE: Careful here: exact floating point comparison in determining
// that we're starting from the right place.  However, because we take
// care to return and add end points exactly, this should actually be
// the correct move.
template <class System>
void Solver<System>::advance_fixed(System& system,
                                   const std::vector<double>& times) {
  if (times.empty()) {
    util::stop("'times' must be vector of at least length 1");
  }
  std::vector<double>::const_iterator t = times.begin();
  if (!util::identical(*t++, time)) {
    util::stop("First element in 'times' must be same as current time");
  }
  while (t != times.end()) {
    step_to(system, *t++);
  }
}

// After `stepper.step()`, the GSL checks to see if the step succeeded
// (some steppers look like they fail for non-user function error),
// and the divides the step size by 2.  If it fails with `EFAULT` or
// `EBADFUNC`, then it aborts.  The only place that errors are
// actually checked in the user function, and the two errors that
// cause abort are the only two that should be thrown there.
//
// There are several different logical step sizes:
//
// 1. this->step_size_last: Size of the last successful step last
//    time, or a suggestion of one.  This will get updated as leave
//    the function only if (1) the step is successful and (2) if we're
//    not in the final step.  It's not actually quite the size of the
//    last step, either -- it's the size that the controller suggested
//    updating the step size too after the last current step.
//
// 2. step_size: The size that the current iteration actually advanced
//    the system (or will) via `stepper.step`.
//
// 3. step_size_next: The size of the proposed next step (or retry of
//    the current step).
template <class System>
void Solver<System>::step(System& system) {
  const double time_orig = time, time_remaining = time_max - time;
  double step_size = step_size_last;

  // Save y in case of failure in a step (recall that stepper.step
  // changes 'y')
  const state_type y_orig = y;
  const size_t size = y.size();

  // Compute the derivatives at the beginning.
  setup_dydt_in(system);

  while (true) {
    // Does this appear to be the last step before reaching `time_max`?
    const bool final_step = step_size > time_remaining;
    if (final_step) {
      step_size = time_remaining;
    }
    stepper.step(system, time, step_size, y, yerr, dydt_in, dydt_out);

    const double step_size_next =
      control.adjust_step_size(size, stepper.order(), step_size,
			       y, yerr, dydt_out);

    if (control.step_size_shrank()) {
        // GSL checks that the step size is actually decreased.
        // Probably we can do this by comparing against hmin?  There are
        // probably loops that this will not catch, but require that
        // hmin << t
         const double time_next = time + step_size_next;
      if (step_size_next < step_size && time_next > time_orig) {
      	// Step was decreased. Undo step (resetting the state y and
        // time), and try again with the new step_size.
      	y         = y_orig;
      	time      = time_orig;
      	step_size = step_size_next;
      } else {
      	// We've reached limits of machine accuracy in differences of
      	// step sizes or time (or both).
      	util::stop("Cannot achive the desired accuracy");
      }
    } else {
      // We have successfully taken a step and will return.  Update
      // time to reflect this, ensuring that if we're on the last step
      // we will end up exactly at time_max.
      //
      // Suggest step size for next time-step. Change of step size is not
      //  suggested in the final step, because that step can be very
      //  small compared to previous step, to reach time_max.
      if (final_step) {
	      time = time_max;
      } else {
	      time += step_size;
	      step_size_last = step_size_next;
      }
      prev_times.push_back(time);
      save_dydt_out_as_in();
      return; // This exits the infinite loop.
    }
  }
}

// This takes a step up to time "time_max_", regardless of what the
// integration error says.  This is used by advance_fixed
template <class System>
void Solver<System>::step_to(System& system, double time_max_) {
  set_time_max(time_max_);
  setup_dydt_in(system);
  stepper.step(system, time, time_max - time, y, yerr, dydt_in, dydt_out);
  save_dydt_out_as_in();

  time = time_max;
  prev_times.push_back(time);
}

template <class System>
void Solver<System>::resize(size_t size_) {
  y.resize(size_);
  yerr.resize(size_);
  dydt_in.resize(size_);
  dydt_out.resize(size_);
  stepper.resize(size_);
}

// TODO: Need to go through and sort out the logic here; it might be
// that we can just access the rates.  Could probably do something
// where we ask what 'y' is in the model, and if it's the same then we
// use the existing rates, otherwise we apply our y.  For now, this
// should be correct, but comes at a cost of an extra evaluation.
template <class System>
void Solver<System>::setup_dydt_in(System& system) {
  if (stepper.can_use_dydt_in && !dydt_in_is_clean) {
    // TODO: Not clear that this is the right thing here; should just
    // be able to look up the correct dydt rates because we've already
    // set state?
    //   system.ode_rates(dydt_in.begin());
    ode::derivs(system, y, dydt_in, time, 0);
    dydt_in_is_clean = true;
  }
}

template <class System>
void Solver<System>::save_dydt_out_as_in() {
  if (stepper.first_same_as_last) {
    dydt_in = dydt_out;
    dydt_in_is_clean = true;
  } else {
    dydt_in_is_clean = false;
  }
}

template <typename System>
void Solver<System>::set_time(double t) {
  const int ulp = 2; // units in the last place (accuracy)
  if (prev_times.size() > 0 &&
      !util::almost_equal(prev_times.back(), t, ulp)) {
    util::stop("Time does not match previous (delta = " +
	       util::to_string(prev_times.back() - t) +
	       "). Reset solver first.");
  }
  time = t;
  if (prev_times.empty()) { // only if first time (avoids duplicate times)
    prev_times.push_back(time);
  }
}

template <class System>
void Solver<System>::set_time_max(double time_max_) {
  if (!util::is_finite(time_max_)) {
    util::stop("time_max must be finite!");
  }
  if (time_max_ < time) {
    util::stop("time_max must be greater than (or equal to) current time");
  }
  time_max = time_max_;
}

}
}

#endif
