// -*-c++-*-
#ifndef TREE_ODE_SOLVER_H_
#define TREE_ODE_SOLVER_H_

#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include <cmath> // fabs

#include "ode_control.h"
#include "ode_step.h"
#include "util.h"

namespace ode {

template <class Problem>
class Solver {
public:
  Solver(Problem *problem_);
  Solver(Problem *problem_, OdeControl control_);

  void set_state(std::vector<double> y_, double t_);
  std::vector<double> get_state() const;
  double get_time() const;
  std::vector<double> get_times() const;

  void set_state_from_problem();

  void step();
  void step_fixed(double step_size);
  void advance(double time_max_);

  void step_to(double time_max_);
  void advance_fixed(std::vector<double> times);

  void reset();

  size_t get_size() const;
  
  Rcpp::NumericMatrix r_run(std::vector<double> times, 
			    std::vector<double> y_);
  OdeControl r_control() const;
  
private:
  void resize(size_t size_);
  void setup_dydt_in();
  void save_dydt_out_as_in();
  void set_time_max(double time_max_);

  Problem *problem;

  size_t size;           // Problem dimension
  int count;             // Number of steps since reset
  int failed_steps;      // Number of failed steps since reset
  double step_size_last; // Size of last successful step (or suggestion)

  double time;     // Current time
  double time_max; // Time we will not go past
  std::vector<double> prev_times; // Vector of previous times.

  std::vector<double> y;        // Vector of current problem state
  std::vector<double> yerr;     // Vector of error estimates
  std::vector<double> dydt_in;  // Vector of dydt at beginning of step
  std::vector<double> dydt_out; // Vector of dydt during step

  // NOTE: Ideas around this may change.
  bool dydt_in_is_clean;

  OdeControl control;
  Step<Problem> stepper;
};

// NOTE I'm setting the initial problem size to 0 here.  In fact, an
// object of class OdeTarget could quite happily initialse completely
// here.
template <class Problem>
Solver<Problem>::Solver(Problem *problem_)
  : problem(problem_),
    size(0),
    stepper(problem_) {
  reset();
}

template <class Problem>
Solver<Problem>::Solver(Problem *problem_, OdeControl control_)
  : problem(problem_),
    size(0),
    control(control_),
    stepper(problem) {
  reset();
}

// The error condition here is to ensure that the prev_times vector
// makes sense.  If the incoming time is different to the last time
// that we were at, the problem should be reset before continuing
// (which will clear `prev_times`).  Not sure if this is the right way
// to work with this.
template <class Problem>
void Solver<Problem>::set_state(std::vector<double> y_, double t_) {
  const bool accuracy = 2; // d.p. of accuracy - see util.h
  if (prev_times.size() > 0 &&
      !util::almost_equal(prev_times.back(), t_, accuracy))
    ::Rf_error("Time does not match previous (delta=%2.5e). Reset solver first.", prev_times.back() - t_);
  y    = y_;
  time = t_;
  if (prev_times.empty()) // only if first time (avoids duplicate times)
    prev_times.push_back(time);
  resize(y.size());
  dydt_in_is_clean = false;
}

template <class Problem>
std::vector<double> Solver<Problem>::get_state() const {
  return y;
}

template <class Problem>
double Solver<Problem>::get_time() const {
  return time;
}

template <class Problem>
std::vector<double> Solver<Problem>::get_times() const {
  return prev_times;
}

template <class Problem>
void Solver<Problem>::set_state_from_problem() {
  std::vector<double> vars(problem->ode_size());
  problem->ode_values(vars.begin());
  set_state(vars, problem->get_time());
  problem->ode_rates(dydt_in.begin());
  dydt_in_is_clean = true;
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

template <class Problem>
void Solver<Problem>::step() {
  const double time_orig = time, time_remaining = time_max - time;
  double step_size = step_size_last;

  // Save y in case of failure in a step (recall that stepper.step
  // changes 'y')
  const std::vector<double> y_orig = y;

  // Compute the derivatives at the beginning.
  setup_dydt_in();

  while (true) {
    // Does this appear to be the last step before reaching `time_max`?
    const bool final_step = step_size > time_remaining;
    if (final_step)
      step_size = time_remaining;
    stepper.step(time, step_size, y, yerr, dydt_in, dydt_out);

    count++;
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
	time      = time_orig; // no-op right now.
	step_size = step_size_next;
	failed_steps++;
      } else {
	// We've reached limits of machine accuracy in differences of
	// step sizes or time (or both).
	::Rf_error("Cannot achive the desired accuracy");
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

      return;
    }
  }
}

// This sets time_max, but I need to be careful about other things
// using it...
template <class Problem>
void Solver<Problem>::advance(double time_max_) {
  if (time_max_ < time)
    ::Rf_error("time_max must be greater than (or equal to) current time");
  set_time_max(time_max_);
  while (time < time_max)
    step();
}

// This takes a step up to time "time_max_", regardless of what the
// error says.  There is duplication here with step_fixed()?
template <class Problem>
void Solver<Problem>::step_to(double time_max_) {
  set_time_max(time_max_);
  setup_dydt_in();
  stepper.step(time, time_max - time, y, yerr, dydt_in, dydt_out);
  save_dydt_out_as_in();

  count++;
  time = time_max;
  prev_times.push_back(time);
}

// NOTE: We take a vector of times {t_0, t_1, ...}.  This vector
// *must* contain a starting time, but can otherwise be empty.  We
// will step exactly to t_1, then to t_2 up to the end point.  No step
// size adjustments will be done.  This is used in the EBT.
//
// NOTE: Careful here: exact floating point comparison in determining
// that we're starting from the right place.  However, because we take
// care to return and add end points exactly, this should actually be
// the correct move.
template <class Problem>
void Solver<Problem>::advance_fixed(std::vector<double> times) {
  if (times.size() < 1)
    ::Rf_error("'times' must be vector of at least length 1");
  std::vector<double>::const_iterator t = times.begin();
  if (!util::identical(*t++, times.front()))
    ::Rf_error("First element in 'times' must be same as current time");
  while (t != times.end())
    step_to(*t++);
}

template <class Problem>
Rcpp::NumericMatrix Solver<Problem>::r_run(std::vector<double> times, 
					   std::vector<double> y_) {
  std::vector<double>::iterator t = times.begin();

  // This sets the initial step size to be small.
  reset();

  // This makes `y` contains mutable state, and `size` contain current
  // problem dimension.
  set_state(y_, *t++);

  Rcpp::NumericMatrix ret(static_cast<int>(size),
			  static_cast<int>(times.size())-1);
  Rcpp::NumericMatrix::iterator out = ret.begin();

  while (t != times.end()) {
    advance(*t++);
    std::copy(y.begin(), y.end(), out);
    out += size;
  }

  return ret;
}

template <class Problem>
OdeControl Solver<Problem>::r_control() const {
  return control;
}

template <class Problem>
void Solver<Problem>::step_fixed(double step_size) {
  if (time + step_size > time_max)
    ::Rf_error("step would push us past time_max");

  // Save y in case of failure in the step
  std::vector<double> y0 = y;

  // Compute the derivatives at the beginning.
  setup_dydt_in();
  stepper.step(time, step_size, y, yerr, dydt_in, dydt_out);
  
  const double step_size_next =
    control.adjust_step_size(size, stepper.order(), step_size, 
			     y, yerr, dydt_out);

  if (control.step_size_shrank()) {
    failed_steps++;
    y = y0;
  } else {
    count++;
    step_size_last = step_size_next;
    time += step_size;
    prev_times.push_back(time);
    save_dydt_out_as_in();
  }
}

// NOTE: This resets *everything* to basically a recreated object.
template <class Problem>
void Solver<Problem>::reset() {
  time = 0;
  prev_times.clear();
  count = 0;
  failed_steps = 0;
  step_size_last = 1e-6; // See ode.md
  time_max = R_PosInf;   // but never allow setting to this value
  size = 0;
  dydt_in_is_clean = false;
  resize(size);
  stepper.reset();
}

template <class Problem>
size_t Solver<Problem>::get_size() const {
  return size;
}

template <class Problem>
void Solver<Problem>::resize(size_t size_) {
  size = size_;
  yerr.resize(size);
  dydt_in.resize(size);
  dydt_out.resize(size);
  stepper.resize(size);
}

template <class Problem>
void Solver<Problem>::setup_dydt_in() {
  if (stepper.can_use_dydt_in() && !dydt_in_is_clean) {
    problem->derivs(time, y.begin(), dydt_in.begin());
    dydt_in_is_clean = true;
  }
}

template <class Problem>
void Solver<Problem>::save_dydt_out_as_in() {
  if (stepper.first_same_as_last()) {
    dydt_in = dydt_out;
    dydt_in_is_clean = true;
  } else {
    dydt_in_is_clean = false;
  }
}

template <class Problem>
void Solver<Problem>::set_time_max(double time_max_) {
  if (!util::is_finite(time_max_))
    ::Rf_error("time_max must be finite!");
  time_max = time_max_;
}

}

#endif
