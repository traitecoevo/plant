// -*-c++-*-
#ifndef TREE_ODE_SOLVER_H_
#define TREE_ODE_SOLVER_H_

#include <Rcpp.h>
#include <vector>
#include <cmath> // fabs

#include "ode_control.h"
#include "ode_step.h"

namespace ode {

template <class Problem>
class Solver {
public:
  Solver(Problem *problem);

  void set_state(std::vector<double> y_, double t_);
  std::vector<double> get_state() const;
  double get_time() const;
  std::vector<double> get_times() const;

  void set_state_from_problem(double time);

  void step();
  void step_fixed(double step_size);
  void advance(double time_max_);

  void step_to(double time_max_);

  void reset();

  size_t get_size() const;
  
  Rcpp::NumericMatrix r_run(std::vector<double> times, 
			    std::vector<double> y_);
  
private:
  Problem *problem;

  void resize(size_t size_);

  size_t size;           // Problem dimension
  int count;             // Number of steps since reset
  int failed_steps;      // Number of failed steps since reset
  double step_size_last; // Size of last successful step (or suggestion)

  double time;     // Current time
  double time_max; // Time we will not go past
  std::vector<double> times; // Vector of previous times.

  std::vector<double> y;        // Vector of current problem state
  std::vector<double> yerr;     // Vector of error estimates
  std::vector<double> dydt_in;  // Vector of dydt at beginning of step
  std::vector<double> dydt_out; // Vector of dydt during step

  // Control parameters
  double step_size_min, step_size_max;
  int    no_steps_max;

  // Used internally.
  OdeControl control;
  Step<Problem> stepper;
};

// NOTE I'm setting the initial problem size to 0 here.  In fact, an
// object of class OdeTarget could quite happily initialse completely
// here.
template <class Problem>
Solver<Problem>::Solver(Problem *problem)
  : problem(problem),
    size(0),
    stepper(problem) {
  reset();
}

template <class Problem>
void Solver<Problem>::set_state(std::vector<double> y_, double t_) {
  times.clear();
  y    = y_;
  time = t_;
  times.push_back(time);
  resize(y.size());
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
  return times;
}

// TODO: Not sure, but we might be able to use get_time here too;
// think about if that should be added to the OdeTarget class?
template <class Problem>
void Solver<Problem>::set_state_from_problem(double time) {
  std::vector<double> y(problem->ode_size());
  problem->ode_values(y.begin());
  set_state(y, time);
}

// TODO: In contrast with GSL, this would make more sense as
// "step()".
// 
// After `stepper.step()`, the GSL checks to see if the step succeeded
// (some steppers look like they fail for non-user function error),
// and the divides the step size by 2.  If it fails with `EFAULT` or
// `EBADFUNC`, then it aborts.  The only place that errors are
// actually checked in the user function, and the two errors that
// cause abort are the only two that should be thrown there.
// 
// TODO: What are we planning on doing with count?  failed_steps?
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
  // Does this appear to be the last step before reaching `time_max`?
  bool final_step = false;
  
  // Save y in case of failure in a step (recall that stepper.step
  // changes 'y')
  const std::vector<double> y_orig = y;

  // Compute the derivatives at the beginning.
  problem->derivs(time, y.begin(), dydt_in.begin());

  while ( true ) {
    // TODO: This allows negative step direction.  Worth it?
    final_step = ((time_remaining >= 0.0 && step_size > time_remaining) || 
		  (time_remaining <  0.0 && step_size < time_remaining));
    if ( final_step )
      step_size = time_remaining;
    stepper.step(time, step_size, y, yerr, dydt_in, dydt_out);

    count++;
    const double step_size_next = 
      control.adjust_step_size(size, stepper.order(), step_size,
			       y, yerr, dydt_out);
    
    if ( control.step_size_shrank() ) {
      // GSL checks that the step size is actually decreased.
      // Probably we can do this by comparing against hmin?  There are
      // probably loops that this will not catch, but require that
      // hmin << t
      const double time_next = time + step_size_next;
      if ( fabs(step_size_next) < fabs(step_size) && 
	   time_next != time_orig ) {
	// Step was decreased. Undo step (resetting the state y and
	// time), and try again with the new step_size.
	y         = y_orig;
	time      = time_orig; // no-op right now.
	step_size = step_size_next;
	failed_steps++;
      } else {
	// We've reached limits of machine accuracy in differences of
	// step sizes or time (or both).
	Rf_error("Cannot achive the desired accuracy");
      }
    } else {
      // We have successfully taken a step and will return.  Update
      // time to reflect this, ensuring that if we're on the last step
      // we will end up exactly at time_max.
      //
      // Suggest step size for next time-step. Change of step size is not
      //  suggested in the final step, because that step can be very
      //  small compared to previous step, to reach time_max. 
      if ( final_step )  {
      	time = time_max;
      } else {
      	time += step_size;
      	step_size_last = step_size_next;
      }
      times.push_back(time);

      return;
    }
  }
}

// This sets time_max, but I need to be careful about other things
// using it...
template <class Problem>
void Solver<Problem>::advance(double time_max_) {
  if ( time_max_ < time )
    Rf_error("time_max must be greater than (or equal to) current time");
  time_max = time_max_;
  while ( time < time_max )
    step();
}

template <class Problem>
void Solver<Problem>::step_to(double time_max_) {
  time_max = time_max_;
  problem->derivs(time, y.begin(), dydt_in.begin());
  stepper.step(time, time_max - time, y, yerr, dydt_in, dydt_out);
  count++;
  time = time_max;
  times.push_back(time);
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

  Rcpp::NumericMatrix ret((int)size, (int)times.size()-1);
  Rcpp::NumericMatrix::iterator out = ret.begin();

  while ( t != times.end() ) {
    advance(*t++);
    std::copy(y.begin(), y.end(), out);
    out += size;
  }

  return ret;
}


// Note that this could push us past time_max!
template <class Problem>
void Solver<Problem>::step_fixed(double step_size) {
  // Save y in case of failure in the step
  std::vector<double> y0 = y;

  // Compute the derivatives at the beginning.
  problem->derivs(time, y.begin(), dydt_in.begin());

  stepper.step(time, step_size, y, yerr, dydt_in, dydt_out);
  
  const double step_size_next =
    control.adjust_step_size(size, stepper.order(), step_size, 
			     y, yerr, dydt_out);

  if ( control.step_size_shrank() ) {
    failed_steps++;
    y = y0;
  } else {
    count++;
    step_size_last = step_size_next;
    time += step_size;
    times.push_back(time);
  }
}

// NOTE: This resets *everything* to basically a recreated object.
//
// TODO: Work out how to deal with control parameters here.
template <class Problem>
void Solver<Problem>::reset() {
  time = 0;
  times.clear();
  count = 0;
  failed_steps = 0;
  step_size_last = 1e-6; // See ode.md
  time_max = DBL_MAX;    // or infinite?
  size = 0;
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

}

#endif
