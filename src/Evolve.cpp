#include "Evolve.h"
#include <cmath>     // fabs

Evolve::Evolve() {
  reset();
}

void Evolve::set_state(std::vector<double> y_, double t_) {
  y    = y_;
  time = t_;
  resize(y.size());
  reset();
}

// TODO: In contrast with GSL, this would make more sense as
// "step()".
// 
// After `s.step()`, the GSL checks to see if the step succeeded
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
//    the system (or will) via `s.step`.
// 
// 3. step_size_next: The size of the proposed next step (or retry of
//    the current step).

void Evolve::step() {
  const double time_orig = time, time_remaining = time_max - time;
  double step_size = step_size_last;
  // Does this appear to be the last step before reaching `time_max`?
  bool final_step = false;
  
  // Save y in case of failure in a step (recall that s.step changes
  // 'y')
  const std::vector<double> y_orig = y;

  // Compute the derivatives at the beginning.
  s.derivs(time, y.begin(), dydt_in.begin());

  while ( true ) {
    // TODO: This allows negative step direction.  Worth it?
    final_step = ((time_remaining >= 0.0 && step_size > time_remaining) || 
		  (time_remaining <  0.0 && step_size < time_remaining));
    if ( final_step )
      step_size = time_remaining;
    s.step(time, step_size, y, yerr, dydt_in, dydt_out);

    count++;
    const double 
      step_size_next = c.adjust_step_size(size, s.order(), step_size,
					  y, yerr, dydt_out);
    
    if ( c.step_size_shrank() ) {
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

      return;
    }
  }
}

// This sets time_max, but I need to be careful about other things
// using it...
void Evolve::advance(double time_max_) {
  if ( time_max_ <= time )
    Rf_error("time_max must be greater than time");
  time_max = time_max_;
  while ( time < time_max )
    step();
}

Rcpp::NumericMatrix Evolve::r_run(std::vector<double> times, 
				  std::vector<double> y_) {
  Rcpp::NumericMatrix ret(size, times.size()-1);
  Rcpp::NumericMatrix::iterator out = ret.begin();

  std::vector<double>::iterator t = times.begin();
  set_state(y_, *t++); // This makes 'y' contains mutable state.

  while ( t != times.end() ) {
    advance(*t++);
    std::copy(y.begin(), y.end(), out);
    out += size;
  }

  return ret;
}


// Note that this could push us past time_max!
void Evolve::step_fixed(double step_size) {
  // Save y in case of failure in the step
  std::vector<double> y0 = y;

  // Compute the derivatives at the beginning.
  s.derivs(time, y.begin(), dydt_in.begin());

  s.step(time, step_size, y, yerr, dydt_in, dydt_out);
  
  const double step_size_next =
    c.adjust_step_size(size, s.order(), step_size, y, yerr, dydt_out);

  if ( c.step_size_shrank() ) {
    failed_steps++;
    y = y0;
  } else {
    count++;
    step_size_last = step_size_next;
    time += step_size;
  }
}

std::vector<double> Evolve::r_derivs() {
  std::vector<double> dydt(size);
  s.derivs(time, y.begin(), dydt.begin());
  return dydt;
}

void Evolve::reset() {
  count = 0;
  failed_steps = 0;
  step_size_last = 1e-6; // See ode.md
  time_max = DBL_MAX;    // or infinite?
  s.reset();
}

void Evolve::resize(int size_) {
  size = size_;
  yerr.resize(size);
  dydt_in.resize(size);
  dydt_out.resize(size);
  s.resize(size);
}
