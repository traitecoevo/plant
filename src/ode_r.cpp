#include "ode_r.h"

namespace ode {

OdeR::OdeR(SEXP fun, SEXP env, SEXP pars)
  : fun(fun),
    env(env),
    pars(pars),
    solver(this) {
}

size_t OdeR::size() const {
  return solver.get_size();
}

void OdeR::derivs(double time, iter_const y, iter dydt) {
  SEXP y_r;
  PROTECT(y_r = Rf_allocVector(REALSXP, (int)size()));
  std::copy(y, y + (int)size(), REAL(y_r));

  // 2. Slot for output
  SEXP dydt_r;

  // 3. Compute derivatives
  PROTECT(dydt_r = target(time, y_r));

  // 4. Copy derivatives from R vector into iterator
  std::copy(REAL(dydt_r), REAL(dydt_r) + size(), dydt);

  // 5. Cleanup
  UNPROTECT(2);
}

void OdeR::set_ode_state(std::vector<double> y, double time) {
  solver.reset(); // reset step sizes, etc.
  solver.set_state(y, time);
}

std::vector<double> OdeR::ode_state() const {
  return solver.get_state();
}

double OdeR::get_time() const {
  return solver.get_time();
}

std::vector<double> OdeR::get_times() const {
  return solver.get_times();
}

void OdeR::step() {
  solver.step();
}

void OdeR::step_fixed(double step_size) {
  solver.step_fixed(step_size);
}

void OdeR::step_to(double time) {
  solver.step_to(time);
}

void OdeR::advance(double time_max) {
  solver.advance(time_max);
}

void OdeR::advance_fixed(std::vector<double> times) {
  solver.advance_fixed(times);
}

// * R interface
std::vector<double> OdeR::r_derivs(double time, std::vector<double> y) {
  set_ode_state(y, time);
  std::vector<double> dydt(size());
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}

Rcpp::NumericMatrix OdeR::r_run(std::vector<double> times,
				      std::vector<double> y) {
  return solver.r_run(times, y);
}

// * Private methods
SEXP OdeR::target(double time, SEXP y) {
  return Rf_eval(Rf_lang4(fun, Rf_ScalarReal(time), y, pars), env);
}

}
