#include "ROde.h"

ROde::ROde(SEXP fun, SEXP env, SEXP pars) :
  fun(fun), env(env), pars(pars), solver(this) {
}

void ROde::derivs(double time,
		  std::vector<double>::const_iterator y,
		  std::vector<double>::iterator dydt) {
  SEXP y_r;
  PROTECT(y_r = Rf_allocVector(REALSXP, size()));
  std::copy(y, y + size(), REAL(y_r));

  // 2. Slot for output
  SEXP dydt_r;

  // 3. Compute derivatives
  PROTECT(dydt_r = target(time, y_r));

  // 4. Copy derivatives from R vector into iterator
  std::copy(REAL(dydt_r), REAL(dydt_r) + size(), dydt);

  // 5. Cleanup
  UNPROTECT(2);
}

std::vector<double> ROde::r_derivs(double time, std::vector<double> y) {
  size_ = y.size();
  std::vector<double> dydt(size());
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}

SEXP ROde::target(double time, SEXP y) {
  return Rf_eval(Rf_lang4(fun, Rf_ScalarReal(time), y, pars), env);
}
