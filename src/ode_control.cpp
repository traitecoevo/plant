#include <tree/ode_control.h>
#include <tree/util.h>

namespace tree {
namespace ode {

OdeControl::OdeControl() : OdeControl(1e-8, 1e-8, 1.0, 0.0,
				      1e-8, 10.0, 1e-6) {
}
OdeControl::OdeControl(double tol_abs_, double tol_rel_,
		       double a_y_, double a_dydt_,
		       double step_size_min_, double step_size_max_,
		       double step_size_initial_)
  : tol_abs(tol_abs_),
    tol_rel(tol_rel_),
    a_y(a_y_),
    a_dydt(a_dydt_),
    step_size_min(step_size_min_),
    step_size_max(step_size_max_),
    step_size_initial(step_size_initial_),
    last_step_size_shrank(false) {
}

double OdeControl::adjust_step_size(size_t dim, size_t ord,
				    double step_size,
				    const state_type& y,
				    const state_type& yerr,
				    const state_type& dydt) {
  double rmax = std::numeric_limits<double>::min();
  const double S = 0.9;

  for (size_t i = 0; i < dim; i++) {
    const double D0 = errlevel(y[i], dydt[i], step_size);
    const double r = std::abs(yerr[i]) / std::abs(D0);
    rmax = std::max(r, rmax);
  }

  if (rmax > 1.1) {
    // decrease step, no more than factor of 5, but a fraction S more
    // than scaling suggests (for better accuracy).
    double r = S / pow(rmax, 1.0 / ord);
    if (r < 0.2) {
      r = 0.2;
    }
    step_size *= r;
    last_step_size_shrank = true;
    if (step_size < step_size_min) {
      step_size = step_size_min;
      util::stop("Step size became too small");
    }
  } else if (rmax < 0.5) {
    // increase step, no more than factor of 5
    double r = S / pow (rmax, 1.0 / (ord + 1.0));
    if (r > 5.0) {
      r = 5.0;
    } else if (r < 1.0) {// Don't allow any decrease caused by S<1
      r = 1.0;
    }
    step_size *= r;
    if (step_size > step_size_max) {
      step_size = step_size_max;
    }
    last_step_size_shrank = false;
  } else {
    // otherwise no change to the size
    last_step_size_shrank = false;
  }

  return step_size;
}

double OdeControl::errlevel(double y, double dydt, double h) const {
  const double errlev = tol_rel * (a_y    * std::abs(y       )  +
				   a_dydt * std::abs(h * dydt)) + tol_abs;
  if (errlev <= 0.0) {
    util::stop("errlev <= zero");
  }
  return errlev;
}

bool OdeControl::step_size_shrank() const {
  return last_step_size_shrank;
}

}
}
