#include "util.h"
#include "Lorenz.h"

namespace ode {

namespace test {

Lorenz::Lorenz(double sigma, double R, double b)
  : sigma(sigma),
    R(R),
    b(b),
    solver(this) {
}

size_t Lorenz::size() const {
  return ode_dimension;
}

void Lorenz::derivs(double time, 
		    std::vector<double>::const_iterator y,
		    std::vector<double>::iterator dydt) {
  const double y0 = *y++;
  const double y1 = *y++;
  const double y2 = *y++;
  
  *dydt++ = sigma * ( y1 - y0 );
  *dydt++ = R * y0 - y1 - y0 * y2;
  *dydt++ = -b * y2 + y0 * y1;
}

void Lorenz::set_ode_state(std::vector<double> y, double t) {
  solver.reset();
  solver.set_state(y, t);
}

std::vector<double> Lorenz::ode_state() const {
  return solver.get_state();
}

double Lorenz::get_time() const {
  return solver.get_time();
}

std::vector<double> Lorenz::get_times() const {
  return solver.get_times();
}

void Lorenz::step() {
  solver.step();
}

void Lorenz::step_fixed(double step_size) {
  solver.step_fixed(step_size);
}

void Lorenz::step_to(double time) {
  solver.step_to(time);
}

void Lorenz::advance(double time_max) {
  solver.advance(time_max);
}

// This is something that can move up to a different level, as it's
// just straight up boilerplate.
std::vector<double> Lorenz::r_derivs(double time,
				     std::vector<double> y) {
  util::check_length(y.size(), size());
  std::vector<double> dydt(size());
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}

Rcpp::NumericMatrix Lorenz::r_run(std::vector<double> times,
				  std::vector<double> y) {
  return solver.r_run(times, y);
}

}

}
