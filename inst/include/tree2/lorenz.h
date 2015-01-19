// -*-c++-*-
#ifndef TREE_LORENZ_H_
#define TREE_LORENZ_H_

#include <Rcpp.h>
#include <tree2/ode_bits.h>
#include <tree2/util.h> // check_length

namespace ode {

namespace test {

// What I'm going to do is to conform to the same requirements that
// rodeint uses, but without the parameters requirement (at least for
// now).
//
//   derivs: (const const T& y, T& dydt, const double t)
//
// The complication I had before was that the nested nature of the
// problem caused some issues: so we want to take a bunch of state
// variables and disbuse them across the system in a really weird
// way.
//
// This is why having the iterator bits was very important.

class Lorenz {
public:
  Lorenz(double sigma_, double R_, double b_)
    : sigma(sigma_), R(R_), b(b_), state(ode_dimension) {}
  static size_t size() { return ode_dimension; }
  
  // These are the big ones:
  std::vector<double> ode_values() const { return state; }
  void set_ode_values(const std::vector<double>& y) {
    util::check_length(y.size(), size());
    state = y;
  }
  void ode_rates(std::vector<double>& dydt) {
    const double y0 = state[0];
    const double y1 = state[1];
    const double y2 = state[2];
    dydt[0] = sigma * (y1 - y0);
    dydt[1] = R * y0 - y1 - y0 * y2;
    dydt[2] = -b * y2 + y0 * y1;
  }

  // Then some stuff for the R interface:
  Rcpp::NumericVector r_get_pars() {
    using namespace Rcpp;
    return NumericVector::create(_["sigma"] = sigma,
                                 _["R"] = R,
                                 _["b"] = b);
  }
  std::vector<double> r_ode_rates() {
    std::vector<double> dydt(ode_dimension);
    ode_rates(dydt);
    return dydt;
  }

private:
  static const int ode_dimension = 3;
  double sigma, R, b;
  std::vector<double> state;
};

}

}

#endif
