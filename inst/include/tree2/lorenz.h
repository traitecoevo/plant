// -*-c++-*-
#ifndef TREE_LORENZ_H_
#define TREE_LORENZ_H_

#include <tree2/ode_interface.h>

namespace ode {

namespace test {

// This is meant to be a simple test case for the ODE solver.  We
// mimick a class that has some nontrivial internal state (this
// doesn't of course, but we pretend that it does).
//
// Probably what I should test at some point is that a *collection*
// of these works as expected.
class Lorenz {
public:
  Lorenz(double sigma_, double R_, double b_)
    : sigma(sigma_), R(R_), b(b_),
      y0(0.0), y1(0.0), y2(0.0),
      dy0dt(0.0), dy1dt(0.0), dy2dt(0.0) {
  }

  // ODE interface.
  size_t ode_size() const {return ode_dimension;}
  ode::const_iterator set_ode_state(ode::const_iterator it) {
    y0 = *it++;
    y1 = *it++;
    y2 = *it++;
    dy0dt = sigma * (y1 - y0);
    dy1dt = R * y0 - y1 - y0 * y2;
    dy2dt = -b * y2 + y0 * y1;
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    *it++ = y0;
    *it++ = y1;
    *it++ = y2;
    return it;
  }
  ode::iterator ode_rates(ode::iterator it) const {
    *it++ = dy0dt;
    *it++ = dy1dt;
    *it++ = dy2dt;
    return it;
  }
  std::vector<double> pars() const {
    std::vector<double> ret;
    ret.push_back(sigma);
    ret.push_back(R);
    ret.push_back(b);
    return ret;
  }

private:
  static const int ode_dimension = 3;
  double sigma, R, b;
  double y0, y1, y2;
  double dy0dt, dy1dt, dy2dt;
};

}

}

#endif
