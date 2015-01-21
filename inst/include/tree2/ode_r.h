// -*-c++-*-
#ifndef TREE_ODE_R_H_
#define TREE_ODE_R_H_

#include <tree2/ode_interface.h>
#include <tree2/util.h>
#include <Rcpp.h>

namespace ode {
namespace test {

// This is not meant to be fast; it's for checking purposes only.
// There's lots of extra copies here, and lots of R->C interface, but
// it should be very polymorphic.
//
// Requirements are a function derivs(y, t) -> dydt
//
// We simulate stored state in much the same way as Lorenz does.
class OdeR {
public:
  OdeR(Rcpp::Function derivs_, Rcpp::Function state_, double time_)
    : derivs(derivs_), state(state_), time(time_) {
    update_state();
  }
  void update_state() {
    state_type res = Rcpp::as<state_type>(state());
    if (y.size() != res.size()) {
      y.resize(res.size());
      dydt.resize(res.size());
    }
    std::copy(res.begin(), res.end(), y.begin());
    update_dydt();
  }
  double ode_size() const {return y.size();}
  double ode_time() const {return time;}
  ode::const_iterator set_ode_values(ode::const_iterator it, double time_) {
    std::copy_n(it, ode_size(), y.begin());
    time = time_;
    update_dydt();
    return it + ode_size();
  }
  ode::iterator ode_values(ode::iterator it) const {
    return std::copy_n(y.begin(), ode_size(), it);
  }
  ode::iterator ode_rates(ode::iterator it) const {
    return std::copy_n(dydt.begin(), ode_size(), it);
  }

private:
  void update_dydt() {
    state_type res = Rcpp::as<state_type>(derivs(y, time));
    util::check_length(res.size(), ode_size());
    std::copy(res.begin(), res.end(), dydt.begin());
  }
  Rcpp::Function derivs;
  Rcpp::Function state;
  state_type y;
  state_type dydt;
  double time;
};

}
}

#endif
