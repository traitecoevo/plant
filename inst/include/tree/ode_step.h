// -*-c++-*-

#ifndef TREE_ODE_STEP_H_
#define TREE_ODE_STEP_H_

#include <vector>
#include <cstddef>
#include "ode_target.h"

namespace ode {

template <class Problem>
class Step {
public:
  Step(Problem *problem_);
  void resize(size_t size_);
  void step(double time, double step_size,
	    std::vector<double> &y,
	    std::vector<double> &yerr,
	    const std::vector<double> &dydt_in,
	    std::vector<double> &dydt_out);
  void reset();
  size_t order() const;
  bool can_use_dydt_in() const;
  bool first_same_as_last() const;
  void derivs(double time, iterator_const y, iterator dydt);

  // Private soon, but not right now (will be once these are never
  // instantiated directly).
  Problem *problem;

private:
  // Intermediate storage, representing state (was GSL rkck_state_t)
  size_t size;
  std::vector<double> k1, k2, k3, k4, k5, k6, ytmp;

  // Cash carp constants, from GSL.
  static const double ah[];
  static const double b21;
  static const double b3[];
  static const double b4[];
  static const double b5[];
  static const double b6[];
  static const double c1;
  static const double c3;
  static const double c4;
  static const double c6;

  // These are the differences of fifth and fourth order coefficients
  // for error estimation
  static const double ec[];
};

template <class Problem>
Step<Problem>::Step(Problem *problem_) : problem(problem_) {
}

template <class Problem>
void Step<Problem>::resize(size_t size_) {
  size = size_;
  k1.resize(size);
  k2.resize(size);
  k3.resize(size);
  k4.resize(size);
  k5.resize(size);
  k6.resize(size);
  ytmp.resize(size);
}

// Think carefully about ownership of data, draw a diagram, and go
// from there.
template <class Problem>
void Step<Problem>::step(double time, double step_size,
			 std::vector<double> &y,
			 std::vector<double> &yerr,
			 const std::vector<double> &dydt_in,
			 std::vector<double> &dydt_out) {
  const double h = step_size; // Historical reasons.

  // k1 step:
  std::copy(dydt_in.begin(), dydt_in.end(), k1.begin());

  for (size_t i = 0; i < size; ++i)
    ytmp[i] = y[i] + b21 * h * k1[i];

  // k2 step:
  derivs(time + ah[0] * h, ytmp.begin(), k2.begin());

  for (size_t i = 0; i < size; ++i)
    ytmp[i] = y[i] + h * (b3[0] * k1[i] + b3[1] * k2[i]);

  // k3 step:
  derivs(time + ah[1] * h, ytmp.begin(), k3.begin());

  for (size_t i = 0; i < size; ++i)
    ytmp[i] = y[i] + h * (b4[0] * k1[i] + b4[1] * k2[i] + b4[2] * k3[i]);

  // k4 step:
  derivs(time + ah[2] * h, ytmp.begin(), k4.begin());

  for (size_t i = 0; i < size; ++i)
    ytmp[i] =
      y[i] + h * (b5[0] * k1[i] + b5[1] * k2[i] + b5[2] * k3[i] +
                  b5[3] * k4[i]);

  // k5 step
  derivs(time + ah[3] * h, ytmp.begin(), k5.begin());

  for (size_t i = 0; i < size; ++i)
    ytmp[i] =
      y[i] + h * (b6[0] * k1[i] + b6[1] * k2[i] + b6[2] * k3[i] +
                  b6[3] * k4[i] + b6[4] * k5[i]);

  // k6 step and final sum
  derivs(time + ah[4] * h, ytmp.begin(), k6.begin());

  for (size_t i = 0; i < size; ++i) {
    // GSL does this in two steps, but not sure why.
    const double d_i = c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i];
    y[i] += h * d_i;
  }

  // Evaluate dydt_out.
  derivs(time + h, y.begin(), dydt_out.begin());

  // Difference between 4th and 5th order, for error calculations
  for (size_t i = 0; i < size; ++i)
    yerr[i] = h * (ec[1] * k1[i] + ec[3] * k3[i] + ec[4] * k4[i] +
		   ec[5] * k5[i] + ec[6] * k6[i]);
}

// Skipping this actually looks totally benign, as these are never
// read within `step` without first being written to.
template <class Problem>
void Step<Problem>::reset() {
  std::fill(k1.begin(),   k1.end(),   0.0);
  std::fill(k2.begin(),   k2.end(),   0.0);
  std::fill(k3.begin(),   k3.end(),   0.0);
  std::fill(k4.begin(),   k4.end(),   0.0);
  std::fill(k5.begin(),   k5.end(),   0.0);
  std::fill(k6.begin(),   k6.end(),   0.0);
  std::fill(ytmp.begin(), ytmp.end(), 0.0);
}

template <class Problem>
size_t Step<Problem>::order() const {
  // In GSL, comment says "FIXME: should this be 4?"
  return 5;
}

template <class Problem>
bool Step<Problem>::can_use_dydt_in() const {
  return true;
}

// This implies can_use_dydt_in -- so this is first_same_as_last AND
// can_use_dydt_in.
template <class Problem>
bool Step<Problem>::first_same_as_last() const {
  return true;
}

template <class Problem>
void Step<Problem>::derivs(double time, iterator_const y, iterator dydt) {
  problem->derivs(time, y, dydt);
}

// RKCK coefficients, from GSL
template <class Problem>
const double Step<Problem>::ah[] = {
  1.0 / 5.0, 0.3, 3.0 / 5.0, 1.0, 7.0 / 8.0 };

template <class Problem>
const double Step<Problem>::b21 = 1.0 / 5.0;
template <class Problem>
const double Step<Problem>::b3[] = { 3.0 / 40.0, 9.0 / 40.0 };
template <class Problem>
const double Step<Problem>::b4[] = { 0.3, -0.9, 1.2 };
template <class Problem>
const double Step<Problem>::b5[] = {
  -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0 };

template <class Problem>
const double Step<Problem>::b6[] = {
  1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0,
  44275.0 / 110592.0, 253.0 / 4096.0 };

template <class Problem>
const double Step<Problem>::c1 = 37.0 / 378.0;
template <class Problem>
const double Step<Problem>::c3 = 250.0 / 621.0;
template <class Problem>
const double Step<Problem>::c4 = 125.0 / 594.0;
template <class Problem>
const double Step<Problem>::c6 = 512.0 / 1771.0;

template <class Problem>
const double Step<Problem>::ec[] = {
  0.0, 37.0 / 378.0 - 2825.0 / 27648.0, 0.0,
  250.0 / 621.0 - 18575.0 / 48384.0,
  125.0 / 594.0 - 13525.0 / 55296.0,
  -277.0 / 14336.0, 512.0 / 1771.0 - 0.25 };

}

#endif
