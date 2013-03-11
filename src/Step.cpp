#include "Step.h"

Step::Step() {
  // Fixed for the lorenz attractor right now.
  derivs = &derivs_lorenz;
}

void Step::resize(int size_) {
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
void Step::step(double time, double step_size,
		std::vector<double> &y,
		std::vector<double> &yerr,
		const std::vector<double> &dydt_in,
		std::vector<double> &dydt_out) {
  const double h = step_size; // Historical reasons.

  // k1 step:
  std::copy(dydt_in.begin(), dydt_in.end(), k1.begin());

  for ( int i = 0; i < size; i++ )
    ytmp[i] = y[i] + b21 * h * k1[i];

  // k2 step:
  derivs(time + ah[0] * h, ytmp.begin(), k2.begin());

  for ( int i = 0; i < size; i++ )
    ytmp[i] = y[i] + h * (b3[0] * k1[i] + b3[1] * k2[i]);

  // k3 step:
  derivs(time + ah[1] * h, ytmp.begin(), k3.begin());

  for ( int i = 0; i < size; i++ )
    ytmp[i] = y[i] + h * (b4[0] * k1[i] + b4[1] * k2[i] + b4[2] * k3[i]);

  // k4 step:
  derivs(time + ah[2] * h, ytmp.begin(), k4.begin());

  for ( int i = 0; i < size; i++ )
    ytmp[i] =
      y[i] + h * (b5[0] * k1[i] + b5[1] * k2[i] + b5[2] * k3[i] +
                  b5[3] * k4[i]);

  // k5 step
  derivs(time + ah[3] * h, ytmp.begin(), k5.begin());

  for ( int i = 0; i < size; i++ )
    ytmp[i] =
      y[i] + h * (b6[0] * k1[i] + b6[1] * k2[i] + b6[2] * k3[i] +
                  b6[3] * k4[i] + b6[4] * k5[i]);

  // k6 step and final sum
  derivs(time + ah[4] * h, ytmp.begin(), k6.begin());
  
  for ( int i = 0; i < size; i++ ) {
    // GSL does this in two steps, but not sure why.
    const double d_i = c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i];
    y[i] += h * d_i;
  }

  // Evaluate dydt_out.
  derivs(time + h, y.begin(), dydt_out.begin());

  // Difference between 4th and 5th order, for error calculations
  for ( int i = 0; i < size; i++ )
    yerr[i] = h * (ec[1] * k1[i] + ec[3] * k3[i] + ec[4] * k4[i] +
		   ec[5] * k5[i] + ec[6] * k6[i]);
}

// Skipping this actually looks totally benign, as these are never
// read within `step` without first being written to.
void Step::reset() {
  std::fill(k1.begin(),   k1.end(),   0.0);
  std::fill(k2.begin(),   k2.end(),   0.0);
  std::fill(k3.begin(),   k3.end(),   0.0);
  std::fill(k4.begin(),   k4.end(),   0.0);
  std::fill(k5.begin(),   k5.end(),   0.0);
  std::fill(k6.begin(),   k6.end(),   0.0);
  std::fill(ytmp.begin(), ytmp.end(), 0.0);
}

unsigned int Step::order() {
  // In GSL, comment says "FIXME: should this be 4?"
  return 5;
}

// RKCK coefficients, from GSL
const double Step::ah[] = {
  1.0 / 5.0, 0.3, 3.0 / 5.0, 1.0, 7.0 / 8.0 };

const double Step::b21 = 1.0 / 5.0;
const double Step::b3[] = { 3.0 / 40.0, 9.0 / 40.0 };
const double Step::b4[] = { 0.3, -0.9, 1.2 };
const double Step::b5[] = { 
  -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0 };

const double Step::b6[] = { 
  1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 
  44275.0 / 110592.0, 253.0 / 4096.0 };

const double Step::c1 = 37.0 / 378.0;
const double Step::c3 = 250.0 / 621.0;
const double Step::c4 = 125.0 / 594.0;
const double Step::c6 = 512.0 / 1771.0;

const double Step::ec[] = { 
  0.0, 37.0 / 378.0 - 2825.0 / 27648.0, 0.0,
  250.0 / 621.0 - 18575.0 / 48384.0,
  125.0 / 594.0 - 13525.0 / 55296.0,
  -277.0 / 14336.0, 512.0 / 1771.0 - 0.25 };

void derivs_lorenz(double t,
		   std::vector<double>::const_iterator y,
		   std::vector<double>::iterator dydt) {
  const double 
    sigma = 10.0,
    R = 28.0,
    b = 8.0 / 3.0;

  const double y0 = *y++;
  const double y1 = *y++;
  const double y2 = *y++;

  *dydt++ = sigma * ( y1 - y0 );
  *dydt++ = R * y0 - y1 - y0 * y2;
  *dydt++ = -b * y2 + y0 * y1;
}


