// -*-c++-*-

#include <vector>

class Step {
public:
  Step();
  void resize(int size_);
  void step(double time, double step_size,
	    std::vector<double> &y,
	    std::vector<double> &yerr,
	    const std::vector<double> &dydt_in,
	    std::vector<double> &dydt_out);
  void reset();
  unsigned int order();

  typedef void (*derivs_fun)(double time,
			     std::vector<double>::const_iterator y,
			     std::vector<double>::iterator dydt);
  derivs_fun derivs;

private:
  // Intermediate storage, representing state (was GSL rkck_state_t)
  double size;
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

// This is the derivatives function for the Lorenz system.
void derivs_lorenz(double t,
		   std::vector<double>::const_iterator y,
		   std::vector<double>::iterator dydt);
