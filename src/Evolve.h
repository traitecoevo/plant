// -*-c++-*-
#ifndef EVOLVE_H_
#define EVOLVE_H_

#include <Rcpp.h>
#include <vector>

#include "Control.h"
#include "Step.h"

#include "Lorenz.h"

class Evolve {
public:
  Evolve();
  ~Evolve();

  void set_state(std::vector<double> y_, double t_);
  std::vector<double> get_state() const { return y; }
  double get_time() const { return time; }

  void step();
  void step_fixed(double step_size);
  void advance(double time_max_);

  void reset();
  
  std::vector<double> r_derivs();
  Rcpp::NumericMatrix r_run(std::vector<double> times, 
			    std::vector<double> y_);
  
private:
  Lorenz *pr;

  void resize(int size_);

  int size;              // Problem dimension
  int count;             // Number of steps since reset
  int failed_steps;      // Number of failed steps since reset
  double step_size_last; // Size of last successful step (or suggestion)

  double time;     // Current time
  double time_max; // Time we will not go past

  std::vector<double> y;        // Vector of current problem state
  std::vector<double> yerr;     // Vector of error estimates
  std::vector<double> dydt_in;  // Vector of dydt at beginning of step
  std::vector<double> dydt_out; // Vector of dydt during step

  // Control parameters
  double step_size_min, step_size_max;
  int    no_steps_max;

  // Used internally.
  Control c;
  Step s;
};

#endif
