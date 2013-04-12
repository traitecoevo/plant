// -*-c++-*-
#ifndef TREE_PLANT_SPLINE_
#define TREE_PLANT_SPLINE_

#include "plant.h"
#include "spline.h"

namespace model {

// The name here may change, as this functions more like a *Strategy*
// than a Plant.  Or a bit like a Species.  ContinuousPlant might be
// better and reflect the fact that this is not part of an inheritance
// heirarcy.

// Also note that this does not inherit from ode::OdeTarget, and so
// cannot be stepped.  It does not have state itself, though the
// plants that it happens to contain do.
class PlantSpline {
public:
  PlantSpline(Strategy s, int n_plants);
  PlantSpline(Strategy *s, int n_plants);

  // Copy constructor, assigment and destructor (rule of three)
  PlantSpline(const PlantSpline &other);
  PlantSpline& operator=(PlantSpline rhs);
  ~PlantSpline();
  friend void swap(PlantSpline &a, PlantSpline &b);

  void compute_vars_phys(spline::Spline *env);

  void r_compute_vars_phys(spline::Spline env);
  Rcpp::List r_get_plants() const;
  Rcpp::NumericMatrix r_get_ode_values() const;
  Rcpp::List r_get_ode_values_approx() const;

private:
  void build(int n_plants);
  size_t ode_size() const;

  bool standalone;
  Strategy *strategy;
  std::vector<double> mass_leaf_log;
  std::vector<Plant> plants;

  std::vector< std::vector<double> > ode_values;
  std::vector< spline::Spline > ode_values_approx;
};

}

#endif
