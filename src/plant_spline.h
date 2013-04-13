// -*-c++-*-
#ifndef TREE_PLANT_SPLINE_
#define TREE_PLANT_SPLINE_

#include "plant.h"
#include "spline.h"       // light environment
#include "multi_spline.h" // approximate plant
#include "ode_target.h"   // iter_const

namespace model {

// The name here may change, as this functions more like a *Strategy*
// than a Plant.  Or a bit like a Species.  ContinuousPlant might be
// better and reflect the fact that this is not part of an inheritance
// hierarchy.

// Also note that this does not inherit from ode::OdeTarget, and so
// cannot be stepped.  It does not have state itself, though the
// plants that it happens to contain do.

// Note that not all aspects of plants are gettable from this; in
// particular, I've focussed on *only* exposing the ODE interface so
// that the underlying plant can change independently of this.
// However, this assumes that plants are identifiable by a single
// axis; any second variable (including location, but also impact of
// history, resprouting, etc) will violate this assumption.

class PlantSpline {
public:
  PlantSpline(Strategy  s, double mass_leaf_max, int n_plants);
  PlantSpline(Strategy *s, double mass_leaf_max, int n_plants);

  // Copy constructor, assigment and destructor (rule of three)
  PlantSpline(const PlantSpline &other);
  PlantSpline& operator=(PlantSpline rhs);
  ~PlantSpline();
  friend void swap(PlantSpline &a, PlantSpline &b);

  double mass_leaf_max() const;
  void compute_vars_phys(spline::Spline *env);
  ode::iter ode_rates(double m, ode::iter it) const;

  void r_compute_vars_phys(spline::Spline env);
  std::vector<double> r_ode_rates(double m) const;
  Rcpp::List r_get_plants() const;
  spline::MultiSpline r_get_plants_approx() const;

private:
  size_t ode_size() const;
  void initialise(double mass_leaf_max, int n_plants);
  void build_plants_approx();

  bool standalone;
  Strategy *strategy;
  Plant seed;
  std::vector<double> mass_leaf;
  std::vector<Plant> plants;
  spline::MultiSpline plants_approx;
};

}

#endif
