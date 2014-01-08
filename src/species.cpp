#include "species.h"

#include "cohort_discrete.h"

// Specialisation for CohortDiscrete, where we add a cohort of the
// appropriate size at the back of the vector.
namespace model {

SpeciesBase::~SpeciesBase() {
}

// This is a bit of a hassle, but a full specialisation will be
// compiled directly so cannot go into the header file (or it breaks
// the One Definition Rule).  We could get around this by inlining,
// but not sure if that is any better.  These functions seem a bit big
// to inline, but I don't know.
template <>
void Species<CohortDiscrete>::add_seeds(int n) {
  if (n > 0) {
    CohortDiscrete p = seed;
    p.set_n_individuals(n);
    plants.push_back(p);
  }
}

template <>
int Species<CohortDiscrete>::r_n_individuals() const {
  int n = 0;
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); ++it)
    n += it->get_n_individuals();
  return n;
}

// NOTE: This is quite a bit trickier than the Plant case.  We have to
// integrate over the end points of the distribution, counting a
// non-existant seed as the left-most point.  So, with at least one
// cohort in the population, the integral is defined.
//
// NOTE: The CohortTop::leaf_area_above() calculation takes into
// account the different densities of different cohorts, so that is
// automatically included in this calculation.  It's possible that it
// makes more sense for the density to be included *within* the
// calculation here, however?
//
// NOTE: This is simply performing numerical integration, via the
// trapezium rule, of the leaf_area_above with respect to plant
// height.  You'd think that this would be nicer to do in terms of a
// call to an external trapezium integration function, but building
// and discarding the intermediate storage ends up being a nontrivial
// cost.  A more general iterator version might be possible, but with
// the fiddliness around the boundary conditions that won't likely be
// useful.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights.  The integral is not defined, but an empty light
// environment seems appropriate.
//
// NOTE: A similar early-exit condition to the Plant version is used;
// once the lower bound of the trazpeium is zero, we stop including
// individuals.  Working with the boundary cohort is tricky here,
// because we might need to include that, too: always in the case of a
// single cohort (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <>
double Species<CohortTop>::leaf_area_above(double height) const {
  if (size() == 0 || height_max() < height)
    return 0.0;
  double tot = 0.0;
  plants_const_iterator it = plants.begin();
  double h1 = it->height(), f_h1 = it->leaf_area_above(height);

  for (++it; it != plants.end(); ++it) {
    const double h0 = it->height(), f_h0 = it->leaf_area_above(height);
    if (!util::is_finite(f_h0))
      Rcpp::stop("Detected non-finite contribution");
    tot += (h1 - h0) * (f_h1 + f_h0);
    // Upper point moves for next time:
    h1   = h0;
    f_h1 = f_h0;
    if (h0 < height)
      break;
  }

  if (size() == 1 || f_h1 > 0) {
    const double h0 = seed.height(), f_h0 = seed.leaf_area_above(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
  }

  return tot / 2;
}

template <>
void Species<CohortTop>::compute_vars_phys(const Environment& environment) {
  if (control().plant_assimilation_approximate_use)
    compute_assimilation_spline(environment);
  for (plants_iterator it = plants.begin();
       it != plants.end(); ++it)
    it->compute_vars_phys(environment);
  seed.compute_initial_conditions(environment);
}

template <>
Rcpp::NumericMatrix Species<CohortTop>::r_get_state() const {
  const size_t nc = size() + 1, nr = seed.state_size();
  // Extra +1 here for the seed.
  state tmp(nc * nr);
  state::iterator it = get_state(tmp.begin());
  seed.get_state(it);

  Rcpp::NumericMatrix ret(static_cast<int>(nr),
			  static_cast<int>(nc));
  std::copy(tmp.begin(), tmp.end(), ret.begin());
  Rcpp::CharacterVector rownames =
    Rcpp::CharacterVector::create("height", "log.mortality",
				  "seeds", "log.density",
				  "pr.survival.birth");
  ret.attr("dimnames") = Rcpp::List::create(rownames, R_NilValue);

  return ret;
}

template <>
void Species<CohortTop>::r_set_state(Rcpp::NumericMatrix x) {
  util::check_dimensions(static_cast<size_t>(x.nrow()),
			 static_cast<size_t>(x.ncol()),
			 seed.state_size(), size() + 1); // + seed
  state tmp(x.begin(), x.end());
  state::const_iterator it = set_state(tmp.begin());
  seed.set_state(it);
}

// TODO: this does not check that the resulting state is sensible (in
// particular the monotonic check).  That should probably be added.
template <>
void Species<CohortTop>::r_force_state(Rcpp::NumericMatrix x) {
  if (x.ncol() == 0)
    Rcpp::stop("Boundary condition for seeds must always be given");
  const size_t n = static_cast<size_t>(x.ncol() - 1);
  plants.clear();
  for (size_t i = 0; i < n; ++i)
    add_seeds(1);
  r_set_state(x);
}

SEXP species(Rcpp::CppClass individual, Strategy s) {
  std::string individual_type =
    util::rcpp_class_demangle(Rcpp::as<std::string>(individual));
  SEXP ret = R_NilValue;
  if (individual_type == "Plant") {
    Species<Plant> obj(s);
    ret = Rcpp::wrap(obj);
  } else if (individual_type == "CohortDiscrete") {
    Species<CohortDiscrete> obj(s);
    ret = Rcpp::wrap(obj);
  } else if (individual_type == "CohortTop") {
    Species<CohortTop> obj(s);
    ret = Rcpp::wrap(obj);
  } else {
    Rcpp::stop("Cannot make Species of " + individual_type);
  }
  return ret;
}

}
