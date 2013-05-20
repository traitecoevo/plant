#include "species.h"

#include "cohort_discrete.h"

// Specialisation for CohortDiscrete, where we add a cohort of the
// appropriate size at the back of the vector.
namespace model {

// This is a bit of a hassle, but a full specialisation will be
// compiled directly so cannot go into the header file (or it breaks
// the One Definition Rule).  We could get around this by inlining,
// but not sure if that is any better.  These functions seem a bit big
// to inline, but I don't know.
template <>
void Species<CohortDiscrete>::add_seeds(int n) {
  if ( n > 0 ) {
    CohortDiscrete p = seed;
    p.r_set_n_individuals(n); // TODO: Temporary
    plants.push_back(p);
  }
}

template <>
int Species<CohortDiscrete>::r_n_individuals() const {
  int n = 0;
  for ( plants_const_iterator it = plants.begin();
	it != plants.end(); it++ )
    n += it->r_n_individuals();
  return n;
}

template<>
double Species<CohortTop>::leaf_area_above(double height) const {
  std::vector<double> x, y;
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); it++ ) {
    const double a = it->leaf_area_above(height);
    if (a > 0) {
      // TODO: Here, it would be nice to abstract away the size
      // dimension, rather than use mass_leaf directly.
      x.push_back(it->mass_leaf());
      y.push_back(a);
    } else {
      break;
    }
  }

  return util::trapezium(x, y);
}

template <>
void Species<CohortTop>::compute_vars_phys(const Environment& environment) {
  for ( plants_iterator it = plants.begin();
	it != plants.end(); it++ )
    it->compute_vars_phys(environment);
  seed.compute_vars_phys(environment);
}

}
