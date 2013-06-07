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
    p.set_n_individuals(n);
    plants.push_back(p);
  }
}

template <>
int Species<CohortDiscrete>::r_n_individuals() const {
  int n = 0;
  for ( plants_const_iterator it = plants.begin();
	it != plants.end(); it++ )
    n += it->get_n_individuals();
  return n;
}

// NOTE: This is quite a bit trickier than the Plant case.  We have to
// integrate over the end points of the distribution, counting a
// non-existant seed as the left-most point.  So, with at least one
// cohort in the population, the integral is defined.
//
// NOTE: In contrast with the Plant version, we stop *after* including
// a y point that is zero, which might be the boundary cohort.  In
// addition, we always need to include at least two points.
//
// NOTE: This is further complicated by the fact that plants are
// stored largest to smallest.  If we used a vector to store x/y
// values we'd have to push_back(), and the resulting integral would
// be *negative* (because the x values would be decreasing).  Using a
// list allows pushing to the front.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights.  The integral is not defined, but an empty light
// environment seems appropriate.
template <>
double Species<CohortTop>::leaf_area_above(double height) const {
  if (size() == 0 || height_max() < height)
    return(0.0);
  std::list<double> x, y;
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); it++) {
    x.push_front(it->height());
    y.push_front(it->leaf_area_above(height));
    if (y.front() == 0)
      break;
  }
  if (y.front() > 0 || y.size() < 2) {
    x.push_front(seed.height());
    y.push_front(seed.leaf_area_above(height));
  }
  return util::trapezium(x, y);
}

template <>
void Species<CohortTop>::compute_vars_phys(const Environment& environment) {
  for ( plants_iterator it = plants.begin();
	it != plants.end(); it++ )
    it->compute_vars_phys(environment);
  seed.compute_initial_conditions(environment);
}

}
