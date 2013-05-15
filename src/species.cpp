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

template <>
ode::iter Species<CohortTop>::ode_rates(ode::iter it) const {
  // This is the base case.
  for ( plants_const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_rates(it);

  // Then tweak the boundary conditions.
  const int n = plants.back().ode_size();
  std::fill(it - n, it, 0.0);

  return it;
}

template <>
void Species<CohortTop>::initialise() {
  add_seeds(1);

  // TODO: Confirm that this is actually needed.
  //
  // Create a blank light environment (so that light value is 1 the
  // whole way down)
  const int blank_env_size = 5;
  spline::Spline env;
  for (int i = 0; i < blank_env_size; i++)
    env.add_point(height_max() * i, 1.0);
  env.init_self();
  // And run the plants in this.  We need this so that the growth rate
  // is correct.
  compute_vars_phys(&env);
}

}
