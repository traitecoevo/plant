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

}
