// -*-c++-*-
#ifndef PLANT_PLANT_CENSUS_GRADIENT_H_
#define PLANT_PLANT_CENSUS_GRADIENT_H_

#include <string>
#include <vector>

namespace plant {

// Why a metric has no gradient, and which species it was found on. A refusal is
// the only thing a number needs said beside it: a parameter no gradient exists
// for cannot be asked for, so it never reaches an answer, and every column that
// does reach one carries a number the sweep computed.
//
// The species is as fine as the grain gets: a row that could not be supplied is
// an intermediate of a recording spanning six stages and every cohort in them, so
// nothing below the species has a component to attribute it to.
struct refusal {
  std::string reason;
  // 1-based where a species' own strategy recorded the refusal, and -1 where
  // what failed spans every cohort in every stage and belongs to no species --
  // the reverse walk's own range refusal is of that kind. A CALLER MUST HANDLE
  // -1 RATHER THAN INDEXING WITH IT.
  int species = -1;

  bool happened() const { return !reason.empty(); }
};

// One census metric's gradient. A refused metric's whole row is not-a-number and
// `why` says what happened: a sum has no defined value with an undefined term, so
// refusal is metric-level and carries no localisation within a metric.
struct census_gradient {
  // metric-major, one column per registered parameter
  std::vector<std::vector<double>> gradient;
  std::vector<refusal> why;

  // What the sweep that produced the rows above did, carried beside them rather
  // than left on the solver. Beside, because a count read through its own
  // accessor is a count that can disagree with the refusal it belongs to -- and
  // did, for a day, while one of two refusal mechanisms zeroed it and the other
  // did not. They are written only where there is an answer to describe, so a
  // refused metric leaves them at the defaults instead of clearing back to them.

  // How many backward ranges were walked: one per width, plus one per requested
  // split that fell inside a range, which is what says a split cut.  Not
  // multiplied by the metric count -- one recording is swept once per metric, so
  // the ranges walked are the same however many rows are carried.
  size_t ranges = 0;

  // The adjoint the walk ended holding: one row per metric swept, one column per
  // entry of the first recorded state. It is d(census)/d(that state), which
  // census_initial_state_tangent computes forwards from the same state over the
  // same steps.
  std::vector<std::vector<double>> at_first_state;
};

}  // namespace plant

#endif
