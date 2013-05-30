// -*-c++-*-
#ifndef TREE_SEED_RAIN_H_
#define TREE_SEED_RAIN_H_

#include <vector>
#include <Rcpp.h>

namespace model {

// This is going to track a bit of information about the seed rain.
// Not sure if it will remain a real independent class, as much of
// what it does could probably be done better with an iterator.  If
// that is the case, we'll remove it and replace with a proper wrapper
// around an iterator.
//
// This really feels badly like a pseudoclass.  However, there are
// some unresolved issues around how I want to have the seed rain play
// out via the Environment.
//
// The constructor here is a bit of a hassle.  Basically -- when
// making a SeedRain object in C++, we'll specify the size at
// initialisation and then fill it with values later on.  But from R,
// it would be nice to give the *actual* values at the start, but
// we'll have to do this in a couple of steps.
class SeedRain {
public:
  SeedRain(size_t n_species);

  size_t size() const;
  double operator()() const;

  void first();
  void next();

  std::vector<double> get_seed_rain() const;
  void set_seed_rain(std::vector<double> x);

private:
  std::vector<double> seed_rain;
  size_t curr_idx;
};

SeedRain seed_rain(std::vector<double> x);

}

RCPP_EXPOSED_CLASS(model::SeedRain)

#endif
