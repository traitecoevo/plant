// -*-c++-*-
#ifndef TREE_SEED_RAIN_H_
#define TREE_SEED_RAIN_H_

#include <vector>

namespace model {

// This is going to track a bit of information about the seed rain.
// Not sure if it will remain a real independent class, as much of
// what it does could probably be done better with an iterator.  If
// that is the case, we'll remove it and replace with a proper wrapper
// around an iterator.

class SeedRain {
public:
  SeedRain(int n);
  SeedRain(std::vector<double> seed_rain);

  double operator()() const;

  void first();
  void next();

  void set(std::vector<double> x);

private:
  std::vector<double> seed_rain;
  size_t curr_idx;
};

}

#endif
