#include "seed_rain.h"

#include <Rcpp.h>
#include "util.h"

namespace model {

// Empty seed rain
SeedRain::SeedRain()
  : curr_idx(0) {
}

// Unless given an initial vector, in which case start with that.
SeedRain::SeedRain(std::vector<double> seed_rain)
  : seed_rain(seed_rain),
    curr_idx(0) {
}

double SeedRain::operator()() const {
  if (curr_idx >= seed_rain.size())
    ::Rf_error("Requested rain out of bounds");
  return seed_rain[curr_idx];
}

void SeedRain::first() {
  curr_idx = 0;
}
void SeedRain::next() {
  curr_idx++;
}

void SeedRain::set(std::vector<double> x) {
  util::check_length(x.size(), seed_rain.size());
  seed_rain = x;
}

}
