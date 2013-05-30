#include "seed_rain.h"

#include <Rcpp.h>
#include "util.h"

namespace model {

// Empty seed rain -- all species initialised to have zero seed rain.
SeedRain::SeedRain(size_t n_species)
  : seed_rain(n_species, 0.0),
    curr_idx(0) {
}

size_t SeedRain::size() const {
  return seed_rain.size();
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

std::vector<double> SeedRain::get_seed_rain() const {
  return seed_rain;
}
void SeedRain::set_seed_rain(std::vector<double> x) {
  util::check_length(x.size(), size());
  seed_rain = x;
}

SeedRain seed_rain(std::vector<double> x) {
  SeedRain rain(x.size());
  rain.set_seed_rain(x);
  return rain;
}

}
