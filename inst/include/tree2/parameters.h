// -*-c++-*-
#ifndef TREE2_PARAMETERS_H_
#define TREE2_PARAMETERS_H_

#include <tree2/control.h>
#include <tree2/strategy.h>
#include <tree2/disturbance.h>
#include <vector>

// TODO: I will possibly move out the "Patch" parameters out into
// their own simple list class at some point, to make this a bit more
// coherent.
//
// TODO: Will require some free functions on the R side:
//   * add_strategy (with flag for mutant/non mutant)

namespace tree2 {

struct Parameters {
  Parameters();
  // Data -- public for now (see github issue #17).
  double c_ext;      // Light extinction coefficient
  double patch_area; // Size of the patch (m^2)
  double Pi_0;       // Probability of survival during dispersal
  size_t n_patches;  // Number of patches in the metacommunity
  std::vector<Strategy> strategies;
  std::vector<double> seed_rain;
  std::vector<bool>   is_resident;

  // Disturbance regime.
  Disturbance disturbance;

  // Algorithm control.
  Control control;

  // Default strategy.
  Strategy strategy_default;

  // Some little query functions for use on the C side:
  size_t size() const;
  size_t n_residents() const;
  size_t n_mutants() const;
  bool validate() const;
};

}

#endif
