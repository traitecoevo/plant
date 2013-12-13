// -*-c++-*-
#ifndef TREE_PARAMETERS_H_
#define TREE_PARAMETERS_H_

#include <Rcpp.h>

#include "lookup.h"
#include "strategy.h"
#include "disturbance.h"
#include "control.h"

// NOTE: Most of the interface decisions here are around interaction
// with R via Rcpp; written in plain C++ this is a really simple class
// with a vector of strategies and a couple of simulation-wide
// parameters.  As such, it might change quite a bit as the model
// evolves.

namespace model {

class Parameters : public util::Lookup {
public:
  typedef util::PtrWrapper<Parameters> ptr;
  // Construction
  Parameters();

  // Querying
  size_t size() const;
  Strategy r_at(size_t idx);
  Rcpp::List r_get_strategies();

  // Setting
  void add_strategy(Strategy s);
  void add_strategy_mutant(Strategy s);

  // Additional control parameters:
  const Control& get_control() const;
  void set_control(Control x);
  void set_control_parameters(Rcpp::List x);

  // Data -- public for now (see github issue #17).
  double c_ext;      // Light extinction coefficient
  double patch_area; // Size of the patch (m^2)
  double Pi_0;       // Probability of survival during dispersal
  size_t n_patches;  // Number of patches in the metacommunity
  double mean_disturbance_interval; // Eponymous, in years.
  std::vector<Strategy> strategies;
  std::vector<double> seed_rain;
  std::vector<bool>   is_resident;

  // Algorithm control.
  Control control;

  // * R interface
  std::vector<double> r_seed_rain() const;
  void r_set_seed_rain(std::vector<double> x);

  std::vector<bool> r_is_resident() const;
  void r_set_is_resident(std::vector<bool> x);

private:
  void do_build_lookup();
  void set_parameters_post_hook();
  void reset();
  void push_control_to_strategies();

  double _n_patches;

  typedef std::vector<Strategy>::iterator strategy_iterator;
  typedef std::vector<Strategy>::const_iterator strategy_const_iterator;
};

}

RCPP_EXPOSED_CLASS(model::Parameters)

#endif
