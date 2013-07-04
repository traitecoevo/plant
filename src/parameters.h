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

// NOTE: The parameter n_patches is stored as a double for use with
// the lookup class, but will be converted to an int on use.
// Non-integer values will be accepted silently but this conversion
// will happen.  Not ideal, but not the end of the world.
//
// TODO: See better solution in Control.
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

  // Additional control parameters:
  const Control& get_control() const;
  void set_control(Control x);

  // Data -- public for now (see github issue #17).
  Disturbance disturbance_regime;
  double c_ext;      // Light extinction coefficient
  double patch_area; // Size of the patch (m^2)
  double Pi_0;       // Probability of survival during dispersal
  size_t n_patches;  // Number of patches in the metacommunity
  std::vector<Strategy> strategies;

  // Algorithm control.
  Control control;

private:
  void do_build_lookup();
  void set_parameters_post_hook();
  void reset();

  double _n_patches;

  typedef std::vector<Strategy>::iterator strategy_iterator;
  typedef std::vector<Strategy>::const_iterator strategy_const_iterator;
};

}

RCPP_EXPOSED_CLASS(model::Parameters)

#endif
