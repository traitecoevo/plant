// -*-c++-*-
#ifndef TREE_PARAMETERS_
#define TREE_PARAMETERS_

#include <Rcpp.h>

#include "lookup.h"
#include "strategy.h"

// NOTE: Most of the interface decisions here are around interaction
// with R via Rcpp; written in plain C++ this is a really simple class
// with a vector of strategies and a couple of simulation-wide
// parameters.  As such, it might change quite a bit as the model
// evolves.

namespace model {

class Parameters : public util::Lookup {
public:
  // Construction
  Parameters();
  void reset();

  // Querying
  size_t size() const {return strategies.size();}
  Strategy r_get_strategy(size_t idx);
  Rcpp::List r_get_strategies();

  // Setting
  void add_strategy(Strategy s);

  // TODO: Support for deletion of strategies, but think about
  // implications in things that use them.

  // Data
  double mean_disturbance_interval;
  double c_ext;      // Light extinction coefficient
  double patch_area; // Size of the patch (m^2, I think?) [TODO]
  double Pi_0;       // Probability of survival during dispersal
  std::vector<Strategy> strategies;
private:
  void do_build_lookup();
};

}

RCPP_EXPOSED_CLASS(model::Parameters)

#endif
