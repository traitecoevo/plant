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
  unsigned int size() const {return strategies.size();}
  Strategy get_strategy(int idx);
  Rcpp::List get_strategies();

  // Setting
  void add_strategy(Strategy s);
  void set_strategy(Rcpp::List x, int idx);

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
