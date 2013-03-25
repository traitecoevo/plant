// -*-c++-*-
#ifndef TREE_PARAMETERS_
#define TREE_PARAMETERS_

#include <Rcpp.h>

#include "lookup.h"
#include "strategy.h"

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
  void add_strategy(Rcpp::List x);
  void set_strategy(Rcpp::List x, int idx);

  // Data
  double mean_disturbance_interval;
  double c_ext;      // Light extinction coefficient
  double patch_area; // Size of the patch (m^2, I think?) [TODO]
  std::vector<Strategy> strategies;
private:
  void do_build_lookup();
};

}

#endif
