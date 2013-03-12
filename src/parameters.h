// -*-c++-*-
#ifndef TREE_PARAMETERS_
#define TREE_PARAMETERS_

#include <Rcpp.h>

#include "strategy.h"

namespace model {

class Parameters {
public:
  // Construction
  Parameters();
  void reset();

  // Querying
  unsigned int size() const {return strategies.size();}
  Rcpp::List get_params() const;
  Rcpp::List get_strategy(int idx) const;
  Rcpp::List get_strategies() const;

  // Setting
  void set_params(Rcpp::List x);
  void add_strategy(Rcpp::List x);
  void set_strategy(Rcpp::List x, int idx);

  // Data
  double mean_disturbance_interval;
  double c_ext; // Light extinction coefficient
  std::vector<Strategy> strategies;
};

}

#endif
