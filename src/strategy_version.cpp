#include <plant.h>

// Scientific version of a model, read from the `scientific_version` constant
// declared on each strategy class (see inst/include/plant/models/*_strategy.h).
// This is the single source of truth; the R accessors model_version() /
// model_id() (R/strategy_support.R) read it through here rather than
// duplicating the number.

// [[Rcpp::export]]
int strategy_scientific_version(std::string type) {
  if (type == "FF16")  return plant::FF16_Strategy::scientific_version;
  if (type == "K93")   return plant::K93_Strategy::scientific_version;
  if (type == "TF24")  return plant::TF24_Strategy::scientific_version;
  if (type == "TF24f") return plant::TF24f_Strategy::scientific_version;
  Rcpp::stop("Unknown type " + type);
}
