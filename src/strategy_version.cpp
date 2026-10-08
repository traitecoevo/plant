#include <plant.h>
#include <string>

// Scientific version of a model, read from the `scientific_version` constant
// declared on each strategy class (see inst/include/plant/models/*_strategy.h).
// Returned as a string so it can carry a compound version: TF24f is a fast
// approximation of TF24, so its version is "<TF24 version>.<revision>" (e.g.
// "2.1") -- the major component auto-tracks TF24 (a TF24 science change
// invalidates TF24f too), the minor is TF24f's own approximation revision.
// This is the single source of truth; the R accessors model_version() /
// model_id() (R/strategy_support.R) read it through here.

// [[Rcpp::export]]
std::string strategy_scientific_version(std::string type) {
  if (type == "FF16")
    return std::to_string(plant::FF16_Strategy::scientific_version);
  if (type == "K93")
    return std::to_string(plant::K93_Strategy::scientific_version);
  if (type == "TF24")
    return std::to_string(plant::TF24_Strategy<double>::scientific_version);
  if (type == "TF24f")
    return std::to_string(plant::TF24_Strategy<double>::scientific_version) + "." +
           std::to_string(plant::TF24f_Strategy<double>::approximation_revision);
  Rcpp::stop("Unknown type " + type);
}
