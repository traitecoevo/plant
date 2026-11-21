#include <plant.h>
#include <string>

// R-facing exports with one entry per model. A new model adds its entries to
// each section below, and its allometry export to src/strategy_expand.cpp.

// --- Node schedules ----------------------------------------------------------

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Model,plant::K93_Environment>& p) {
   return plant::schedule_default<plant::Parameters<plant::K93_Model,plant::K93_Environment> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Model,plant::K93_Environment>& p) {
   return plant::make_schedule<plant::Parameters<plant::K93_Model,plant::K93_Environment> >(p);
}

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___TF24__TF24_Env(const plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___TF24__TF24_Env(const plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___TF24f__TF24_Env(const plant::Parameters<plant::TF24f_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::TF24f_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___TF24f__TF24_Env(const plant::Parameters<plant::TF24f_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::TF24f_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___TF24t__TF24_Env(const plant::Parameters<plant::TF24t_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::TF24t_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___TF24t__TF24_Env(const plant::Parameters<plant::TF24t_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::TF24t_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___TF24w__TF24_Env(const plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___TF24w__TF24_Env(const plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>> >(p);
}

// --- Scientific versions -----------------------------------------------------

// Scientific version of a model, read from the `scientific_version` constant
// declared on each model class (see inst/include/plant/models/*_model.h).
// Returned as a string so it can carry a compound version: TF24f is an
// approximation of TF24, so its version is "<TF24 version>.<revision>" (e.g.
// "2.1") -- the major component auto-tracks TF24 (a TF24 science change
// invalidates TF24f too), the minor is TF24f's own approximation revision.
// TF24t is compound the same way: "<TF24 version>.<thermal revision>".
// This is the single source of truth; the R accessors model_version() /
// model_id() (R/strategy_support.R) read it through here.

// [[Rcpp::export]]
std::string strategy_scientific_version(std::string type) {
  if (type == "FF16")
    return std::to_string(plant::FF16_Model<double>::scientific_version);
  if (type == "K93")
    return std::to_string(plant::K93_Model::scientific_version);
  if (type == "TF24")
    return std::to_string(plant::TF24_Model<double>::scientific_version);
  if (type == "TF24f")
    return std::to_string(plant::TF24_Model<double>::scientific_version) + "." +
           std::to_string(plant::TF24f_Model<double>::approximation_revision);
  if (type == "TF24t")
    return std::to_string(plant::TF24_Model<double>::scientific_version) + "." +
           std::to_string(plant::TF24t_Model<double>::thermal_revision);
  if (type == "TF24w")
    return std::to_string(plant::TF24_Model<double>::scientific_version) + "." +
           std::to_string(plant::TF24w_Model<double>::private_water_revision);
  Rcpp::stop("Unknown type " + type);
}
