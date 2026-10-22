#include <plant.h>

// The census exports stand_gradient() and stand_census() call.

// One body per question, instantiated per model below. A metric is added to a
// strategy's own census_metrics(); nothing here names one, and the codomain is
// that list's length. The R side picks the model's entry points off the SCM's
// own class, so a model declaring no census is refused there by name.

// Defined below the helpers that call it, with external linkage: the ladder's
// split export (gradient_ladder.cpp) builds its answer through it too.
Rcpp::List census_gradient_to_r(const plant::census_gradient& g);

namespace {

// The metrics' names, each carrying its kind ("census" or "fitness") as the
// R name, so the R side splits the two without parsing the metric's name.
template <class SCM>
Rcpp::CharacterVector metric_names_of(SCM& scm) {
  const auto metrics = scm.stand_metrics();
  Rcpp::CharacterVector ret(metrics.size());
  Rcpp::CharacterVector kinds(metrics.size());
  for (size_t m = 0; m < metrics.size(); ++m) {
    ret[m] = metrics[m].name;
    kinds[m] = SCM::metric_kind_name(metrics[m]);
  }
  ret.names() = kinds;
  return ret;
}

template <class SCM>
std::vector<std::vector<double>> state_adjoint_of(SCM& scm) {
  return scm.census_state_and_trait_rows().state.to_rows();
}

// The parameters no gradient exists for, each with the sentence saying why, so
// asking for one is refused in the model's own words. Without this a caller who
// asks for `eta` is told the name is unknown, which it is not.
template <class SCM>
Rcpp::CharacterVector undifferentiable_of() {
  const auto reasons = SCM::strategy_type::undifferentiable_reasons();
  Rcpp::CharacterVector ret(reasons.size());
  Rcpp::CharacterVector names(reasons.size());
  for (size_t i = 0; i < reasons.size(); ++i) {
    names[i] = reasons[i].first;
    ret[i] = reasons[i].second;
  }
  ret.names() = names;
  return ret;
}

// One row per census metric, one column per trait in the model's trait-name
// order. The active scalar lives inside this call and only doubles leave it.
//
// Absent `which_metrics` means every metric, which is what a caller that does
// not know the census asks for. Naming a subset is what makes a single metric
// cost one, and a name is what crosses: a position means a different metric as
// soon as the model's list changes.
template <class SCM>
Rcpp::List trait_gradient_of(SCM& scm,
                             Rcpp::Nullable<Rcpp::CharacterVector> which_metrics) {
  std::vector<std::string> wanted;
  if (which_metrics.isNotNull()) {
    const Rcpp::CharacterVector v(which_metrics);
    for (int i = 0; i < v.size(); ++i) {
      wanted.push_back(Rcpp::as<std::string>(v[i]));
    }
  }
  return census_gradient_to_r(scm.census_trait_gradient({}, wanted));
}

// The same rule as the metrics above: a name is what crosses. The entries arrive
// paired from SCM::gradient_control(), so nothing here chooses an order.
template <class SCM>
Rcpp::NumericVector control_of(SCM& scm) {
  const auto entries = scm.gradient_control();
  Rcpp::NumericVector out(entries.size());
  Rcpp::CharacterVector names(entries.size());
  for (size_t i = 0; i < entries.size(); ++i) {
    out[i] = entries[i].second;
    names[i] = entries[i].first;
  }
  out.names() = names;
  return out;
}

}  // namespace

// The gradient and, per metric, whether it was refused. External linkage:
// gradient_ladder.cpp's split export builds its answer through this too. Refusal is metric-level
// -- a sum has no defined value with an undefined term -- so one description
// serves the whole row, and a metric that answered carries none rather than an
// empty one, so the two are not the same object in R.
Rcpp::List census_gradient_to_r(const plant::census_gradient& g) {
  const size_t n_metric = g.gradient.size();
  Rcpp::List gradient(n_metric), refusal(n_metric);
  for (size_t m = 0; m < n_metric; ++m) {
    gradient[m] = Rcpp::wrap(g.gradient[m]);
    const plant::refusal& why = g.why[m];
    if (why.happened()) {
      refusal[m] = Rcpp::List::create(Rcpp::_["reason"] = why.reason,
                                      Rcpp::_["species"] = why.species);
    } else {
      refusal[m] = R_NilValue;
    }
  }
  return Rcpp::List::create(Rcpp::_["gradient"] = gradient,
                            Rcpp::_["refusal"] = refusal,
                            Rcpp::_["segments"] =
                                static_cast<double>(g.segments),
                            Rcpp::_["at_first_state"] =
                                Rcpp::wrap(g.at_first_state));
}

// ---- TF24 ----------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::CharacterVector census_metric_names_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return metric_names_of(*obj_);
}

// [[Rcpp::export]]
std::vector<double> census_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census();
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
census_state_adjoint_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return state_adjoint_of(*obj_);
}

// The gradient's columns: each species' differentiable parameters in
// ad_parameters() order, species-major, which is the order the trait adjoints
// accumulate in. Each name carries its species index, as "1.lma".
// [[Rcpp::export]]
std::vector<std::string>
census_trait_names_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->r_patch().trait_adjoint_names();
}

// [[Rcpp::export]]
Rcpp::CharacterVector census_undifferentiable_tf24() {
  return undifferentiable_of<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > >();
}

// [[Rcpp::export]]
Rcpp::List
census_trait_gradient_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                           Rcpp::Nullable<Rcpp::CharacterVector> which_metrics = R_NilValue) {
  return trait_gradient_of(*obj_, which_metrics);
}

// [[Rcpp::export]]
Rcpp::NumericVector gradient_control_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return control_of(*obj_);
}

// ---- FF16 ----------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::CharacterVector census_metric_names_ff16(plant::RcppR6::RcppR6<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > > obj_) {
  return metric_names_of(*obj_);
}

// [[Rcpp::export]]
std::vector<double> census_ff16(plant::RcppR6::RcppR6<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > > obj_) {
  return obj_->census();
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
census_state_adjoint_ff16(plant::RcppR6::RcppR6<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > > obj_) {
  return state_adjoint_of(*obj_);
}

// [[Rcpp::export]]
std::vector<std::string>
census_trait_names_ff16(plant::RcppR6::RcppR6<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > > obj_) {
  return obj_->r_patch().trait_adjoint_names();
}

// [[Rcpp::export]]
Rcpp::CharacterVector census_undifferentiable_ff16() {
  return undifferentiable_of<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > >();
}

// [[Rcpp::export]]
Rcpp::List
census_trait_gradient_ff16(plant::RcppR6::RcppR6<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > > obj_,
                           Rcpp::Nullable<Rcpp::CharacterVector> which_metrics = R_NilValue) {
  return trait_gradient_of(*obj_, which_metrics);
}

// [[Rcpp::export]]
Rcpp::NumericVector gradient_control_ff16(plant::RcppR6::RcppR6<plant::SCM<plant::FF16_Strategy<double>, plant::FF16_Environment<double> > > obj_) {
  return control_of(*obj_);
}
