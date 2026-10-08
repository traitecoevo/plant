#include <plant.h>

// The census exports stand_gradient() and stand_census() call.

// The census metrics and their trait gradient, for the one strategy that
// carries an active scalar. A metric is added to the strategy's own
// census_metrics(); nothing here names one, and the codomain is that list's
// length.

// [[Rcpp::export]]
std::vector<std::string> census_metric_names_tf24() {
  return plant::census_metric_names<plant::TF24_Strategy<double>>();
}

// [[Rcpp::export]]
std::vector<double> census_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census();
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
census_state_adjoint_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_state_and_trait_rows().state.to_rows();
}

// The gradient's columns: each species' differentiable parameters in
// ad_parameters() order, species-major, which is the order the trait adjoints
// accumulate in. Each name carries its species index, as "1.lma".
// [[Rcpp::export]]
std::vector<std::string>
census_trait_names_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->r_patch().trait_adjoint_names();
}

// The parameters no gradient exists for, each with the sentence saying why, so
// asking for one is refused in the model's own words. Without this a caller who
// asks for `eta` is told the name is unknown, which it is not.
// [[Rcpp::export]]
Rcpp::CharacterVector census_undifferentiable_tf24() {
  const auto reasons =
      plant::TF24_Strategy<double>::undifferentiable_reasons();
  Rcpp::CharacterVector ret(reasons.size());
  Rcpp::CharacterVector names(reasons.size());
  for (size_t i = 0; i < reasons.size(); ++i) {
    names[i] = reasons[i].first;
    ret[i] = reasons[i].second;
  }
  ret.names() = names;
  return ret;
}

// The gradient and, per metric, whether it was refused. Refusal is metric-level
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
                            Rcpp::_["ranges"] =
                                static_cast<double>(g.ranges),
                            Rcpp::_["at_first_state"] =
                                Rcpp::wrap(g.at_first_state));
}

// One row per census metric, one column per trait in census_trait_names_tf24()
// order. The active scalar lives inside this call and only doubles leave it.
// [[Rcpp::export]]
Rcpp::List
census_trait_gradient_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                           Rcpp::Nullable<Rcpp::CharacterVector> which_metrics = R_NilValue) {
  // Absent means every metric, which is what a caller that does not know the
  // census asks for. Naming a subset is what makes a single metric cost one, and
  // a name is what crosses: a position means a different metric as soon as the
  // model's list changes.
  std::vector<std::string> wanted;
  if (which_metrics.isNotNull()) {
    const Rcpp::CharacterVector v(which_metrics);
    for (int i = 0; i < v.size(); ++i) {
      wanted.push_back(Rcpp::as<std::string>(v[i]));
    }
  }
  return census_gradient_to_r(obj_->census_trait_gradient({}, wanted));
}

// The same rule as the metrics above: a name is what crosses. The entries arrive
// paired from SCM::gradient_control(), so nothing here chooses an order.
// [[Rcpp::export]]
Rcpp::NumericVector gradient_control_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const auto entries = obj_->gradient_control();
  Rcpp::NumericVector out(entries.size());
  Rcpp::CharacterVector names(entries.size());
  for (size_t i = 0; i < entries.size(); ++i) {
    out[i] = entries[i].second;
    names[i] = entries[i].first;
  }
  out.names() = names;
  return out;
}
