#include <plant.h>

// The census exports stand_gradient() and stand_census() call, and the
// diagnostics that say what regimes a run's leaf solves met.

// One body per question, and one export per question that takes the runner and
// its model's name, read on the R side off the runner's own class. A metric is
// added to a strategy's own census_metrics(); nothing here names one, and the
// codomain is that list's length.

// Defined below the helpers that call it, with external linkage: the ladder's
// split export (gradient_ladder.cpp) builds its answer through it too.
Rcpp::List census_gradient_to_r(const plant::census_gradient& g);

namespace {

// The metrics' names, each carrying its kind ("census" or "fitness") as the
// R name, so the R side splits the two without parsing the metric's name.
template <class Runner>
Rcpp::CharacterVector metric_names_of(Runner& runner) {
  const auto metrics = runner.stand_metrics();
  Rcpp::CharacterVector ret(metrics.size());
  Rcpp::CharacterVector kinds(metrics.size());
  for (size_t m = 0; m < metrics.size(); ++m) {
    ret[m] = metrics[m].name;
    kinds[m] = Runner::metric_kind_name(metrics[m]);
  }
  ret.names() = kinds;
  return ret;
}

template <class Runner>
std::vector<std::vector<double>> state_adjoint_of(Runner& runner) {
  return runner.census_state_and_trait_rows().state.to_vectors();
}

// The parameters no gradient exists for, each with the sentence saying why, so
// asking for one is refused in the model's own words. Without this a caller who
// asks for `eta` is told the name is unknown, which it is not.
// The ones no equation reads go in the attribute "unread": a trade-off may
// derive one of those and the chain rule loses nothing by its missing column.
template <class Runner>
Rcpp::CharacterVector undifferentiable_of() {
  const auto reasons = Runner::strategy_type::undifferentiable_reasons();
  Rcpp::CharacterVector ret(reasons.size());
  Rcpp::CharacterVector names(reasons.size());
  std::vector<std::string> unread;
  for (size_t i = 0; i < reasons.size(); ++i) {
    names[i] = reasons[i].name;
    ret[i] = reasons[i].why;
    if (!reasons[i].read) {
      unread.push_back(reasons[i].name);
    }
  }
  ret.names() = names;
  ret.attr("unread") = unread;
  return ret;
}

// One row per census metric, one column per trait in the model's trait-name
// order. The active scalar lives inside this call and only doubles leave it.
//
// Absent `which_metrics` means every metric, which is what a caller that does
// not know the census asks for. Naming a subset is what makes a single metric
// cost one, and a name is what crosses: a position means a different metric as
// soon as the model's list changes.
template <class Runner>
Rcpp::List trait_gradient_of(Runner& runner,
                             Rcpp::Nullable<Rcpp::CharacterVector> which_metrics) {
  std::vector<std::string> wanted;
  if (which_metrics.isNotNull()) {
    const Rcpp::CharacterVector v(which_metrics);
    for (int i = 0; i < v.size(); ++i) {
      wanted.push_back(Rcpp::as<std::string>(v[i]));
    }
  }
  return census_gradient_to_r(runner.census_trait_gradient({}, wanted));
}

// The same rule as the metrics above: a name is what crosses. The entries arrive
// paired from Runner::gradient_control(), so nothing here chooses an order.
template <class Runner>
Rcpp::NumericVector control_of(Runner& runner) {
  const auto entries = runner.gradient_control();
  Rcpp::NumericVector out(entries.size());
  Rcpp::CharacterVector names(entries.size());
  for (size_t i = 0; i < entries.size(); ++i) {
    out[i] = entries[i].second;
    names[i] = entries[i].first;
  }
  out.names() = names;
  return out;
}

// How often each clamp site fired, one row per species and one column per site
// named, the environment's sites on the first row. On the forward run, or with
// `differentiated`, where the sweep ran, which is the only path a clamp severs a
// row on.
template <class Runner>
Rcpp::NumericMatrix clamp_counts_of(Runner& runner, bool differentiated) {
  const std::vector<std::vector<size_t>> counts =
      differentiated ? runner.clamp_counts_differentiated()
                     : runner.clamp_counts();
  Rcpp::NumericMatrix ret(static_cast<int>(counts.size()),
                          plant::CLAMP_SITE_COUNT);
  for (size_t i = 0; i < counts.size(); ++i) {
    for (int k = 0; k < plant::CLAMP_SITE_COUNT; ++k) {
      ret(static_cast<int>(i), k) = static_cast<double>(counts[i][k]);
    }
  }
  Rcpp::CharacterVector sites(plant::CLAMP_SITE_COUNT);
  for (int k = 0; k < plant::CLAMP_SITE_COUNT; ++k) {
    sites[k] = plant::clamp_site_name(k);
  }
  ret.attr("dimnames") = Rcpp::List::create(R_NilValue, sites);
  return ret;
}

using TF24_runner = plant::DeterministicPatchRunner<plant::TF24_Model<double>,
                                                    plant::TF24_Environment<double>>;
using TF24f_runner = plant::DeterministicPatchRunner<plant::TF24f_Model<double>,
                                                     plant::TF24_Environment<double>>;
using FF16_runner = plant::DeterministicPatchRunner<plant::FF16_Model<double>,
                                                    plant::FF16_Environment<double>>;

template <class Runner>
struct runner_tag {
  using type = Runner;
};

// Calls f with the runner type `model` names, as a tag: the models that
// declare a census, and the one place a model is added to them. Any other is
// refused here by name.
template <class F>
SEXP with_census_model(const std::string& model, F&& f) {
  if (model == "TF24") return Rcpp::wrap(f(runner_tag<TF24_runner>{}));
  if (model == "TF24f") return Rcpp::wrap(f(runner_tag<TF24f_runner>{}));
  if (model == "FF16") return Rcpp::wrap(f(runner_tag<FF16_runner>{}));
  Rcpp::stop("No census is declared for the " + model + " model, so it has no "
             "stand gradient; FF16, TF24 and TF24f declare one");
}

// Calls f with the runner `x` holds, typed by `model`.
template <class F>
SEXP with_census_runner(SEXP x, const std::string& model, F&& f) {
  return with_census_model(model, [&](auto tag) {
    using Runner = typename decltype(tag)::type;
    plant::RcppR6::RcppR6<Runner> runner = Rcpp::as<plant::RcppR6::RcppR6<Runner>>(x);
    return f(*runner);
  });
}

// A model whose strategy classifies each leaf solve's operating point and
// records the profit curvature its gradient met. Read off
// the strategy rather than listed by name, so a model that declares them is
// served without an entry here.
template <class Runner>
concept reports_incidence = requires(typename Runner::strategy_type& s) {
  s.operating_point_counts;
  s.curvature_margin;
};

// Calls f with the runner `x` holds, for a model that reports incidence, and
// refuses any other census model by name.
template <class F>
SEXP with_incidence_runner(SEXP x, const std::string& model, F&& f) {
  return with_census_runner(x, model, [&](auto& runner) -> SEXP {
    using Runner = std::decay_t<decltype(runner)>;
    if constexpr (reports_incidence<Runner>) {
      return f(runner);
    } else {
      Rcpp::stop("The " + model + " model solves no leaf, so it has no "
                 "operating points or curvature margin to report; TF24 and "
                 "TF24f report them");
      return R_NilValue;
    }
  });
}

// A per-species tally, as R reads it: one numeric vector per species.
SEXP counts_to_r(const std::vector<std::vector<size_t>>& counts) {
  std::vector<std::vector<double>> ret;
  ret.reserve(counts.size());
  for (const std::vector<size_t>& row : counts) {
    ret.push_back(std::vector<double>(row.begin(), row.end()));
  }
  return Rcpp::wrap(ret);
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
      // A refusal not localised to one species is NA, never -1: R reads a
      // negative index as "every species but".
      refusal[m] = Rcpp::List::create(
          Rcpp::_["reason"] = why.reason,
          Rcpp::_["species"] = why.species < 0 ? NA_INTEGER : why.species);
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

// `x` is a deterministic runner and `model` its model, as the runner's class
// names it.

// The metrics' names, each carrying its kind as the R name.
// [[Rcpp::export]]
SEXP census_metric_names(SEXP x, std::string model) {
  return with_census_runner(x, model,
                            [](auto& runner) { return metric_names_of(runner); });
}

// [[Rcpp::export]]
SEXP census_values(SEXP x, std::string model) {
  return with_census_runner(x, model, [](auto& runner) { return runner.census(); });
}

// [[Rcpp::export]]
SEXP census_state_adjoint(SEXP x, std::string model) {
  return with_census_runner(x, model,
                            [](auto& runner) { return state_adjoint_of(runner); });
}

// The gradient's columns: each species' differentiable parameters in
// ad_parameters() order, species-major, which is the order the trait adjoints
// accumulate in. Each name carries its species index, as "1.lma".
// [[Rcpp::export]]
SEXP census_trait_names(SEXP x, std::string model) {
  return with_census_runner(x, model, [](auto& runner) {
    return runner.r_patch().trait_adjoint_names();
  });
}

// [[Rcpp::export]]
SEXP census_undifferentiable(std::string model) {
  return with_census_model(model, [](auto tag) {
    return undifferentiable_of<typename decltype(tag)::type>();
  });
}

// [[Rcpp::export]]
SEXP census_trait_gradient(SEXP x, std::string model,
                           Rcpp::Nullable<Rcpp::CharacterVector> which_metrics = R_NilValue) {
  return with_census_runner(x, model, [&](auto& runner) {
    return trait_gradient_of(runner, which_metrics);
  });
}

// [[Rcpp::export]]
SEXP census_gradient_control(SEXP x, std::string model) {
  return with_census_runner(x, model,
                            [](auto& runner) { return control_of(runner); });
}

// The run's incidence, for every model whose leaf solve reports it: the
// regimes its solves met. The clamp counts have readers of their own.

// The forward run's classification tally: one row per species, one column per
// operating-point kind in census_operating_point_names() order. The leaf
// classifies by the branch taken and the next plant overwrites it, so a run's
// incidence is not recoverable afterwards from anything but this. It is what
// says whether a regime the gradient refuses is rare or is most of the run.
// [[Rcpp::export]]
SEXP census_operating_point_counts(SEXP x, std::string model) {
  return with_incidence_runner(x, model, [](auto& runner) {
    return counts_to_r(runner.operating_point_counts());
  });
}

// The kinds, in the order the counts are reported.
// [[Rcpp::export]]
std::vector<std::string> census_operating_point_names() {
  std::vector<std::string> ret;
  ret.reserve(plant::Leaf::operating_point_kind_count);
  for (size_t k = 0; k < plant::Leaf::operating_point_kind_count; ++k) {
    ret.push_back(plant::Leaf::operating_point_kind_name(
        static_cast<plant::Leaf::OperatingPointKind>(k)));
  }
  return ret;
}

// The smallest profit curvature the differentiated path met, one per species,
// or -1 where it met none. The guard on it refuses a row rather than returning
// amplification, and a guard that held reports the same green as a guard
// nothing reached -- so the distance to the floor is reported rather than
// assumed. gradient_control publishes the floor and stand_gradient_compare
// refuses two gradients taken at different ones, and this is the only thing
// that says how close a run came to it.
// [[Rcpp::export]]
SEXP census_curvature_margin(SEXP x, std::string model) {
  return with_incidence_runner(x, model, [](auto& runner) {
    return Rcpp::wrap(runner.curvature_margins());
  });
}

// [[Rcpp::export]]
SEXP census_clamp_counts(SEXP x, std::string model, bool differentiated = false) {
  return with_census_runner(x, model, [&](auto& runner) {
    return clamp_counts_of(runner, differentiated);
  });
}

// The clamp sites, in the order the counts' columns are reported. Read from
// the enum so a site cannot be counted under its neighbour's name.
// [[Rcpp::export]]
std::vector<std::string> census_clamp_names() {
  std::vector<std::string> ret;
  ret.reserve(plant::CLAMP_SITE_COUNT);
  for (int i = 0; i < plant::CLAMP_SITE_COUNT; ++i) {
    ret.push_back(plant::clamp_site_name(i));
  }
  return ret;
}

// Every tally the run reports, back to zero: a count carried in from an
// earlier run of the same runner reads as this run's.
// [[Rcpp::export]]
void census_clear_diagnostics(SEXP x, std::string model) {
  with_census_runner(x, model, [](auto& runner) {
    runner.clear_diagnostics();
    return R_NilValue;
  });
}
