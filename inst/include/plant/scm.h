// -*-c++-*-
#ifndef PLANT_PLANT_SCM_H_
#define PLANT_PLANT_SCM_H_

#include <plant/node_schedule.h>
#include <plant/census_gradient.h>
#include <odelia/ode_solver.hpp>
#include <plant/patch.h>
#include <plant/scm_utils.h>

#include <odelia/sweep.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <tuple>

using namespace Rcpp;

namespace plant {

// A census differentiated on both of its input sets, one row per metric each.
// Kept as a pair because they come out of one recording: handing them back
// separately is what let two callers take two recordings of one function.
struct census_rows {
  odelia::ode::adjoint_rows state;  // one column per ODE state entry
  odelia::ode::adjoint_rows trait;  // one column per registered trait
};

// The refusal any of a patch's species recorded, with the one thing a species
// knows about itself filled in. Read after a recording rather than during one:
// what a leaf could not record is not something a void interface can return, so
// it is left where the leaf is and collected here.
template <typename Patch>
refusal recorded_refusal(Patch& patch) {
  for (size_t i = 0; i < patch.size(); ++i) {
    refusal why = *patch.at_species(i).strategy_ptr()->recorded_refusal;
    if (why.happened()) {
      why.species = static_cast<int>(i + 1);
      return why;
    }
  }
  return refusal{};
}

// It latches for as long as the storage behind it lives, so it is cleared where
// the call that reads it starts rather than trusted to be clean.
template <typename Patch>
void clear_recorded_refusal(Patch& patch) {
  for (size_t i = 0; i < patch.size(); ++i) {
    *patch.at_species(i).strategy_ptr()->recorded_refusal = refusal{};
  }
}

// SCM: the "Solver for Characteristics Method" driver.
//
// Owns a Patch (the population being integrated), a NodeSchedule (when each
// species' nodes are introduced), and an ODE Solver. It steps the patch
// forward by repeatedly: introducing all nodes due at the current time, then
// integrating the patch ODE system up to the next introduction time.
//
// The patch owns all of the ecology (fitness, offspring, competition, error
// computations); the SCM is the time-stepping/scheduling layer on top of it.
// Most r_* members are thin facades that expose the C++ API to R via RcppR6.
template <typename T, typename E> class SCM {
public:
  // ---- Type aliases ------------------------------------------------------
  typedef T                strategy_type;
  typedef E                environment_type;
  typedef Individual<T, E> individual_type;
  typedef Node<T, E>       node_type;
  typedef Species<T, E>    species_type;
  typedef Patch<T, E>      patch_type;
  typedef Parameters<T, E> parameters_type;

  // ---- Construction ------------------------------------------------------
  // An empty `ev` means "no events supplied": the schedule then comes from
  // p.node_schedule_times exactly as it always has (#522).
  SCM(parameters_type p, environment_type e, plant::Events ev, plant::Control c);

  // ---- Simulation lifecycle ----------------------------------------------

  // Run the whole schedule from t = 0 to completion.
  void run();

  // Integrate `p`'s strategies as invaders against the field the resident run
  // stood in: every one of them, the resident's own included, sees a field it
  // does not move. Destructive -- `p` becomes this SCM's parameters, and its
  // outputs are `p`'s. Needs a finished resident run to take a program from.
  void run_mutant(parameters_type p);

  // Run, keeping the state at each accepted step, and return one record per step.
  //
  // A run pinned to this run's times and sizes DOES reproduce these states, bit for
  // bit over 3381 steps, and pays 15% less because it attempts no step it will
  // reject. This said the opposite -- that a rejected attempt moves patch state
  // which is not ODE state, and a pinned run makes none -- and the states are
  // measurably identical, so whatever a rejected attempt leaves behind is rebuilt
  // from the state before anything reads it.
  // The record the run kept, read in place off the solver: one row per accepted
  // step carrying its time, the size that reached it and the state there. Not a
  // copy and not a re-pairing -- the solver holds exactly this.
  using trajectory = std::span<const odelia::ode::step_record<patch_type>>;
  trajectory store_trajectory();

  // Set before run() to keep the state at every accepted step. The reverse pass
  // needs those states and cannot recover them from a finished run, so a run that
  // did not keep them has to be repeated -- one whole forward integration. Off by
  // default because the store is one double per state entry per step and a forward
  // run has no use for it.
  bool record_trajectory = false;

  // Adaptively refine the node-introduction schedule entirely in C++:
  // repeatedly run, flag nodes whose combined error exceeds schedule_eps,
  // and bisect the interval below each flagged node (upwind scheme), up to
  // schedule_nsteps times. Replaces the R build_schedule loop.
  void refine_schedule();

  // Return patch, schedule and solver to their t = 0 state; clear history.
  void reset();

  // The resident run an invasion stands in: one row per accepted step, carrying
  // the field each of its stages was taken in. Kept whole rather than re-derived,
  // because the point of an invasion sweep is many invaders against ONE resident.
  // Its program is not held beside it: a recording IS a program, row by row, and
  // two containers that have to agree can be made to disagree.
  std::vector<odelia::ode::step_record<patch_type>> resident_recording;

  // True once every scheduled node introduction has been consumed.
  bool complete() const;

  // Current patch time.
  double time() const;

  // ---- Outputs -----------------------------------------------------------
  // Total (not per-capita) offspring. These delegate to the patch, which owns
  // the fitness/offspring computations.
  std::vector<double> net_reproduction_ratios() const { return patch.net_reproduction_ratios(); }
  std::vector<double> offspring_production() const { return patch.offspring_production(); }

  // Every metric the strategy declares, summed over the species, in the order it
  // declares them. The codomain is the list's length.
  std::vector<double> census() const;

  // What a census of this model reads, at whatever scalar `P` carries. The one
  // place a strategy is asked, and the one place the question is refused for a
  // strategy that does not answer it: a census over a list nothing declared
  // would otherwise fail inside whichever loop reached for it.
  template <class P>
  static const auto& metrics_of() {
    static_assert(Censusable<typename P::strategy_type>,
                  "a census of this model is being taken, so its strategy must "
                  "declare census_metrics()");
    return P::strategy_type::census_metrics();
  }

  // One metric summed over every species of `p`. Templated on the patch so the
  // value and its derivative are the same reduction at two scalars.
  template <class P>
  static typename P::value_type census_sum(
      const P& p, const census_metric<typename P::strategy_type>& metric) {
    typename P::value_type tot = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
      tot += p.at_species(i).census_integral(metric);
    }
    return tot;
  }

  // The reverse pass runs on the birth-date coordinate only, and refuses the
  // other one here rather than answering it. On the height coordinate the
  // abscissa is state, so the quadrature weights carry a derivative nothing
  // supplies and the density rate carries a compression term the recorded step
  // does not compute: the sweep is then the transpose of a function the forward
  // model is not evaluating. Nothing about the arithmetic
  // complains, and the two coordinates are different functions rather than two
  // discretisations of one -- one census metric's trait sensitivity changes
  // sign between them -- so the answer would be finite, plausible and wrong.
  void require_birth_date_coordinate(const char* entry) const {
    if (!control.node_density_in_birth_date) {
      util::stop(std::string(entry) +
                 ": the reverse-mode gradient runs on the birth-date "
                 "size-density coordinate only. Set "
                 "control$node_density_in_birth_date = TRUE and re-run.");
    }
  }

  // d(census)/d(ODE state) and d(census)/d(trait) at the current time, one row
  // per metric each, from one recording. The state half is what the reverse pass
  // is seeded with; the trait half is what no sweep produces, because a metric
  // reads the traits itself and the boundary node's own quantities are rebuilt
  // when the state is set. Columns as ode_state writes them, and as
  // census_trait_gradient reports them.
  census_rows census_state_and_trait_rows() const;

  // d(census)/d(trait), one row per metric and one column per trait in each
  // strategy's ad_parameters() order, species-major. Requires an adaptive run to
  // have resolved the schedule this replays.
  //
  // `extra_stops` names recorded steps at which the sweep stops and resumes. The
  // adjoint recursion is linear in the step, so composition over steps is
  // associative and any split must give the same numbers bit for bit; a
  // difference is something carried across a step boundary that is not the
  // adjoint. Splits outside a range's interior are ignored, so a caller may
  // pass a boundary index without special-casing it.
  //
  // Returns the numbers AND what each of them is. The two travel together
  // because a row of doubles cannot say whether it is an answer: a refusal
  // anywhere in a metric's sweep makes that metric's whole gradient undefined,
  // and not-a-number is the only mark a bare matrix carries. Within a metric
  // that answered, every column is a number the sweep computed -- an exact zero
  // included, which is why a parameter no gradient exists for is refused by name
  // rather than given a column.
  census_gradient
  // `which_metrics` names the rows to sweep, empty meaning every one. A metric
  // not asked for is not seeded and not swept, so asking for one costs one --
  // which is what a caller differentiating a single census wants and what
  // computing all of them and subsetting the answer does not give.
  census_trait_gradient(const std::vector<size_t>& extra_stops = {},
                        const std::vector<std::string>& which_metrics = {});

  // The four references the ladder checks census_trait_gradient against.

  // The same quantity differenced, by moving the prepared strategy exactly where
  // the recording seeds it. It referees the trait half above while sharing none of
  // it: that one records the census and sweeps a tape, this one evaluates the
  // census twice. A difference that rebuilt from Parameters would re-run
  // preparation and carry the birth-size channel the differentiated path imposes
  // to zero, so this one perturbs in place.
  std::vector<std::vector<double>> census_trait_difference(double rel);

  // One exact directional derivative of the census, by a forward tangent of the
  // same trajectory stepped at the sizes the run recorded. `direction` carries
  // one weight per trait, species-major in each strategy's ad_parameters()
  // order; a coordinate direction gives one Jacobian column and a mixed one a
  // contraction. Returns one tangent per metric, and writes the metrics the
  // replay itself reached: a reference whose value disagrees with the model is a
  // reference to a different function, and the gap is this check's own floor.
  std::vector<double> census_trait_tangent(const std::vector<double>& direction,
                                           std::vector<double>& value);

  // One exact directional derivative of the census with respect to the first
  // recorded state, by a forward tangent stepped at the sizes the run recorded.
  // `direction` carries one weight per component of that state.
  //
  // No trait, no derived quantity and no census direct term is on this path, so
  // it isolates how the trajectory carries a perturbation to a state a cohort
  // starts at, and census_initial_state_replay is what referees it.
  //
  // `range` picks where the seeding happens, and only a state a sweep can be
  // RE-ENTERED at will do -- one whose width matches a piece's start, since a
  // walk resuming anywhere else is a width apart from the rows it would carry.
  // `j >= 1` is the state reached just after the jth introduction, which no
  // record holds and the map rebuilds. 0 is the first recorded state, and what
  // it carries is the run's setup: nothing, an introduction the schedule placed
  // at the initial time, or the cohorts a resumed patch was seeded with, which
  // can be several per species at birth dates before the run began.
  std::vector<double>
  census_initial_state_tangent(const std::vector<double>& direction,
                               std::vector<double>& value,
                               size_t range = 0);

  // The census a plain-double replay of the recorded steps reaches from
  // `state0`. Differencing it moves the state the tangent above seeds, through
  // the same steps and the same introductions, so the two differentiate one
  // function and a disagreement is the propagation's own.
  std::vector<double>
  census_initial_state_replay(const std::vector<double>& state0,
                              size_t range = 0);

  // The replay both entry points above run. `seed` fills the scalar state the
  // replay starts from, given the recorded one.
  //
  // Storing a trajectory runs the model, so the seeding is handed in rather than
  // applied by the caller: a caller reading the recorded state for itself would
  // store twice and run twice.
  template <class Scalar, class Seed>
  std::vector<Scalar> replay_initial_state(size_t from_range, Seed seed);

  // The Control entries that move the trajectory or move which states answer,
  // and so move the gradient. The curvature floor is here for the second reason
  // rather than the first: it changes no forward number and still decides which
  // rows exist.
  //
  // ⚠️ EACH NAME BESIDE ITS VALUE, because two of the five are 1e-3 at the
  // defaults. The names used to be attached positionally in R, one file away from
  // the order built here -- so transposing ci_abs_tol and gradient_curvature_floor
  // left both readable, both plausible, and stand_gradient_compare() refusing on
  // the wrong pair. census_gradient.cpp states the rule this now follows: a name
  // is what crosses, because a position means a different entry as soon as the
  // list changes.
  std::vector<std::pair<std::string, double>> gradient_control() const {
    return {{"GSS_tol_abs", control.GSS_tol_abs},
            {"ci_abs_tol", control.ci_abs_tol},
            {"node_gradient_eps", control.node_gradient_eps},
            {"schedule_eps", control.schedule_eps},
            {"gradient_curvature_floor", control.gradient_curvature_floor}};
  }

  // ---- R interface -------------------------------------------------------

  // Run / parameters / state access
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }

  // Every diagnostic below reads the LIVE system rather than r_patch(), and the
  // reason is narrower than r_patch() being a snapshot: it is live for the tallies
  // a strategy owns and stale for the ones the environment owns.
  //
  // Species::strategy_ptr() hands back the shared_ptr by value, so a copied Patch
  // shares one strategy with the live one and every tally on it aliases -- which
  // is why the ladder reads leaf_placements() straight off r_patch(). The
  // environment is a Patch member by value, so its half is a real copy, and a
  // sweep advances the live one afterwards through be_at_step's
  // compute_environment. Reading the live system is what covers both halves under
  // one rule.

  // The classification tally the FORWARD run built, one row per species and one
  // column per operating-point kind.
  std::vector<std::vector<size_t>> operating_point_counts() {
    const patch_type& live = solver.get_system_ref();
    std::vector<std::vector<size_t>> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      ret.push_back(live.at_species(i).strategy_ptr()->operating_point_counts);
    }
    return ret;
  }
  // The same for the clamp sites, which are counted for the same reason: a
  // clamped row and a true zero are the same number.
  //
  // The sites live in two objects and the list is one, so the environment's
  // tally is ADDED to the species' rather than reported beside it -- a caller
  // asking how often a site fired is asking about the site, not about which
  // object happened to reach it. The environment is one, so its counts land on
  // the first row.
  std::vector<std::vector<size_t>> clamp_counts() {
    const patch_type& live = solver.get_system_ref();
    std::vector<std::vector<size_t>> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      const auto s = live.at_species(i).strategy_ptr();
      std::vector<size_t> row = s->clamps.forward;
      // The leaf's own sites keep ONE tally across both paths, because the leaf
      // solves in double on both and every rebound copy shares its storage. So the
      // forward share is the total less what the sweep was measured to take.
      const std::vector<std::size_t> leaf_total = s->leaf.clamp_counts();
      const std::vector<size_t>& swept = *s->clamps.differentiated;
      for (std::size_t k = 0; k < leaf_total.size(); ++k) {
        const std::size_t at = CLAMP_LEAF_FIRST + k;
        if (at >= CLAMP_SITE_COUNT) { break; }
        row[at] += leaf_total[k] > swept[at] ? leaf_total[k] - swept[at] : 0;
      }
      ret.push_back(row);
    }
    add_environment_clamps(ret, live.environment_clamps().forward);
    return ret;
  }
  // And the same sites counted where the sweep runs, which is the only path a
  // clamp severs a row on. The forward tally cannot stand in for this one: the
  // sweep visits the recorded steps rather than every solve.
  std::vector<std::vector<size_t>> clamp_counts_differentiated() {
    const patch_type& live = solver.get_system_ref();
    std::vector<std::vector<size_t>> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      ret.push_back(*live.at_species(i).strategy_ptr()->clamps.differentiated);
    }
    add_environment_clamps(ret, *live.environment_clamps().differentiated);
    return ret;
  }
  // The smallest profit curvature the differentiated path met, one per species.
  // A guard that held and a guard nothing reached report the same green, so the
  // distance to the floor is carried out beside the refusals it did not raise.
  std::vector<double> curvature_margins() {
    const patch_type& live = solver.get_system_ref();
    std::vector<double> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      ret.push_back(*live.at_species(i).strategy_ptr()->curvature_margin);
    }
    return ret;
  }

  // Every diagnostic this SCM reports, back to where a run starts from: the
  // operating-point tally, the strategy's clamp counts on both paths, the leaf's
  // and its root model's, the curvature margin, and the environment's clamps.
  //
  // One call rather than five, because a caller wanting one of them before a run
  // wants all of them: a count carried in from an earlier run reads as this run's.
  void clear_diagnostics() {
    const patch_type& live = solver.get_system_ref();
    for (size_t i = 0; i < live.size(); ++i) {
      std::vector<size_t>& c =
        live.at_species(i).strategy_ptr()->operating_point_counts;
      c.assign(c.size(), 0);
      live.at_species(i).strategy_ptr()->clamps.clear();
      live.at_species(i).strategy_ptr()->leaf.clear_clamp_counts();
      *live.at_species(i).strategy_ptr()->curvature_margin = -1.0;
    }
    live.environment_clamps().clear();
  }
  const std::vector<patch_type> &r_history() const { return history; }

  // How many times the inflow boundary was evaluated over one sweep. It is
  // evaluated once per rate evaluation, and a step is recorded once and swept
  // per metric, so the count is six per step and does not scale with the metrics
  // asked for. A row that acts once per stage is multiplied by that count, so
  // the count is part of the row and belongs beside its value.
  size_t boundary_condition_evaluations() { return solver.recorded_rates(); }
  void clear_boundary_condition_evaluations() { solver.clear_recorded_rates(); }
  Rcpp::List r_get_state() const { return patch.r_get_state(); };

  // Fitness / reproduction
  double r_net_reproduction_ratio_for_species(util::index species_index) const;
  std::vector<std::vector<double>> r_net_reproduction_ratio_errors() const;

  // Schedule-refinement error signals.
  // Per-node refinement error: element-wise max of the competition error
  // (sampled during the run) and the reproduction error (computed at the end).
  // This is the signal that drives refine_schedule().
  std::vector<std::vector<double>> refinement_error_by_node() const;
  std::vector<double>
  r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const;

  // Node schedule access
  NodeSchedule r_node_schedule() const { return node_schedule; }
  // The schedule read back out in the R-facing wire format, so that a run
  // whose schedule was refined can be re-run, or inspected, as events.
  Events r_events() const {
    return events_from_schedule_events(node_schedule.get_events());
  }
  // What the events actually did, in the order they were applied. Cleared by
  // reset(), so it always describes the run you are looking at.
  EventLog r_event_log() const { return event_log_from_records(event_log); }
  void r_set_node_schedule(NodeSchedule x);
  void r_set_node_schedule_times(std::vector<std::vector<double>> x);

  // ODE times: the step times the solver actually used on the last run. Carried
  // back on the parameters (`p$ode_times`) they make the next run a replay --
  // holding them is what says so, and no flag says it separately.
  std::vector<double> r_ode_times() const;

  // The size of the step that reached each of r_ode_times(), NaN first. Carried
  // BESIDE the times the next run repeats this one exactly; with the times
  // alone it steps TO each of them and chooses its own sub-steps, which is what
  // a difference of two runs wants -- one grid, each side free to reach it.
  std::vector<double> r_ode_step_sizes() const;

  // The trajectory as a list of records, each a time, the step size that reached it,
  // and the state there.
  Rcpp::List r_store_trajectory();

  // ---- Public state ------------------------------------------------------
  // The two toggles are exposed to R directly (access: field), so they need
  // no getter/setter wrappers.
  bool collect;                    // record a patch snapshot after each step
  bool collect_refinement_errors;  // accumulate competition errors during run
  std::vector<patch_type> history; // per-step patch snapshots when collect

private:
  // One environment serves every species, so its tally has no species of its own;
  // it lands on the first row rather than being dropped or duplicated. Private
  // because the two readers above are the only callers.
  static void add_environment_clamps(std::vector<std::vector<size_t>>& ret,
                                     const std::vector<size_t>& env) {
    if (ret.empty()) {
      ret.push_back(env);
      return;
    }
    for (size_t s = 0; s < env.size() && s < ret[0].size(); ++s) {
      ret[0][s] += env[s];
    }
  }

  // Upwind bisection: insert the midpoint of the interval below each flagged
  // node. Mirrors split_times() in build_schedule.R.
  static std::vector<double> bisect_flagged_intervals(const std::vector<double>& times,
                                                      const std::vector<bool>& split);

  // Uniform grid for fixed-step forward-Euler integration (control.fixed_time_step).
  static std::vector<double> uniform_euler_times(double t0, double t1, double dt);

  // Advance one event: introduce every node due at the current time, then
  // integrate to the next introduction (or over the pinned ode times). Returns
  // the species introduced. The solver owns the patch system, so the live state
  // is solver.get_system_ref(); run() refreshes the `patch` snapshot once, after
  // its loop, rather than per event.
  std::vector<size_t> run_next();

  parameters_type parameters;
  Control control;
  patch_type patch;
  NodeSchedule node_schedule;
  // Lives on the runner rather than the patch: the patch is copied into
  // `history` once per step, and a log that grew with the run would be copied
  // with it every time.
  std::vector<EventRecord> event_log;
  odelia::ode::Solver<patch_type> solver;
};

// ---- Construction --------------------------------------------------------

template <typename T, typename E>
SCM<T, E>::SCM(parameters_type p, environment_type e, Events ev, Control c)
    : parameters(p), control(c), patch(parameters, e, c),
      node_schedule(make_node_schedule(parameters, ev)),
      solver(patch, make_ode_control(c)) {

  parameters.validate();

  collect = false;
  collect_refinement_errors = false;
  solver.set_collect(false);

  if (!util::identical(parameters.patch_area, 1.0)) {
    util::warning("We recommened keeping patch_area = 1 for the SCM, as need to check units for all other sizes");
  }

}

// Build a uniform grid {t0, t0 + dt, ..., t1} with spacing dt, starting exactly
// at t0 and ending exactly at t1 (the final interval may be shorter than dt).
// Used to drive forward-Euler integration between schedule events.
template <typename T, typename E>
std::vector<double> SCM<T, E>::uniform_euler_times(double t0, double t1,
                                                   double dt) {
  std::vector<double> times;
  times.push_back(t0);
  if (t1 <= t0) {
    return times;
  }
  // Number of (mostly dt-sized) intervals; the small tolerance avoids spawning
  // a spurious tiny final interval when (t1 - t0) is an FP-near multiple of dt.
  const size_t n =
      static_cast<size_t>(std::ceil((t1 - t0) / dt - 1e-10));
  for (size_t i = 1; i < n; ++i) {
    times.push_back(t0 + static_cast<double>(i) * dt);
  }
  times.push_back(t1); // exact endpoint
  return times;
}

// ---- Simulation lifecycle ------------------------------------------------

template <typename T, typename E> void SCM<T, E>::run() {
  // Set before reset(), which records the state the run starts from. One flag, and
  // it is the solver's: the choices this run's rate evaluations make are the same
  // recording as its states, so what keeps the states is what says the choices are
  // kept too.
  solver.set_keep_states(record_trajectory);
  reset();
  // The solver owns the live patch system; operate on it directly during the
  // run and avoid per-step copies into the `patch` member.
  if (collect) {
    history.push_back(solver.get_system_ref());
  }

  while (!complete()) {
    std::vector<size_t> added = run_next();
    if (collect_refinement_errors) {
      solver.get_system_ref().collect_competition_errors(added);
    }
    if (collect) {
      history.push_back(solver.get_system_ref());
    }
  }

  // Expose the final state through the `patch` accessor after the loop.
  patch = solver.get_system_ref();
}

template <typename T, typename E>
std::vector<size_t> SCM<T, E>::run_next() {
  std::vector<size_t> ret;
  const double t0 = time();
  // The live patch system is owned by the solver; mutate it in place.
  auto &sys = solver.get_system_ref();

  const schedule_entry& intro = node_schedule.next();

  // Resume support: if the next scheduled introduction is in the future,
  // integrate the gap up to it without introducing any node. This happens on
  // the first step of a run resumed from an exported state -- the patch is
  // already populated (in reset()) and starts at parameters.initial_time, which
  // falls before the first residual schedule entry. It never happens for an
  // empty patch, whose schedule always starts at t0 = 0, so the normal path
  // below is unchanged. The next call will then introduce at that time.
  if (intro.time > t0) {
    solver.set_state_from_system();
    if (node_schedule.using_ode_steps()) {
      util::stop("Resuming from an initial state is not supported for "
                 "replaying a recorded run");
    } else if (control.fixed_time_step > 0.0) {
      solver.advance_euler(
          uniform_euler_times(t0, intro.time, control.fixed_time_step));
    } else {
      solver.advance_adaptive({solver.time(), intro.time});
    }
    return ret; // empty: nothing introduced this step
  }

  if (!util::identical(t0, intro.time)) {
    util::stop("Start time not what was expected");
  }
  // The species this introduction names, which the schedule grouped when it was
  // set rather than the run regrouping them by walking equal times.
  ret = intro.species;
  const double t_end = node_schedule.time_end();
  node_schedule.pop();

  // Every action at this instant applies before the introductions, so a node
  // introduced here sees the post-event environment. The schedule holds
  // them already ordered by event_type_rank, so this is a walk and not a sort --
  // and `intro` stays valid across the pop, which only moved the cursor.
  for (const auto& a : intro.actions) {
    event_log.push_back(sys.apply_event(a));
  }
  sys.introduce_nodes(ret, intro.time);
  solver.set_state_from_system();
  // The insertion, as its own row: it holds the wider state the introduction just
  // reached, which is what the next step runs from and which no step reached.
  // Recorded here because this is where the width changes.
  solver.push_insertion();

  // Three integration modes:
  //  - pinned ode times (resident replay for a mutant): step exactly to the
  //    cached times via the full RKCK stepper, by their recorded step sizes
  //    when the schedule carries them;
  //  - fixed-step forward Euler (control.fixed_time_step > 0): walk a uniform
  //    sub-grid between this introduction and the next;
  //  - otherwise: adaptive, error-controlled RKCK to the next introduction.
  if (node_schedule.using_ode_steps()) {
    if (control.fixed_time_step > 0.0) {
      // A recorded field is kept per rate evaluation, and forward Euler makes
      // one per step where RKCK makes six -- so a program recorded under one
      // cannot be replayed under the other. Refuse rather than mis-integrate.
      util::stop("fixed_time_step (forward Euler) is not supported for a pinned "
                 "ODE schedule");
    }
    // Each recorded step carries the size it took and the time it reached, so a
    // replay lands where the run landed rather than a rounding short of it.
    // Stepping to the times instead would take different steps, because a size
    // differenced back out of two recorded times is not the size that was taken.
    //
    // A size is NaN where the schedule is a grid rather than a recording, and
    // the step is taken TO that time instead. The schedule holds no step at an
    // interval's own end -- a run stops there by clamping its last step to the
    // boundary, and the step below does the same arithmetic.
    // An interval with no step inside it is one the schedule crosses in a single
    // step, so there is nothing to replay before the step to its end. Deciding
    // that per interval rather than per schedule is what silently integrated
    // those intervals adaptively.
    const std::vector<odelia::ode::instruction> inside =
        node_schedule.program_within(t0, t_end);
    if (!inside.empty()) {
      solver.advance_recorded(inside);
    }
    solver.advance_fixed({solver.time(), t_end});
  } else if (control.fixed_time_step > 0.0) {
    solver.advance_euler(
        uniform_euler_times(t0, t_end, control.fixed_time_step));
  } else {
    solver.advance_adaptive({solver.time(), t_end});
  }

  return ret;
}

// An invader integrates against a field it does not move: every strategy in `p`
// is an invader, the resident's own among them, and that one coming back with the
// resident's fitness is what says the machinery is sound.
//
// Two passes, because the field lives at the Runge-Kutta stages and no recorded
// step boundary carries one. The first replays the resident over its own recorded
// program, keeping the field at each rate evaluation; the second walks that
// recording with `p`'s strategies in the patch, and each evaluation LOADS the
// field the resident stood in rather than building its own. Which of the two a
// pass is doing is the constness of what the walk hands over, so there is no mode
// here and none to go stale.
//
// ⚠️ THE RECORDING PASS IS PINNED AND NOT ADAPTIVE. An adaptive run evaluates
// inside attempts it then rejects; pinned to its own program it rejects nothing,
// so what it keeps is exactly the sequence a replay of that program makes.
//
// The replay needs no schedule of its own: the recording carries its own
// insertions, and the patch answers each one from ITS introduction times, so `p`'s
// species are introduced where the resident's were. That is the same map the sweep
// transposes rather than a second spelling of it.
//
// The record outlives the call, so a second invasion against the same resident
// pays for one pass and not two -- and, as on develop, it is still the FIRST
// resident's field, because this call has overwritten `parameters` with `p`.
template <typename T, typename E>
void SCM<T, E>::run_mutant(parameters_type p) {
  if (resident_recording.empty()) {
    const std::vector<double> times = r_ode_times();
    const std::vector<double> sizes = r_ode_step_sizes();
    // Two, not one: a run that never stepped still reports the instant it
    // started at, and a program of one entry is a start with nothing after it.
    if (times.size() < 2) {
      util::stop("Run a resident first to generate a competitive landscape");
    }
    // Pinned to its own program, so it rejects nothing and what it keeps is
    // exactly the sequence a replay of that program makes.
    node_schedule.r_set_ode_steps(times, sizes);
    node_schedule.reset();
    patch.set_keep_field(true);
    const bool kept = record_trajectory;
    record_trajectory = true;
    run();
    record_trajectory = kept;
    patch.set_keep_field(false);
    const trajectory rec = solver.recording();
    resident_recording.assign(rec.begin(), rec.end());
  }

  // Destructive, as it has always been: the mutants become this SCM's community,
  // and its outputs are theirs.
  parameters = p;
  patch.overwrite_strategies(parameters.strategies);
  node_schedule = make_node_schedule(parameters);
  reset();
  // ⚠️ A REPLAY KEEPS NOTHING. It is the recording's consumer, not another one of
  // them, and a recording of a replay is a copy of its own input that nobody
  // reads -- a state vector and six rows per step, for every invader in a sweep.
  solver.set_keep_states(false);
  solver.advance_recorded(resident_recording);
  patch = solver.get_system_ref();
}

template <typename T, typename E>
typename SCM<T, E>::trajectory SCM<T, E>::store_trajectory() {
  // Run only if this run kept no states to read. Reading does not consume them,
  // and a walk puts the patch back on the last recorded step when it is done, so
  // a second consumer reads the same record rather than repeating the run.
  //
  // Asked of this flag and not of the solver: run() is what sets the solver's
  // from this one, so asking the solver is asking it what it was just told.
  if (!record_trajectory) {
    record_trajectory = true;
    run();
  }

  // Handed back as the solver holds it. The state, the time it was reached at and
  // the size that reached it are one record there, so nothing here pairs them and
  // nothing copies them -- projecting the three apart and rebuilding the same row
  // was work whose only product was a chance to mispair.
  return solver.recording();
}

// Upwind bisection of flagged intervals. For each flagged node j (j >= 1; the
// first and last nodes are never flagged), insert the midpoint of the interval
// (t[j-1], t[j]). Equivalent to sort(c(times, times[i] - dt[i-1]/2)) in R.
template <typename T, typename E>
std::vector<double> SCM<T, E>::bisect_flagged_intervals(const std::vector<double>& times,
                                                        const std::vector<bool>& split) {
  std::vector<double> ret = times;
  for (size_t j = 1; j < split.size(); ++j) {
    if (split[j]) {
      ret.push_back(0.5 * (times[j] + times[j - 1]));
    }
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void SCM<T, E>::refine_schedule() {
  collect_refinement_errors = true;
  const double eps = control.schedule_eps;

  for (size_t step = 0; step < control.schedule_nsteps; ++step) {
    run(); // resets, then runs with collect_refinement_errors set

    std::vector<std::vector<double>> node_error = refinement_error_by_node();

    // Flag nodes whose refinement error exceeds the threshold.
    std::vector<std::vector<bool>> split(node_error.size());
    bool any = false;
    for (size_t i = 0; i < node_error.size(); ++i) {
      split[i].assign(node_error[i].size(), false);
      for (size_t j = 0; j < node_error[i].size(); ++j) {
        if (node_error[i][j] > eps) {
          split[i][j] = true;
          any = true;
        }
      }
    }
    if (!any) {
      break; // converged: no interval needs refining
    }

    // Bisect flagged intervals and install the denser schedule.
    std::vector<std::vector<double>> times = node_schedule.get_times();
    for (size_t i = 0; i < times.size(); ++i) {
      times[i] = bisect_flagged_intervals(times[i], split[i]);
    }
    node_schedule.set_times(times);
  }

  // Leave Parameters self-describing: record the refined schedule and the
  // ode times from the final run (mirrors build_schedule.R).
  parameters.node_schedule_times = node_schedule.get_times();
  // Carried out in the parameters so a later run of them replays this schedule,
  // which is what these two fields are for on the way IN. Read off the solver's
  // one record rather than through two accessors: a time and the size that
  // reached it are paired there, and pairing them again here is a chance to
  // pair them wrong.
  const std::vector<odelia::ode::instruction> taken = solver.schedule();
  parameters.ode_times.clear();
  parameters.ode_step_sizes.clear();
  parameters.ode_times.reserve(taken.size());
  parameters.ode_step_sizes.reserve(taken.size());
  for (const odelia::ode::instruction& step : taken) {
    parameters.ode_times.push_back(step.time);
    parameters.ode_step_sizes.push_back(step.step_size);
  }
}

// NOTE: solver.reset() sets the solver's internal time to zero. There is
// currently no other way to set that time; it might be cleaner to add an
// odelia::ode::Solver::set_time and call set_time(0) explicitly here.
template <typename T, typename E> void SCM<T, E>::reset() {
  // The schedule may have been changed since the patch was built, and a
  // reconciliation during the sweep reads it to work out the shape at a step.
  // Refreshed where the run that uses it begins, which is the one place both
  // the patch and the solver are put back to t = 0.
  patch.set_introduction_times(parameters.node_schedule_times);
  patch.reset();
  node_schedule.reset();
  // Seed the solver's owned system from the freshly reset patch, then reset
  // the solver's time/step state and sync the snapshot back.
  solver.get_system_ref() = patch;
  solver.reset();
  patch = solver.get_system_ref();
  history.clear();
  event_log.clear();
}


template <typename T, typename E> bool SCM<T, E>::complete() const {
  return node_schedule.remaining() == 0;
}

template <typename T, typename E> double SCM<T, E>::time() const {
  return solver.time();
}

// ---- R interface ---------------------------------------------------------
//
// The fitness/offspring and per-node error computations live on the patch
// (patch.h); the SCM methods below are thin facades that preserve the R API.
//
// Several of these are diagnostic/inspection hooks rather than part of the
// production run path: outside the C++ refinement loop they are only called
// from the test suite and the node_spacing vignette (noted per method below).

// Per-species fitness: the net reproduction ratio (expected offspring per seed)
// for one species. A genuine biological quantity, not just a diagnostic.
template <typename T, typename E>
double SCM<T, E>::r_net_reproduction_ratio_for_species(
    util::index species_index) const {
  const size_t idx = species_index.check_bounds(patch.size());
  auto scalars = std::vector<double>(patch.at_species(idx).size(), 1.0);
  return patch.net_reproduction_ratio_for_species(idx, scalars);
}

// Diagnostic: per-node discretisation error in the reproduction integral. One
// of the two components of refinement_error_by_node; exposed for inspection
// and validation (tests / vignette).
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::r_net_reproduction_ratio_errors() const {
  return patch.net_reproduction_ratio_errors();
}

// The combined per-node refinement error. Used internally by refine_schedule();
// also exposed to R so tests / the vignette can inspect the signal that drives
// schedule refinement.
template <typename T, typename E>
std::vector<std::vector<double>> SCM<T, E>::refinement_error_by_node() const {
  return patch.refinement_error_by_node();
}

// Diagnostic probe: the per-node competition (light) error for one species --
// the per-step sample that collect_competition_errors() accumulates. Exposed
// mainly so tests / the vignette can reconstruct the error signal by hand.
template <typename T, typename E>
std::vector<double>
SCM<T, E>::r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const {
  // The per-node error is scaled by the total competition effect inside the
  // patch-level call below: it computes compute_competition(0.0) -- which
  // already divides by patch area -- and passes it through as the scaling
  // argument. The live schedule-refinement collector
  // (Patch::collect_competition_errors) reconstructs the signal via this same
  // path, so no extra area scaling is needed here to keep them consistent
  // (resolves the scaling question in #478).
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_compute_competition_effect_error_by_node_for_species_i(idx);
}

template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule(NodeSchedule x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  node_schedule = x;

  // Update here so that extracting Parameters reflects the new schedule,
  // keeping Parameters self-sufficient.
  parameters.node_schedule_times = node_schedule.get_times();
}

template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule_times(
    std::vector<std::vector<double>> x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  node_schedule.set_times(x);
  parameters.node_schedule_times = x;
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::r_ode_times() const {
  return solver.times();
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::r_ode_step_sizes() const {
  return solver.step_sizes();
}

template <typename T, typename E>
Rcpp::List SCM<T, E>::r_store_trajectory() {
  const trajectory rec = store_trajectory();
  Rcpp::List ret(rec.size());
  for (size_t i = 0; i < rec.size(); ++i) {
    ret[i] = Rcpp::List::create(
        Rcpp::_["time"] = rec[i].time,
        Rcpp::_["step_size"] = rec[i].step_size,
        // Which of the two kinds a row is.
        // odelia records an insertion; to a reader of plant it is the
        // introduction that made it.
        Rcpp::_["introduction"] = rec[i].insertion,
        Rcpp::_["state"] = rec[i].state);
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::census() const {
  const auto& metrics = metrics_of<patch_type>();
  std::vector<double> ret;
  ret.reserve(metrics.size());
  for (const census_metric<T>& metric : metrics) {
    ret.push_back(odelia::util::to_passive(census_sum(patch, metric)));
  }
  return ret;
}

// The census differentiated with respect to the state AND the traits, from ONE
// recording of one metric algebra, so the seam between the halves is the solver's
// own rather than a second one written here.
//
// The order the halves are written in is the reason they are one recording rather
// than two calls: the traits are placed first and the state loaded after, so a
// quantity the state determines is derived at the traits the recording registered.
// Loading the state first derives it at the values they had before.
//
// set_state_and_boundary rebuilds the environment and the boundary node from the state
// it is given. Both are on the census's path -- the boundary node is the
// reduction's lower grid point and is not ODE state -- so the recording must carry
// that rebuild. Loading the state without it leaves the boundary node at the values
// it was copied with, and its whole contribution to the seed is then exactly zero
// with nothing thrown. Loading it with set_ode_state alone leaves the condition at
// its first evaluation, which is not the one census() reads.
//
// One patch, one tape, one recording, and a seed per metric. The recording does
// not depend on which metric is being asked for -- it writes every metric into y
// and only the seed picks one out -- so a recording per metric was a recording
// repeated.
//
// What made that repetition look necessary is real and is worth stating, because
// it is the trap next door. Clearing a tape returns its derivative-slot counter
// to zero, so an active value built outside a sweep loop and read inside it
// refers, after the first clear, to a slot that now belongs to something else.
// Measured when that was live: the second and third metrics' seeds were wrong by
// three orders and their heartwood columns read exactly zero -- the first metric
// correct and lending its credibility to the rest. The answer is not a patch per
// metric, it is to clear once and record once, which is what sweeping a batch
// does: the clear happens before the recording, and between sweeps only the
// derivative slots are returned to zero.
template <typename T, typename E>
census_rows SCM<T, E>::census_state_and_trait_rows() const {
  require_birth_date_coordinate("census_state_and_trait_rows");
  using scalar = odelia::ode::active_scalar<double>;

  std::vector<double> state(patch.ode_size());
  patch.ode_state(state.begin());

  const size_t n_metric = metrics_of<patch_type>().size();
  census_rows ret;
  ret.trait.assign(n_metric, patch.trait_adjoint_size());

  auto reduce = [&](auto& active, typename std::vector<scalar>::const_iterator x,
                    std::vector<scalar>& y) -> void {
    // The traits carry their derivative from where they sit on the strategy; this
    // buffer is the state and nothing else.
    active.set_state_and_boundary(x, time());
    const auto& metrics = metrics_of<std::decay_t<decltype(active)>>();
    for (size_t m = 0; m < metrics.size(); ++m) {
      y[m] = census_sum(active, metrics[m]);
    }
  };
  // One recording, so the tape and the rebind are odelia's to make. A sweep taking
  // many hands in the tape it holds for the descent.
  odelia::ode::state_and_parameter_adjoints(
      patch, state, odelia::ode::adjoint_rows::all_rows(n_metric), reduce,
      ret.state, ret.trait);
  return ret;
}

// The census twice per trait, at the state held, with the strategy moved in place.
// See the declaration for why it perturbs rather than rebuilds.
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::census_trait_difference(double rel) {
  require_birth_date_coordinate("census_trait_difference");
  const auto& metrics = metrics_of<patch_type>();
  const size_t n_metric = metrics.size();
  const size_t n_state = patch.ode_size();

  std::vector<double> state(n_state);
  patch.ode_state(state.begin());
  const double time_ = time();

  // The patch answers for this order; walking the species here would be free to
  // walk it differently.
  const std::vector<typename T::value_type*> pars = patch.ad_parameters();

  // The state is re-set on every evaluation, which is what makes the moved trait
  // reach the quantities a state determines -- the boundary node among them.
  auto census_at = [&](std::vector<double>& out) -> void {
    patch.set_state_and_boundary(state.begin(), time_);
    out.clear();
    for (const census_metric<T>& metric : metrics) {
      out.push_back(odelia::util::to_passive(census_sum(patch, metric)));
    }
  };

  std::vector<std::vector<double>> ret(n_metric,
                                       std::vector<double>(pars.size(), 0.0));
  std::vector<double> up, dn;
  for (size_t c = 0; c < pars.size(); ++c) {
    const double base = odelia::util::to_passive(*pars[c]);
    const double h = std::max(std::abs(base) * rel, rel);
    *pars[c] = base + h;
    census_at(up);
    *pars[c] = base - h;
    census_at(dn);
    *pars[c] = base;
    for (size_t m = 0; m < n_metric; ++m) {
      ret[m][c] = (up[m] - dn[m]) / (2.0 * h);
    }
  }
  // Leave the patch where it was found, so this call is repeatable beside the
  // recording that shares its state.
  census_at(up);
  return ret;
}

// Seed lambda on the states the census reads at T, then run the reverse pass
// back over the recorded steps. The trait adjoints accumulate across every
// cohort and every step, so the accumulator is cleared once per metric and read
// once the sweep is done. It is the solver's system that accumulates: `patch` is
// a snapshot the run copies out, and reading its accumulator gives zeros.
template <typename T, typename E>
census_gradient
SCM<T, E>::census_trait_gradient(const std::vector<size_t>& extra_stops,
                                 const std::vector<std::string>& which_metrics) {
  require_birth_date_coordinate("census_trait_gradient");
  // Which rows to sweep, resolved against the strategy's own list before
  // anything runs, so the shape of the answer is known on the refusal path too.
  // Named rather than positional: a caller indexing by position gets a different
  // metric's gradient when the list changes, and nothing says so.
  const auto& metrics = metrics_of<patch_type>();
  std::vector<size_t> rows;
  if (which_metrics.empty()) {
    rows.resize(metrics.size());
    for (size_t m = 0; m < rows.size(); ++m) {
      rows[m] = m;
    }
  } else {
    for (const std::string& want : which_metrics) {
      size_t at = metrics.size();
      for (size_t m = 0; m < metrics.size(); ++m) {
        if (want == metrics[m].name) {
          at = m;
          break;
        }
      }
      if (at == metrics.size()) {
        std::string known;
        for (const census_metric<T>& metric : metrics) {
          known += known.empty() ? "" : ", ";
          known += metric.name;
        }
        util::stop("census_trait_gradient: this model has no census metric `" +
                   want + "`; it has " + known);
      }
      rows.push_back(at);
    }
  }
  // The sweep needs the state at every accepted step, and reads them off the
  // solver itself. store_trajectory() repeats the run to get them unless
  // record_trajectory kept them the first time, and either way it may run, so the
  // seeds below are taken after it.
  store_trajectory();

  patch_type& live = solver.get_system_ref();

  // Cleared here so that a previous call's degeneracy cannot refuse this one.
  clear_recorded_refusal(live);
  // The seeds and the direct term are on the gradient path too, so a refusal is
  // reachable before the sweep as well as inside it. Polled at both, and the
  // shape of the answer is known from what the caller asked for and the patch's
  // trait width either way -- so a refusal before the sweep costs the sweep and
  // nothing else.
  census_rows both = census_state_and_trait_rows();
  const size_t width = live.trait_adjoint_size();
  refusal why = recorded_refusal(live);
  odelia::ode::adjoint_rows trait_adjoint;
  census_gradient ret;
  if (!why.happened()) {
    // Every metric's sweep visits the same trajectory and differs only in its
    // seed, so they are carried TOGETHER: a block is recorded once and swept once
    // per metric, where the loop this replaces recorded it once per metric. The
    // recording is a model evaluation and a sweep is arithmetic, so the second and
    // third metrics were costing what the first did and now cost almost nothing.
    odelia::ode::adjoint_rows lambda = both.state.select(rows);

    // The accumulator the sweep adds into, owned here for the length of the sweep
    // rather than kept on the patch between calls. Every writer reaches it through
    // the driver, so a row that does not match the seeds is a length mismatch.
    //
    // It STARTS at the direct term, which is what the total derivative is: the
    // census reads the traits as well as the state, so a metric's gradient is that
    // reading plus the sum over the trajectory. Seeded rather than added at the end
    // because a term added last is a term that can be left out, and a gradient
    // missing it is a plausible number rather than an error.
    trait_adjoint = both.trait.select(rows);
    util::check_length(trait_adjoint.width(), width);

    // One range per width, highest first, narrowing across each introduction and
    // transposing the map that took it. The solver owns that walk: what is left
    // here is the census the sweep is seeded from.
    //
    // ⚠️ CAUGHT BY TYPE, AND NOT BY `runtime_error`. A descent that leaves the
    // range a double holds is not a bug in the sweep and not a state the model
    // has no meaning for, and this is the one place that can say so: the numbers
    // it would otherwise hand back are the overflow's, one per input, with
    // nothing to distinguish them from an answer. A broader catch here would read
    // a genuine length mismatch in the walk the same way.
    //
    // The refusal names NO SPECIES, and that is the honest grain: what overflowed
    // is an intermediate of one recording spanning every cohort in every stage,
    // so nothing finer has a component to attribute it to -- which is why it is
    // built here rather than recorded on a strategy the way a leaf's is.
    size_t ranges = 0;
    try {
      ranges = solver.solve_adjoint(lambda, trait_adjoint, extra_stops);
    } catch (const odelia::util::AdjointRangeError& e) {
      why = refusal{std::string("TF24 gradient: ") + e.what(), -1};
    }
    if (!why.happened()) {
      std::vector<std::vector<double>> at_first_state = lambda.to_rows();

      // Polled again, because the sweep is where most refusals are raised. The
      // two diagnostics cross to the answer only here, on the one path that has
      // a sweep to describe -- so a refusal leaves them at their defaults rather
      // than being cleared back to them.
      why = recorded_refusal(live);
      if (!why.happened()) {
        ret.ranges = ranges;
        ret.at_first_state = std::move(at_first_state);
      }
    }
  }

  // A refusal costs every metric, and the grain is forced rather than chosen: the
  // row that could not be supplied is an intermediate of a recording spanning six
  // stages and every cohort in them, so no seed carries a component to attribute
  // it to.
  ret.gradient.reserve(rows.size());
  ret.why.reserve(rows.size());
  for (size_t m = 0; m < rows.size(); ++m) {
    if (why.happened()) {
      // A sum has no defined value with an undefined term, so the whole metric
      // goes -- not the cohort's column and not the parameter's entry. The
      // numbers are not-a-number rather than absent so that a caller indexing by
      // position still finds the shape it expects.
      ret.gradient.push_back(
          std::vector<double>(width, std::numeric_limits<double>::quiet_NaN()));
      ret.why.push_back(why);
      continue;
    }
    // Nothing is read against a declaration. A parameter no gradient exists for
    // cannot be asked for, so every column here carries a number the sweep
    // computed, and an exact zero in one is the sweep's answer.
    ret.gradient.emplace_back(trait_adjoint[m].begin(), trait_adjoint[m].end());
    ret.why.emplace_back();
  }

  return ret;
}


// The reference the trajectory sweep is checked against. It is a tangent of the
// same forward source: exact, with no step size of its own and no truncation,
// and it traverses both reductions and the introduction boundary while none of
// the transposes under test are on its path.
//
// The recorded step sizes are replayed rather than the times, and rather than a
// controller of its own. A tangent run left to choose its own steps
// differentiates the controller, which the model does not contain -- and a size
// differenced back out of two recorded times is not the size that was taken,
// since fl(fl(t + h) - t) != h.
template <typename T, typename E>
std::vector<double>
SCM<T, E>::census_trait_tangent(const std::vector<double>& direction,
                                std::vector<double>& value) {
  require_birth_date_coordinate("census_trait_tangent");

  const trajectory rec = store_trajectory();
  patch_type& live = solver.get_system_ref();
  odelia::ode::be_at_step(live, rec, 0);
  auto active = live.template rebind_from<tangent>();

  // Seeded before the state is set: the quantities a state determines read the
  // parameters, and would otherwise be derived at the unseeded values.
  size_t at = 0;
  for (tangent* p : active.ad_parameters()) {
    if (at >= direction.size()) {
      util::stop("census_trait_tangent: one weight per trait, species-major");
    }
    seed_direction(*p, direction[at++]);
  }
  util::check_length(direction.size(), at);
  std::vector<tangent> x0(rec[0].state.size());
  for (size_t i = 0; i < x0.size(); ++i) {
    x0[i] = rec[0].state[i];
  }
  active.set_ode_state(x0.begin(), rec[0].time);

  odelia::ode::Solver<decltype(active)> forward(active, make_ode_control(control));
  forward.set_collect(false);
  forward.set_state_from_system();

  // From row 0: the state seeded above is that row's, at its own width, and any
  // introduction the run made is a row of its own above it.
  forward.advance_recorded(odelia::ode::program_from(
      rec, 0, {forward.time(), std::numeric_limits<double>::quiet_NaN()}));

  // Leave the double system where the run left it, so this call is repeatable
  // beside the sweep that shares its trajectory.
  odelia::ode::be_at_step(solver.get_system_ref(), rec, rec.size() - 1);

  const auto& reached = forward.get_system_ref();
  const auto& metrics = metrics_of<std::decay_t<decltype(reached)>>();
  std::vector<double> ret;
  ret.reserve(metrics.size());
  value.clear();
  value.reserve(metrics.size());
  // Reduced ONCE per metric, and its value read off the same tangent as its
  // derivative. Two reductions would be two evaluations of one function, which is
  // the shape this whole check exists to catch elsewhere.
  for (const auto& metric : metrics) {
    const tangent reached_metric = census_sum(reached, metric);
    ret.push_back(derivative_along(reached_metric));
    value.push_back(odelia::util::to_passive(reached_metric));
  }
  return ret;
}


template <typename T, typename E>
template <class Scalar, class Seed>
std::vector<Scalar> SCM<T, E>::replay_initial_state(size_t from_range,
                                                    Seed seed) {
  const trajectory rec = store_trajectory();
  patch_type& live = solver.get_system_ref();

  std::vector<double> base;
  size_t start = 0;
  const double t0 = odelia::ode::state_at_range(live, rec,
                                                  from_range, base, start);

  auto active = live.template rebind_from<Scalar>();
  std::vector<Scalar> x0(base.size());
  seed(x0, base);
  active.set_ode_state(x0.begin(), t0);

  odelia::ode::Solver<decltype(active)> forward(active, make_ode_control(control));
  forward.set_collect(false);
  forward.set_state_from_system();

  // `start` is the row the System was put on, so the program is what follows it
  // whether that row is the beginning of the recording or an introduction already
  // applied.
  forward.advance_recorded(odelia::ode::program_from(
      rec, start, {forward.time(), std::numeric_limits<double>::quiet_NaN()}));

  // Leave the double system where the run left it, so this call is repeatable
  // beside the sweep that shares its trajectory.
  odelia::ode::be_at_step(live, rec, rec.size() - 1);

  const auto& reached = forward.get_system_ref();
  const auto& metrics = metrics_of<std::decay_t<decltype(reached)>>();
  std::vector<Scalar> out;
  out.reserve(metrics.size());
  for (const auto& metric : metrics) {
    out.push_back(census_sum(reached, metric));
  }
  return out;
}

template <typename T, typename E>
std::vector<double>
SCM<T, E>::census_initial_state_tangent(const std::vector<double>& direction,
                                        std::vector<double>& value,
                                        size_t range) {
  require_birth_date_coordinate("census_initial_state_tangent");

  const std::vector<tangent> reached =
    replay_initial_state<tangent>(range,
      [&](std::vector<tangent>& x0,
          const std::vector<double>& base) -> void {
        util::check_length(direction.size(), base.size());
        for (size_t i = 0; i < x0.size(); ++i) {
          x0[i] = base[i];
          seed_direction(x0[i], direction[i]);
        }
      });

  std::vector<double> ret;
  ret.reserve(reached.size());
  value.clear();
  value.reserve(reached.size());
  for (const tangent& metric : reached) {
    ret.push_back(derivative_along(metric));
    value.push_back(odelia::util::to_passive(metric));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double>
SCM<T, E>::census_initial_state_replay(const std::vector<double>& state0,
                                       size_t range) {
  require_birth_date_coordinate("census_initial_state_replay");
  return replay_initial_state<double>(range,
    [&](std::vector<double>& x0, const std::vector<double>& base) -> void {
      util::check_length(state0.size(), base.size());
      x0 = state0;
    });
}


} // namespace plant

#endif
