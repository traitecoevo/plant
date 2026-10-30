// -*-c++-*-
#ifndef SPECIES
#define SPECIES

#include <vector>
#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <plant/util.h>
#include <plant/canopy_shape.h>
#include <plant/environment.h>
#include <odelia/ode_interface.hpp>
#include <plant/node.h>
#include <plant/species_base.h>
#include <plant/census.h>
#include <odelia/drivers.hpp>
#include <plant/with_slope.h>

namespace plant {

// This is purely for running the deterministic model. It shares its storage and
// ODE plumbing with the stochastic species through SpeciesBase (species_base.h);
// the size-density-specific machinery below (density-weighted competition,
// survival-weighted rates, lifetime fitness, schedule-refinement error) stays
// here.

// A strategy that says where its height growth stops: what the absorbing top
// layer of a continuous patch is pinned at (Control::absorb_epsilon, #519).
template <typename T>
concept HasAsymptoticHeight = requires(const T& s, double epsilon) {
  { s.asymptotic_height(epsilon) } -> std::convertible_to<double>;
};

// A strategy whose cohorts each carry a private root-zone soil column (TF24w,
// #681): the patch rates it in two extra passes, applies rain pulses to every
// column, and refuses the height coordinate, the gradient path and the
// stochastic patch for it.
template <typename T>
concept HasPrivateWater = requires {
  { T::n_water_layers } -> std::convertible_to<size_t>;
};

template <typename T, typename E>
class Species : public SpeciesBase<Species<T, E>, T, E, Node<T, E>> {
  typedef SpeciesBase<Species<T, E>, T, E, Node<T, E>> base_type;
public:
  using value_type = typename T::value_type;

  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef Node<T,E> node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Species(strategy_type s);
  // Build on a strategy that is already prepared; see SpeciesBase.
  explicit Species(strategy_type_ptr s);

  // ODE plumbing and the per-element serialisers are inherited from SpeciesBase
  // and iterate all nodes (the deterministic model has no notion of "dead").
  // The nodes' blocks, then the top layer's when it is on (top_ode_size()).
  size_t ode_size() const { return base_type::ode_size() + top_ode_size(); }
  template <typename It> It ode_state(It it) const {
    return top_ode_state(base_type::ode_state(it));
  }
  template <typename It> It ode_rates(It it) const {
    return top_ode_rates(base_type::ode_rates(it));
  }
  using base_type::get_node_state;
  using base_type::get_node_aux;
  typename std::vector<node_type>::iterator node_begin() { return nodes.begin(); }
  typename std::vector<node_type>::iterator node_end() { return nodes.end(); }
  typename std::vector<node_type>::const_iterator node_begin() const { return nodes.begin(); }
  typename std::vector<node_type>::const_iterator node_end() const { return nodes.end(); }

  // Every soil column this species holds, for the private-water passes (#681):
  // each node's individual with the ground it stands for -- its trapezium
  // weight on the birth-date grid times its density, individuals per m2 of
  // patch -- then the boundary node, which closes the grid as it does in every
  // other reduction here, then the absorbed top layer at its density. Retired
  // nodes come through at weight zero. The grid is the one consumption_rate()
  // integrates over, so a column's share of the patch is the share of the water
  // it draws.
  //
  // The third argument is the weight's own relative rate of change, which the
  // land-transfer pass adds to a share's rate. The boundary node sits at the
  // current time, so the youngest cohort's half-interval grows at one half a
  // year per year; every other weight is fixed between introductions, and an
  // introduction leaves every weight continuous (the new boundary stands where
  // the old one did). The top layer is a point mass.
  //
  // The fourth says which kind of column it is: 0 a node, 1 the boundary node,
  // 2 the top layer. The boundary node's state is re-seeded every evaluation
  // rather than integrated, so it owns no ground: the passes count its uptake
  // (it closes the quadrature grid, as it does in consumption_rate) but charge
  // that uptake to the pool, and give it no share.
  template <class F> void for_each_water_column(F&& f) {
    using odelia::util::to_passive;
    const bool birth_date = density_in_birth_date();
    const size_t m = nodes.size();
    if (m > 0) {
      std::vector<double> x;
      x.reserve(m + 1);
      for (auto& c : nodes) {
        x.push_back(abscissa_of(c, birth_date));
      }
      x.push_back(abscissa_of(new_node, birth_date));
      for (size_t k = 0; k < m; ++k) {
        const double left = k > 0 ? x[k] - x[k - 1] : 0.0;
        const double right = x[k + 1] - x[k];
        const double w = (left + right) / 2.0;
        const double w_dt = (k + 1 == m) ? 0.5 : 0.0;
        f(nodes[k].individual, w * to_passive(nodes[k].get_density()),
          w > 0.0 ? w_dt / w : 0.0, 0);
      }
      const double w = (x[m] - x[m - 1]) / 2.0;
      f(new_node.individual, w * to_passive(new_node.get_density()), 0.0, 1);
    }
    if (top_on() && top_density_ > 0.0) {
      f(top_, to_passive(top_density_), 0.0, 2);
    }
  }

  size_t size() const;
  void clear();
  void introduce_new_node();
  // Introduce a node, stamping it with the introduction time and patch-age
  // density at birth (called by Patch, which knows the time and disturbance).
  // The canonical insertion: a node carries three numbers the ODE state does not,
  // and all three are set here so no caller can push a node missing one.
  void introduce_new_node(double time, double patch_density,
                          double pr_patch_survival);
  // Drop the node introduce_new_node pushed last, which is the newest: the width
  // a reverse sweep needs before an introduction.
  void remove_newest_node();

  value_type height_max() const;
  // The query height is a knot position on the interpolant's own grid, which is
  // double (see ResourceSpline::rebuild_spline); what the cohorts put into the
  // field arrives through the node contributions summed over below.
  value_type compute_competition(const value_type& height) const;

  // The reduction and its vertical derivative from one traversal, so each node's
  // u^eta is evaluated once and the two sums add their terms in the same order.
  // The first entry equals compute_competition(height) bit for bit.
  with_slope<value_type>
  compute_competition_and_slope(const value_type& height) const;

  // The reduction stopped before its closing trapezium: the sums so far and the
  // last sample, which is everything one more interval needs. Holding it lets the
  // field WITH the boundary interval be formed from the field without it in one
  // operation rather than by a second reduction over every node -- the two differ
  // only in that trapezium, and the boundary node it needs is not known until the
  // field without it exists.
  struct competition_split {
    // The running sums, and the sample at x1 that one more interval needs.
    with_slope<value_type> sum{value_type(0.0), value_type(0.0)};
    with_slope<value_type> at_x1{value_type(0.0), value_type(0.0)};
    // The abscissa is a position (abscissa_of takes it passive), so a width the
    // closing forms from it is a position too and not state.
    double x1{0.0};
    bool closes{false};

    // The field the nodes were rated in. close_competition_and_slope() is the
    // same halved sums with the boundary interval added first.
    with_slope<value_type> without_boundary() const {
      return {sum.value / 2, sum.slope / 2};
    }

    template <class F>
    void for_each_active(F&& f) {
      odelia::ode::visit_active(f, sum, at_x1);
    }
  };
  // A split at every height of a set, from ONE pass over the nodes.
  //
  // Q is a polynomial in w = (z / H)^eta and w separates into z^eta times H^-eta,
  // so the trapezium's intervals carry three running sums -- one per power -- and
  // every height reads a prefix of them. That makes the build linear in nodes plus
  // heights where a height-by-height walk is their product, and inside a recording
  // the operation count IS the tape.
  void field_splits(const std::vector<value_type>& heights,
                    std::vector<competition_split>& out) const;
  // The top layer as a point mass, added to a split's running sums (twice its
  // contribution, since every reader halves them).
  void add_top(competition_split& c, const value_type& height) const;
  // The inclusive reduction, from a split taken at the same height. Bit-identical
  // to compute_competition_and_slope(height) at the boundary node it is closed
  // with.
  with_slope<value_type>
  close_competition_and_slope(const competition_split& c,
                              const value_type& height) const;

  // Evaluate the inflow boundary condition in the environment passed. Split out
  // of compute_rates() so the field build owns it and the field stops reading a
  // density carried from the previous evaluation.
  void compute_boundary_node(const environment_type& environment,
                             double pr_patch_survival,
                             const value_type& birth_rate) {
    new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
  }

  // The trapezium integral of n_k psi(state_k) over the size distribution, with
  // n_k = exp(l_k). A census is a quadrature of a density, so the grid is the
  // coordinate that density is carried in; taking gaps in any other variable
  // integrates one density against another's spacing. The inflow boundary node
  // closes the grid, being the lower limit of the distribution.
  // The metric integrated over this species' size distribution.
  value_type census_integral(const census_metric<T>& metric) const;

  // Every active value this species holds: its nodes, the boundary node, the
  // cached height scan, and the strategy it owns -- which is reached here and
  // nowhere else, because every node shares it.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, nodes, new_node, height_scan_cache.h_max,
                              top_, top_density_, top_weight_,
                              top_weight_birth_, top_seeds_, top_seeds_birth_,
                              top_density_dt_, top_weight_dt_,
                              top_weight_birth_dt_, top_seeds_dt_,
                              top_seeds_birth_dt_, strategy);
  }

  // Whether the decreasing-height node ordering still holds (see height_max()).
  bool heights_are_decreasing() const;

  // The tallest node height and whether the heights are ordered, from a single
  // pass. compute_competition() needs both on every call, and walking the heights
  // twice was measurably slower on FF16 (~5% on the SCM benchmark) than the
  // O(1) nodes.front() it replaced.
  // h_max is the canopy top, a position; decreasing is a comparison outcome and
  // structural, so it stays bool whatever the heights are made of.
  struct HeightScan { value_type h_max; bool decreasing; };
  // Cached: heights change only when the ODE state is set or a node is
  // introduced/cleared, whereas compute_competition() is called once per spline
  // knot, so this is hundreds of calls per change. Every mutator invalidates.
  const HeightScan& scan_heights() const;

  // Setting the ODE state rewrites every node's height, so the cached scan goes
  // with it. Shadows (rather than uses) the SpeciesBase version for that reason.
  template <typename It> It set_ode_state(It it) {
    invalidate_height_scan();
    return set_top_ode_state(base_type::set_ode_state(it));
  }

  // * The absorbing top layer (Control::absorb_epsilon, #519)
  //
  // One cohort per species, pinned at the strategy's asymptotic height, into
  // which a cohort that has grown there is absorbed at a node introduction. The
  // absorbed cohort stays in the node list as a zero-density grid point, exactly
  // as a retired one does, and the share of every trapezium it held -- its
  // quadrature weight times its density -- moves into the top layer. Keeping
  // the grid point is what makes the transfer exact: the neighbours' intervals
  // are untouched, so total density is conserved to rounding at the instant,
  // and the node structure never changes, so a pinned run or an export needs no
  // record of which cohorts went.
  //
  // Its block is the pinned individual's state (height fixed, the rest
  // integrated: TF24's storage, the heartwood), then five numbers: the density
  // it carries, the fitness weight it carries without and with the birth-rate
  // scaling (sum of quadrature weight x patch density at birth x S_D x survival
  // at absorption / patch survival at birth, and the same x birth rate at
  // birth), and the seed output each has accrued since. The lifetime output an
  // absorbed cohort has already accrued stays on its node.
  size_t top_ode_size() const {
    return top_on() ? individual_type::ode_size() + 5 : 0;
  }
  bool top_on() const { return control().absorb_epsilon > 0.0; }
  // Every number the layer carries is finite. Its readers test density > 0, so
  // a NaN layer would otherwise read as an empty one.
  bool top_finite() const {
    using odelia::util::to_passive;
    for (size_t i = 0; i < individual_type::ode_size(); ++i) {
      if (!util::is_finite(to_passive(top_.state(i)))) {
        return false;
      }
    }
    for (const value_type* v : {&top_density_, &top_weight_, &top_weight_birth_,
                                &top_seeds_, &top_seeds_birth_}) {
      if (!util::is_finite(to_passive(*v))) {
        return false;
      }
    }
    return true;
  }
  // The height the layer sits at, from the strategy's allocation to reproduction.
  double absorb_height() const;
  // Absorb every live cohort, bar the youngest, at or above absorb_height(), and
  // say how many. Under an arrival weighting a cohort the window weights is
  // never absorbed, so every cohort the fitness integral reads keeps its own
  // lifetime output.
  size_t absorb_tall(bool arrival_weighting);
  // The seed output the layer has accrued, for the fitness integral.
  value_type top_seeds(bool birth_scaled) const {
    if (!top_on()) {
      return value_type(0.0);
    }
    return birth_scaled ? top_seeds_birth_ : top_seeds_;
  }
  double top_density() const {
    return top_on() ? odelia::util::to_passive(top_density_) : 0.0;
  }
  double top_height() const {
    return odelia::util::to_passive(top_.state(HEIGHT_INDEX));
  }
  Rcpp::List r_top_layer() const;
  // The layer's density and weights decay linearly, and at the bottom of the
  // subnormal range a decay step rounds back to the value it started from, so a
  // layer whose inflow has stopped never empties by itself: it would keep the
  // canopy top at the absorbing height and a pinned tree rated at every stage.
  // Below the smallest normal double a value carries nothing, so it is put to
  // zero. The seeds the layer accrued are output, not state that decays, and
  // stay. Called at each introduction instant, with absorption.
  // An event that removes a fraction of the stand at the layer's height removes
  // that fraction of the trees the layer holds: its density and both fitness
  // weights scale together. Its seeds are past output and stay. Returns the
  // density removed.
  double scale_top(double phi) {
    const double before = odelia::util::to_passive(top_density_);
    top_density_ *= phi;
    top_weight_ *= phi;
    top_weight_birth_ *= phi;
    return before - odelia::util::to_passive(top_density_);
  }
  void flush_top() {
    const double tiny = std::numeric_limits<double>::min();
    for (value_type* v : {&top_density_, &top_weight_, &top_weight_birth_}) {
      if (odelia::util::to_passive(*v) < tiny) {
        *v = value_type(0.0);
      }
    }
  }
  void compute_rates(const environment_type& environment, double pr_patch_survival, const value_type& birth_rate);
  void compute_top_rates(const environment_type& environment, double pr_patch_survival);

  std::vector<double> net_reproduction_ratio_by_node() const;
  // Per-node lifetime offspring, weighted by patch-age density and S_D.
  std::vector<double> net_reproduction_ratio_by_node_weighted() const;
  // The same, at the scalar the nodes carry: what a differentiated fitness
  // reads. The double form above is this with the derivative dropped.
  std::vector<value_type> weighted_fecundity_by_node() const;
  // Introduction times of each node (the integration x-axis for fitness).
  std::vector<double> node_times() const;

  // The boundary node's birth date is *now*, but compute_initial_conditions()
  // (which stamps it) runs inside compute_rates(), which the ODE stepper calls
  // after set_ode_state() has already rebuilt the environment. Reading the stamp
  // during that rebuild would therefore pick up the previous derivs call's time
  // and shorten the boundary trapezium by one Runge-Kutta stage. Patch refreshes
  // it before building the profile. Harmless on the height path, where the
  // boundary abscissa is the constant initial height.
  void set_new_node_birth_date(double time) {
    new_node.set_introduction_time(time);
  }

  // Two nodes sharing a birth date give a zero-width trapezium interval, so the
  // birth-date quadrature silently loses them. Cannot happen for a scheduled
  // run (introduction times are distinct by construction) but can for a patch
  // initialised or imported without per-node times.
  bool birth_dates_are_distinct() const;

  // Which coordinate this species' size distribution is carried in: Control's
  // node_density_coordinate with "auto" resolved against the model's own
  // requirement (see control.h). Exposed so Patch can check every species agrees
  // before summing their contributions.
  bool density_in_birth_date() const {
    return plant::density_in_birth_date(control(),
                                        T::density_in_birth_date_default);
  }

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  // (ode_size/set_ode_state/ode_state/ode_rates come from SpeciesBase.)
  size_t aux_size() const;

  void resize_consumption_rates(int i);
  value_type consumption_rate(int i) const;

  template <typename It> It ode_aux(It it) const;

  Rcpp::NumericMatrix r_get_state() const;

  // * R interface
  std::vector<double> r_heights() const;
  std::vector<double> r_heights_rev() const;
  void r_set_heights(std::vector<double> heights);
  const node_type& node_at(size_t i) const {return nodes[i];}
  const node_type& r_new_node() const {return new_node;}
  std::vector<node_type> r_nodes() const {return nodes;}
  const node_type& r_node_at(util::index idx) const {
    return nodes[idx.check_bounds(size())];
  }

  // Do this with set_ode_state, using an iterator?
  /* double state(int i) const { return vars.state(i); } */

  /* double rate(int i) const { return vars.rate(i); } */

  /* void set_state(int i, double v) { */
  /*   vars.set_state(i, v); */
  /* } */


  // These are used to determine the degree of node refinement.
  std::vector<double> r_compute_competition_effect_by_nodes() const;
  std::vector<double> r_compute_competition_effect_by_nodes_error(double scal) const;

  // Per-node size density, **always** as a density in height whichever
  // coordinate the solver carried it in, so downstream code (tidy_outputs.R's
  // `density`, interpolate_to_heights(), the plots) keeps its meaning. In
  // birth-date coordinates that means dividing the carried quantity by the
  // Jacobian; see height_jacobian(). NA for a node where the Jacobian vanishes,
  // which is where the height density genuinely does not exist.
  std::vector<double> r_log_densities() const;
  // The quantity actually integrated: the density in height on the height path,
  // and the density in birth date (nu) on the birth-date one. Reported
  // alongside r_log_densities() rather than instead of it, because it is the
  // thing whose ODE the solver solves and the thing that cannot go non-finite.
  std::vector<double> r_log_densities_state() const;
  // Per-node rate of change of log density; used to guard against initial
  // conditions whose densities would explode to non-finite values. This is the
  // rate of the *carried* quantity, so on the birth-date path it is -mortality.
  std::vector<double> r_log_density_rates() const;

  // |dh/dtau| per node, with the boundary node appended last: the Jacobian of
  // the change of variables between the two coordinates, N = nu / |dh/dtau|.
  std::vector<double> height_jacobian() const;

  // Per-node birth bookkeeping, exposed so an exported patch state can be
  // re-imported faithfully (see node.h::set_birth_state). node_times() above
  // already returns the per-node introduction times.
  std::vector<double> r_patch_densities() const;
  std::vector<double> r_pr_patch_survival_at_birth() const;
  // Restore birth bookkeeping for imported nodes (resume); the argument lengths
  // must each match the current node count.
  void set_birth_state(const std::vector<double>& times,
                       const std::vector<double>& patch_density,
                       const std::vector<double>& pr_patch_survival);

  std::vector<bool> r_retired() const;

  ExtrinsicDrivers extrinsic_drivers() const {return strategy->extrinsic_drivers;}

  // The prepared strategy this species and its nodes share.
  strategy_type_ptr strategy_ptr() const {return this->strategy;}

private:
  // The reduction, over the nodes in ascending abscissa. `order` names that order
  // where the node list is not already in it and is empty where it is; the early
  // exit is the decreasing heights', not this parameter's.
  competition_split reduce_competition(const value_type& height,
                                       const std::vector<std::size_t>& order) const;
  // The node positions in ascending height abscissa, tallest first, for where the
  // heights are no longer ordered: the height path's quadrature grid and the
  // birth-date path's prefix.
  // Positions only: nothing here evaluates a contribution, so no width the
  // reduction forms out of it can carry a derivative.
  std::vector<std::size_t> ascending_by_abscissa() const;
  // Whether the closing trapezium to the boundary node is taken. Below a node
  // that no longer contributes the density is zero and so is the interval; a
  // single node has no interval but that one, and a birth-date abscissa does not
  // order with the support at all.
  bool closes_on(const value_type& value_at_x1) const {
    return size() == 1 || density_in_birth_date() ||
           value_at_x1 > 0;
  }
  competition_split compute_competition_and_slope_split(const value_type& height) const;

  // Cache for scan_heights(). ⚠️ EVERY PATH THAT CAN CHANGE A NODE HEIGHT MUST
  // CALL invalidate_height_scan(): a stale cache here reports the wrong ordering
  // and so the wrong competition profile, with nothing raised. Check the coverage
  // by asserting cache == freshly-computed on every call across the suite and the
  // scenario gateway.
  HeightScan compute_height_scan() const;
  void invalidate_height_scan() { height_scan_valid = false; }
  mutable HeightScan height_scan_cache{value_type(0.0), true};
  mutable bool height_scan_valid = false;

  // Storage (strategy, nodes) and control() live in SpeciesBase; the
  // using-declarations let the unqualified references below resolve through the
  // dependent base.
  using base_type::nodes;
  using base_type::strategy;
  using base_type::control;
  node_type new_node;
  // The top layer (see top_ode_size()) and its rates.
  individual_type top_;
  value_type top_density_{0.0}, top_weight_{0.0}, top_weight_birth_{0.0},
    top_seeds_{0.0}, top_seeds_birth_{0.0};
  value_type top_density_dt_{0.0}, top_weight_dt_{0.0},
    top_weight_birth_dt_{0.0}, top_seeds_dt_{0.0}, top_seeds_birth_dt_{0.0};
  template <typename It> It top_ode_state(It it) const {
    if (!top_on()) {
      return it;
    }
    for (size_t i = 0; i < individual_type::ode_size(); ++i) {
      util::write_iterator_scalar(it, top_.state(i));
    }
    util::write_iterator_scalar(it, top_density_);
    util::write_iterator_scalar(it, top_weight_);
    util::write_iterator_scalar(it, top_weight_birth_);
    util::write_iterator_scalar(it, top_seeds_);
    util::write_iterator_scalar(it, top_seeds_birth_);
    return it;
  }
  template <typename It> It set_top_ode_state(It it) {
    if (!top_on()) {
      return it;
    }
    for (size_t i = 0; i < individual_type::ode_size(); ++i) {
      top_.set_state(i, *it++);
    }
    top_density_ = *it++;
    top_weight_ = *it++;
    top_weight_birth_ = *it++;
    top_seeds_ = *it++;
    top_seeds_birth_ = *it++;
    return it;
  }
  template <typename It> It top_ode_rates(It it) const {
    if (!top_on()) {
      return it;
    }
    for (size_t i = 0; i < individual_type::ode_size(); ++i) {
      util::write_iterator_scalar(it, top_.rate(i));
    }
    util::write_iterator_scalar(it, top_density_dt_);
    util::write_iterator_scalar(it, top_weight_dt_);
    util::write_iterator_scalar(it, top_weight_birth_dt_);
    util::write_iterator_scalar(it, top_seeds_dt_);
    util::write_iterator_scalar(it, top_seeds_birth_dt_);
    return it;
  }
  // The nodes' part of field_splits().
  void field_splits_nodes(const std::vector<value_type>& heights,
                          std::vector<competition_split>& out) const;

  // The abscissa every reduction over the size distribution is taken over,
  // increasing as the node list is walked from the tallest down. Heights are
  // negated so that both coordinates increase in the same direction; negation is
  // exact, so the height branch's trapezium widths are bit-identical to
  // differencing the heights themselves. Callers in hot loops read the
  // coordinate once and pass it in.
  //
  // It is a position, so it is read at its value even where the height it comes
  // from is active: a quadrature grid is structure. On the birth-date coordinate
  // that is exact, the date being fixed at birth. On the height coordinate it
  // drops the weights' own channel and nothing supplies it, which is why the
  // reverse pass refuses that coordinate rather than answering on it.
  static double abscissa_of(const node_type& n, bool birth_date) {
    using odelia::util::to_passive;
    return birth_date ? to_passive(n.introduction_time())
                      : -to_passive(n.height());
  }
  std::vector<double> quadrature_abscissae() const {
    std::vector<double> ret;
    ret.reserve(size());
    const bool birth_date = density_in_birth_date();
    for (auto& c : nodes) {
      ret.push_back(abscissa_of(c, birth_date));
    }
    return ret;
  }

  typedef typename std::vector<node_type>::iterator nodes_iterator;
  typedef typename std::vector<node_type>::const_iterator nodes_const_iterator;
};

template <typename T, typename E>
Species<T,E>::Species(strategy_type s)
  : base_type(s),
    new_node(this->strategy),
    top_(this->strategy) {
}

template <typename T, typename E>
Species<T,E>::Species(strategy_type_ptr s)
  : base_type(s),
    new_node(this->strategy),
    top_(this->strategy) {
}

template <typename T, typename E>
size_t Species<T,E>::size() const {
  return nodes.size();
}

template <typename T, typename E>
void Species<T,E>::clear() {
  invalidate_height_scan();
  nodes.clear();
  // Reset the new_node to a blank new_node, too.
  new_node = node_type(strategy);
  // And the top layer: it belongs to a run, so a run that starts again starts
  // without one. A resumed patch loads the exported layer after this.
  top_ = individual_type(strategy);
  top_density_ = top_weight_ = top_weight_birth_ = value_type(0.0);
  top_seeds_ = top_seeds_birth_ = value_type(0.0);
  top_density_dt_ = top_weight_dt_ = top_weight_birth_dt_ = value_type(0.0);
  top_seeds_dt_ = top_seeds_birth_dt_ = value_type(0.0);
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node() {
  invalidate_height_scan();
  // new_node already holds the initial conditions computed against the current
  // environment by the most recent compute_rates() call (see compute_rates ->
  // new_node.compute_initial_conditions above), and the member is refreshed
  // again on the next compute_rates() ready for the following introduction.
  // Recomputing it here would be redundant, and would (wrongly) recompute against
  // the post-introduction environment rather than the environment at the
  // node's introduction time (resolves the recompute question in #478).
  nodes.push_back(new_node);
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual, which will be at least as tall as a seed.
//
// ⚠️ DO NOT RETURN nodes.front() ON THE DECREASING-HEIGHT ORDERING. That ordering
// holds only while height growth is a function of height and the shared
// environment, and TF24's reserve-gated growth makes dh/dt depend on a cohort's
// own storage, so two cohorts born moments apart into a rapidly changing
// environment can cross in height. The front is then a height below the tallest
// living cohort, which truncates the light spline's domain. The cached scan
// already answers this, so asking it here makes the tallest height and the
// ordering one walk rather than two.
template <typename T, typename E>
typename Species<T,E>::value_type Species<T,E>::height_max() const {
  value_type h = nodes.empty() ? new_node.height() : scan_heights().h_max;
  // The top layer is not a node, so the scan does not see it; the canopy top
  // the field spans has to.
  if (top_on() && top_density_ > 0.0 && top_.state(HEIGHT_INDEX) > h) {
    h = top_.state(HEIGHT_INDEX);
  }
  return h;
}

// Are the node heights ordered largest to smallest? See height_max() above for
// why this cannot be assumed. Heights only, so this is cheap relative to the
// per-node crown integrals it guards.
template <typename T, typename E>
bool Species<T,E>::heights_are_decreasing() const {
  return scan_heights().decreasing;
}

template <typename T, typename E>
const typename Species<T,E>::HeightScan& Species<T,E>::scan_heights() const {
  if (!height_scan_valid) {
    height_scan_cache = compute_height_scan();
    height_scan_valid = true;
  }
  return height_scan_cache;
}

// Tallest height and orderedness in one pass over the heights.
template <typename T, typename E>
typename Species<T,E>::HeightScan Species<T,E>::compute_height_scan() const {
  HeightScan ret{value_type(-std::numeric_limits<double>::infinity()), true};
  value_type h_prev = std::numeric_limits<double>::infinity();
  // ⚠️ WITH THE TOP LAYER ON, A FROZEN COHORT IS LEFT OUT OF THE SCAN. Its height
  // froze when it was absorbed (or reached zero density) and younger cohorts
  // grow past it, so counted it breaks the ordering and every field build falls
  // back to a walk per height: measured, a continuous FF16 stand ran 1.6x SLOWER
  // with frozen cohorts counted than with none frozen at all. It contributes
  // exactly zero to every reduction, so the ordering of the cohorts that
  // contribute is the one the fast build and the walk's early exit rest on.
  // Gated on the top layer, so a run that has not opted in scans exactly as it
  // always has.
  const bool skip_retired = freezes_retired(control());
  bool any = false;
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    if (skip_retired && it->retired()) {
      continue;
    }
    any = true;
    const value_type h = it->height();
    if (h > h_prev) {
      ret.decreasing = false;
    }
    if (h > ret.h_max) {
      ret.h_max = h;
    }
    h_prev = h;
  }
  // Every cohort retired: the canopy top is the seed's, as for an empty species.
  if (skip_retired && !any) {
    ret.h_max = new_node.height();
  }
  return ret;
}

// Because of nodes are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all nodes further down
// the list.
//
// NOTE: This is simply performing numerical integration,  via the
// trapezium rule, of the compute_competition with respect to plant
// height.  You'd think that this would be nicer to do in terms of a
// call to an external trapezium integration function, but building
// and discarding the intermediate storage ends up being a nontrivial
// cost.  A more general iterator version might be possible, but with
// the fiddliness around the boundary conditions that won't likely be
// useful.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights.  The integral is not defined, but an empty light
// environment seems appropriate.
//
// NOTE: A similar early-exit condition to the Plant version is used;
// once the lower bound of the trazpeium is zero, we stop including
// individuals.  Working with the boundary node is tricky here,
// because we might need to include that, too: always in the case of a
// single node (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::compute_competition(const value_type& height) const {
  // The value is the fused reduction's first entry, and taking it from there is
  // what makes them equal rather than a test's business. The two walked the same
  // grid with the same early exit and the same closing trapezium, and a value
  // and a slope from sums that associate differently disagree in their last
  // bits, so the agreement had to be asserted over a grid of heights and crown
  // shapes. One walk cannot disagree with itself.
  //
  // The slope costs an extra evaluation per node. Nothing on the field build's
  // path arrives here -- it takes compute_competition_and_slope_split() and
  // closes it -- so this serves the accessors, where the pair was already being
  // computed one call away.
  return compute_competition_and_slope(height).value;
}

// The same trapezium integral as compute_competition(), and alongside it the
// integral of the vertical derivative, from one traversal of the nodes. Both
// sums visit the same nodes in the same order and associate identically, so the
// first entry is compute_competition(height) bit for bit; a check that it is
// lives in test-canopy-methods.R. The early exit and the closing boundary
// trapezium are driven by the value, as they are there.
template <typename T, typename E>
with_slope<typename Species<T,E>::value_type>
Species<T,E>::compute_competition_and_slope(const value_type& height) const {
  competition_split c = compute_competition_and_slope_split(height);
  add_top(c, height);
  return close_competition_and_slope(c, height);
}

template <typename T, typename E>
void Species<T,E>::add_top(competition_split& c, const value_type& height) const {
  if (!top_on() || !(top_density_ > 0.0)) {
    return;
  }
  const with_slope<value_type> fs = top_.compute_competition_and_slope(height);
  c.sum.value += 2.0 * top_density_ * fs.value;
  c.sum.slope += 2.0 * top_density_ * fs.slope;
}

template <typename T, typename E>
void Species<T,E>::field_splits(const std::vector<value_type>& heights,
                                std::vector<competition_split>& out) const {
  field_splits_nodes(heights, out);
  for (std::size_t k = 0; k < heights.size(); ++k) {
    add_top(out[k], heights[k]);
  }
}

// One pass over the nodes for the whole height set. Where the prefix form does not
// hold -- a broken ordering on the height coordinate, or a profile that is not a
// polynomial in w -- every height takes the walk, which is the same branch the
// walk's own early exit rests on.
template <typename T, typename E>
void Species<T,E>::field_splits_nodes(const std::vector<value_type>& heights,
                                      std::vector<competition_split>& out) const {
  out.assign(heights.size(), competition_split());
  if (size() == 0) {
    return;
  }
  const HeightScan& scan = scan_heights();
  const bool birth_date = density_in_birth_date();
  const std::size_t n_moments = strategy->canopy_shape.n_moments();
  if ((!scan.decreasing && !birth_date) || n_moments == 0) {
    for (std::size_t k = 0; k < heights.size(); ++k) {
      out[k] = compute_competition_and_slope_split(heights[k]);
    }
    return;
  }

  using moments = std::array<value_type, CanopyShape<value_type>::max_moments>;
  const std::size_t n = size();

  // Per node, read once for the whole height set rather than once per height: the
  // factor multiplying Q -- which is the node's contribution at height zero,
  // because Q(0) is exactly one for every profile that has this form -- and the
  // powers of its own inverse height.
  std::vector<value_type> scale(n);
  std::vector<moments> mom(n);
  std::vector<double> abscissa(n);
  {
    std::size_t i = 0;
    for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it, ++i) {
      scale[i] = it->compute_competition(0.0);
      if (!util::is_finite(scale[i])) {
        util::stop("Detected non-finite contribution");
      }
      strategy->canopy_shape.crown_moments(1.0 / it->height(), mom[i]);
      abscissa[i] = abscissa_of(*it, birth_date);
    }
  }
  const bool skip_retired = freezes_retired(control());
  auto frozen = [&](std::size_t i) -> bool {
    return skip_retired && nodes[i].retired();
  };
  moments weight, weight_slope;

  if (!scan.decreasing) {
    // Crossed crowns on the birth-date coordinate, which TF24's reserve-gated
    // growth allows. The widths are in birth date and the heights do not move
    // them, so each node carries its own weight, the widths of the two intervals
    // it bounds, and a height reads the prefix of the nodes that reach it,
    // tallest first: below its own height a node adds an exact zero. The walk's
    // sum, associated differently.
    const std::vector<std::size_t> order = ascending_by_abscissa();
    std::vector<moments> prefix(n + 1);
    for (std::size_t j = 0; j < n_moments; ++j) {
      prefix[0][j] = value_type(0.0);
    }
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t o = order[i];
      const double width = (o > 0 ? abscissa[o] - abscissa[o - 1] : 0.0) +
                           (o + 1 < n ? abscissa[o + 1] - abscissa[o] : 0.0);
      for (std::size_t j = 0; j < n_moments; ++j) {
        prefix[i + 1][j] = prefix[i][j] + width * scale[o] * mom[o][j];
      }
    }
    // The heights ascend and the nodes in `order` descend, so the count reaching
    // a height only falls: one merge over both. A retired cohort is passed over
    // where it sits; its density is zero, so its moments are too.
    std::size_t reach = n;
    for (std::size_t k = 0; k < heights.size(); ++k) {
      const value_type& height = heights[k];
      if (scan.h_max < height) {
        continue;  // no node reaches it; the empty split stands
      }
      while (reach > 0 && (frozen(order[reach - 1]) ||
                           nodes[order[reach - 1]].height() < height)) {
        --reach;
      }
      strategy->canopy_shape.height_weights(height, weight);
      strategy->canopy_shape.height_weight_slopes(height, weight_slope);
      competition_split& c = out[k];
      for (std::size_t j = 0; j < n_moments; ++j) {
        c.sum.value += weight[j] * prefix[reach][j];
        c.sum.slope += weight_slope[j] * prefix[reach][j];
      }
      // The closing trapezium runs from the youngest node, wherever its height.
      c.x1 = abscissa[n - 1];
      c.at_x1 = nodes[n - 1].compute_competition_and_slope(height);
      c.closes = closes_on(c.at_x1.value);
    }
    return;
  }

  // The trapezium's intervals, accumulated in the order the walk accumulates them
  // so that a height reads a prefix rather than a re-association across nodes.
  // prefix[i] is the sum over intervals 1..i; prefix[0] is empty.
  std::vector<moments> prefix(n);
  for (std::size_t j = 0; j < n_moments; ++j) {
    prefix[0][j] = value_type(0.0);
  }
  for (std::size_t i = 1; i < n; ++i) {
    const double width = abscissa[i] - abscissa[i - 1];
    for (std::size_t j = 0; j < n_moments; ++j) {
      prefix[i][j] = prefix[i - 1][j] +
                     width * (scale[i - 1] * mom[i - 1][j] + scale[i] * mom[i][j]);
    }
  }

  // `crossing` is the first node the height is above -- the one the walk stops at.
  // The heights ascend and the nodes descend, so it only ever moves one way: one
  // merge over both, rather than a search per height.
  //
  // ⚠️ A FROZEN COHORT IS PASSED OVER, as the walk passes over it. Its height is
  // outside the ordering the scan certified, so it can sit above the height after
  // a live cohort below it; it carries density zero, so its moments are exact
  // zeros in the prefix wherever it sits. The walk stops at the first LIVE cohort
  // below the height, and so does this.
  std::size_t crossing = n;
  for (std::size_t k = 0; k < heights.size(); ++k) {
    const value_type& height = heights[k];
    if (scan.h_max < height) {
      continue;  // no node reaches it; the empty split stands
    }
    while (crossing > 0 &&
           (frozen(crossing - 1) || nodes[crossing - 1].height() < height)) {
      --crossing;
    }
    // Over the frozen cohorts that sit between the last live one above the height
    // and the first live one below it.
    std::size_t stop = crossing;
    while (stop < n && frozen(stop)) {
      ++stop;
    }
    // The walk stops AFTER the interval whose upper node is the first below the
    // height, so the last node it visited is that one -- or the last node of all,
    // where none is below.
    const std::size_t last = stop < n ? stop : n - 1;
    const std::size_t summed = last > 0 ? last - (stop < n ? 1 : 0) : 0;

    strategy->canopy_shape.height_weights(height, weight);
    strategy->canopy_shape.height_weight_slopes(height, weight_slope);
    with_slope<value_type> sum{value_type(0.0), value_type(0.0)};
    for (std::size_t j = 0; j < n_moments; ++j) {
      sum.value += weight[j] * prefix[summed][j];
      sum.slope += weight_slope[j] * prefix[summed][j];
    }

    // The interval the support crosses is the one place the separated form does
    // not hold: its lower node reaches the height and its upper node does not, so
    // the polynomial would evaluate (z / H)^eta above one there rather than the
    // zero the profile has. It is the interval the walk closes by hand too.
    competition_split& c = out[k];
    const with_slope<value_type> fs_last =
      nodes[last].compute_competition_and_slope(height);
    if (stop < n) {
      const with_slope<value_type> fs_above =
        nodes[stop - 1].compute_competition_and_slope(height);
      const double width = abscissa[stop] - abscissa[stop - 1];
      sum.value += width * (fs_above.value + fs_last.value);
      sum.slope += width * (fs_above.slope + fs_last.slope);
    }

    c.sum = sum;
    c.x1 = abscissa[last];
    c.at_x1 = fs_last;
    c.closes = closes_on(c.at_x1.value);
  }
}

// One reduction over the nodes in ascending abscissa. `order` names that order
// where the node list is not in it; empty means in place, which is what the
// decreasing heights buy -- abscissa_of negates height so that the two coincide.
// The early exit belongs to those heights and not to this parameter: below the
// query height a crown contributes an exact zero, so the reduction can stop
// there only while the heights are known to keep falling.
template <typename T, typename E>
typename Species<T,E>::competition_split
Species<T,E>::reduce_competition(const value_type& height,
                                 const std::vector<std::size_t>& order) const {
  const bool birth_date = density_in_birth_date();
  const HeightScan& scan = scan_heights();
  const bool skip_retired = freezes_retired(control());
  const bool permuted = !order.empty();
  const std::size_t n = size();
  auto node_at = [&](std::size_t k) -> const node_type& {
    return nodes[permuted ? order[k] : k];
  };

  const node_type& first = node_at(0);
  const with_slope<value_type> fs1 =
    first.compute_competition_and_slope(height);
  if (!util::is_finite(fs1.value) || !util::is_finite(fs1.slope)) {
    util::stop("Detected non-finite contribution");
  }
  with_slope<value_type> sum{value_type(0.0), value_type(0.0)};
  double x1 = abscissa_of(first, birth_date);
  with_slope<value_type> at_x1 = fs1;

  for (std::size_t k = 1; k < n; ++k) {
    const node_type& node = node_at(k);
    const with_slope<value_type> fs0 =
      node.compute_competition_and_slope(height);
    const double x0 = abscissa_of(node, birth_date);
    if (!util::is_finite(fs0.value) || !util::is_finite(fs0.slope)) {
      util::stop("Detected non-finite contribution");
    }
    sum.value += (x0 - x1) * (at_x1.value + fs0.value);
    sum.slope += (x0 - x1) * (at_x1.slope + fs0.slope);
    x1    = x0;
    at_x1 = fs0;
    // A retired cohort is outside the ordering the scan certified (see
    // compute_height_scan), so its frozen height says nothing about the cohorts
    // after it; it contributes zero, and the walk goes on past it.
    if (scan.decreasing && node.height() < height &&
        !(skip_retired && node.retired())) {
      break;
    }
  }

  competition_split c;
  c.sum = sum;
  c.x1 = x1;
  c.at_x1 = at_x1;
  c.closes = closes_on(at_x1.value);
  return c;
}

// Sorted on the abscissae rather than the nodes, so the scratch holds positions
// and the contributions are evaluated by the one reduction, in the order this
// hands it. Sorting the nodes' own heights instead put a subtraction of two live
// scalars in every width, and the reduction then carried a weight derivative that
// the walk it stands in for structurally cannot have.
//
// Dropping the zero-density nodes instead would be wrong. A node whose density
// has collapsed to exactly zero contributes f = 0, and that zero is meaningful --
// it is the reconstruction saying density vanishes at that size. Removing those
// grid points would interpolate live density straight across the band and
// overestimate it, so they stay in and the grid gets sorted.
template <typename T, typename E>
std::vector<std::size_t> Species<T,E>::ascending_by_abscissa() const {
  std::vector<double> at(size());
  std::vector<std::size_t> order(size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    at[i] = abscissa_of(nodes[i], false);
    order[i] = i;
  }
  std::sort(order.begin(), order.end(),
            [&at](std::size_t a, std::size_t b) -> bool { return at[a] < at[b]; });
  return order;
}

// The reduction up to its closing trapezium. Everything the closing term needs
// travels in the result, so a caller that has to wait for the boundary node can
// close it later without walking the nodes again.
template <typename T, typename E>
typename Species<T,E>::competition_split
Species<T,E>::compute_competition_and_slope_split(const value_type& height) const {
  const HeightScan& scan = scan_heights();
  if (size() == 0 || scan.h_max < height) {
    return competition_split();
  }
  // Introduction times ascend by construction; -height ascends only while the
  // heights fall, and TF24's reserve-gated growth lets two cohorts cross.
  if (density_in_birth_date() || scan.decreasing) {
    return reduce_competition(height, {});
  }
  return reduce_competition(height, ascending_by_abscissa());
}

template <typename T, typename E>
with_slope<typename Species<T,E>::value_type>
Species<T,E>::close_competition_and_slope(const competition_split& c,
                                          const value_type& height) const {
  if (!c.closes) {
    return c.without_boundary();
  }
  const with_slope<value_type> fs0 =
    new_node.compute_competition_and_slope(height);
  const double x0 =
    abscissa_of(new_node, density_in_birth_date());
  return {(c.sum.value + (x0 - c.x1) * (c.at_x1.value + fs0.value)) / 2,
          (c.sum.slope + (x0 - c.x1) * (c.at_x1.slope + fs0.slope)) / 2};
}

template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment, double pr_patch_survival, const value_type& birth_rate) {
  for (auto& c : nodes) {
    c.compute_rates(environment, pr_patch_survival);
  }
  // The boundary condition, evaluated in the field the nodes above were just
  // rated in. This is not the evaluation the field itself reads -- that one is in
  // a field excluding the boundary interval, and Patch::compute_environment owns
  // it -- so the two are the same function at different arguments rather than one
  // computed twice. This value is the one an introduced node inherits.
  new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
  if (top_on()) {
    compute_top_rates(environment, pr_patch_survival);
  }
}

// The pinned individual's rates, height fixed; its density and fitness weights
// decay at its mortality rate, and the seed output accrues at its fecundity rate
// times the patch survival, as a node's does.
template <typename T, typename E>
void Species<T,E>::compute_top_rates(const E& environment,
                                     double pr_patch_survival) {
  if (!(top_density_ > 0.0) && !(top_weight_ > 0.0) && !(top_weight_birth_ > 0.0)) {
    top_.zero_rates();
    top_density_dt_ = top_weight_dt_ = top_weight_birth_dt_ = value_type(0.0);
    top_seeds_dt_ = top_seeds_birth_dt_ = value_type(0.0);
    return;
  }
  top_.compute_rates(environment);
  top_.hold_state(HEIGHT_INDEX);
  const value_type mortality = top_.rate(MORTALITY_INDEX);
  const value_type fecundity = top_.rate(FECUNDITY_INDEX);
  top_density_dt_ = -mortality * top_density_;
  top_weight_dt_ = -mortality * top_weight_;
  top_weight_birth_dt_ = -mortality * top_weight_birth_;
  top_seeds_dt_ = fecundity * pr_patch_survival * top_weight_;
  top_seeds_birth_dt_ = fecundity * pr_patch_survival * top_weight_birth_;
}

template <typename T, typename E>
double Species<T,E>::absorb_height() const {
  if constexpr (HasAsymptoticHeight<T>) {
    return strategy->asymptotic_height(control().absorb_epsilon);
  } else {
    util::stop("control$absorb_epsilon needs a strategy that declares an "
               "asymptotic height, and this one does not");
    return 0.0;
  }
}

template <typename T, typename E>
size_t Species<T,E>::absorb_tall(bool arrival_weighting) {
  if (!top_on()) {
    return 0;
  }
  flush_top();
  if (size() < 2) {
    return 0;
  }
  using odelia::util::to_passive;
  const double h_abs = absorb_height();
  const double S_D = to_passive(strategy->pars.S_D);
  size_t count = 0;
  // The youngest cohort is never absorbed: its quadrature weight is not final
  // until the next one is introduced.
  for (size_t k = 0; k + 1 < size(); ++k) {
    node_type& c = nodes[k];
    if (!(c.get_density() > 0.0) || to_passive(c.height()) < h_abs) {
      continue;
    }
    if (arrival_weighting && c.patch_density() > 0.0) {
      continue;
    }
    // The cohort's trapezium weight on the birth-date grid, retired neighbours
    // included: they are zero-valued grid points and keep their intervals.
    const double left = k > 0 ?
      nodes[k].introduction_time() - nodes[k - 1].introduction_time() : 0.0;
    const double right = nodes[k + 1].introduction_time() - nodes[k].introduction_time();
    const double weight = (left + right) / 2.0;
    const value_type added = weight * c.get_density();
    if (!(top_density_ > 0.0)) {
      top_ = c.individual;
    } else {
      // The states the layer integrates are the density-weighted mean of what
      // it holds and what it absorbs.
      for (size_t i = 0; i < individual_type::ode_size(); ++i) {
        if (static_cast<int>(i) == HEIGHT_INDEX) {
          continue;
        }
        top_.set_state(i, (top_density_ * top_.state(i) +
                           added * c.individual.state(i)) /
                          (top_density_ + added));
      }
    }
    top_.set_state(HEIGHT_INDEX, value_type(h_abs));
    top_density_ += added;
    const value_type w = weight * c.patch_density() * S_D *
      c.survival_individual() / c.get_pr_patch_survival_at_birth();
    top_weight_ += w;
    top_weight_birth_ +=
      w * strategy->birth_rate_at(c.introduction_time());
    c.retire();
    ++count;
  }
  if (count > 0) {
    invalidate_height_scan();
  }
  return count;
}

template <typename T, typename E>
Rcpp::List Species<T,E>::r_top_layer() const {
  using odelia::util::to_passive;
  return Rcpp::List::create(
      Rcpp::_["on"] = top_on(),
      Rcpp::_["height"] = top_on() && top_density_ > 0.0 ? top_height() : NA_REAL,
      Rcpp::_["density"] = to_passive(top_density_),
      Rcpp::_["weight"] = to_passive(top_weight_),
      Rcpp::_["weight_birth"] = to_passive(top_weight_birth_),
      Rcpp::_["seeds"] = to_passive(top_seeds_),
      Rcpp::_["seeds_birth"] = to_passive(top_seeds_birth_));
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node(double time, double patch_density,
                                      double pr_patch_survival) {
  invalidate_height_scan();
  // Stamp the pushed copy (not new_node) so the member stays pristine for
  // the no-arg introduction paths.
  nodes.push_back(new_node);
  nodes.back().set_birth_state(time, patch_density, pr_patch_survival);
}

template <typename T, typename E>
void Species<T,E>::remove_newest_node() {
  if (nodes.empty()) {
    util::stop("no node to remove from this species");
  }
  invalidate_height_scan();
  nodes.pop_back();
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(odelia::util::to_passive(c.fecundity()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node_weighted() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(
      odelia::util::to_passive(c.weighted_fecundity(strategy->pars.S_D)));
  }
  return ret;
}

template <typename T, typename E>
std::vector<typename Species<T,E>::value_type>
Species<T,E>::weighted_fecundity_by_node() const {
  std::vector<value_type> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.weighted_fecundity(strategy->pars.S_D));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::node_times() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.introduction_time());
  }
  return ret;
}

template <typename T, typename E>
bool Species<T,E>::birth_dates_are_distinct() const {
  for (size_t i = 1; i < nodes.size(); ++i) {
    if (nodes[i].introduction_time() == nodes[i - 1].introduction_time()) {
      return false;
    }
  }
  return true;
}

template <typename T, typename E>
void Species<T,E>::resize_consumption_rates(int r) {
  new_node.resize_consumption_rates(r);
  top_.resize_consumption_rates(r);
}

template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::consumption_rate(int i) const {
  if (size() == 0) {
    return value_type(0.0);
  }
  // The grid is abscissa_of, as the competition walk's is, so a width is a
  // position and not state on either coordinate and no reduction has a route to
  // an active one. Ascending in the abscissa is oldest first in birth date and
  // tallest first in height, and the closing node is both the newest and the
  // shortest, so it ends the grid either way.
  const bool birth_date = density_in_birth_date();
  std::vector<double> x;
  std::vector<value_type> rates;
  x.reserve(size() + 1);
  rates.reserve(size() + 1);
  for (auto& c : nodes) {
    x.push_back(abscissa_of(c, birth_date));
    rates.push_back(c.consumption_rate(i));
  }
  x.push_back(abscissa_of(new_node, birth_date));
  rates.push_back(new_node.consumption_rate(i));

  // Birth dates are strictly increasing by construction, so only the height
  // coordinate can arrive crossed -- and there neighbouring trapezia would cancel
  // instead of accumulating.
  if (!birth_date && !std::is_sorted(x.begin(), x.end())) {
    std::vector<size_t> order(x.size());
    for (size_t j = 0; j < order.size(); ++j) {
      order[j] = j;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) -> bool {
      return x[a] < x[b];
    });
    std::vector<double> x_sorted;
    std::vector<value_type> r_sorted;
    x_sorted.reserve(x.size());
    r_sorted.reserve(rates.size());
    for (size_t j : order) {
      x_sorted.push_back(x[j]);
      r_sorted.push_back(rates[j]);
    }
    x.swap(x_sorted);
    rates.swap(r_sorted);
  }
  value_type ret = util::trapezium(x, rates);
  if (top_on() && top_density_ > 0.0) {
    ret += top_density_ * top_.consumption_rate(i);
  }
  return ret;
}


template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::census_integral(const census_metric<T>& metric) const {
  if (size() == 0) {
    return value_type(0.0);
  }
  // The abscissa ascends as the nodes are walked in storage order and the
  // boundary node closes the grid: it is the youngest, and the shortest.
  const bool birth_date = density_in_birth_date();
  if (!birth_date && !heights_are_decreasing()) {
    util::stop("The census needs the node heights in decreasing order; on a"
               " crossed grid neighbouring trapezia cancel instead of"
               " accumulating");
  }
  std::vector<double> x;
  std::vector<value_type> weighted;
  x.reserve(size() + 1);
  weighted.reserve(size() + 1);
  for (auto& c : nodes) {
    x.push_back(abscissa_of(c, birth_date));
    weighted.push_back(c.get_density() * c.individual.census_value(metric));
  }
  x.push_back(abscissa_of(new_node, birth_date));
  weighted.push_back(new_node.get_density() *
                     new_node.individual.census_value(metric));
  value_type ret = util::trapezium(x, weighted);
  if (top_on() && top_density_ > 0.0) {
    ret += top_density_ * top_.census_value(metric);
  }
  return ret;
}

// bit clunky...
template <typename T, typename E>
size_t Species<T,E>::aux_size() const {
  return size() * strategy->aux_size();
}

template <typename T, typename E>
template <typename It>
It Species<T,E>::ode_aux(It it) const {
  return odelia::ode::ode_aux(nodes.begin(), nodes.end(), it);
}

template <typename T, typename E>
Rcpp::NumericMatrix Species<T, E>::r_get_state() const {

  size_t ode_size = node_type::ode_size(), n_nodes = size();
  size_t aux_size = strategy->aux_size();

  // On the birth-date path `log_density` is converted to the density in height
  // before it leaves C++, so every downstream consumer of this matrix keeps its
  // meaning, and the quantity actually integrated is reported alongside it as
  // `log_density_state`. Note export_patch_state() resumes from patch$ode_state,
  // not from here, so the raw state is what a resume reloads.
  const size_t extra = density_in_birth_date() ? 1 : 0;

  // Set output size. // +1 is seed
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size + aux_size + extra), n_nodes + 1);
  Rcpp::NumericMatrix::iterator it = ret.begin();

  for (size_t i = 0; i < n_nodes; ++i)
  {
    it = get_node_state(nodes[i], it);
    it = get_node_aux(nodes[i], it);
    it += extra;
  }

  it = get_node_state(new_node, it);
  it = get_node_aux(new_node, it);
  it += extra;

  // Combine ode_names and aux_names into a single vector for dimnames
  std::vector<std::string> names = node_type::ode_names();
  std::vector<std::string> aux = strategy->aux_names();
  names.insert(names.end(), aux.begin(), aux.end());

  if (extra > 0) {
    const int ld = static_cast<int>(
      std::find(names.begin(), names.end(), "log_density") - names.begin());
    const int st = static_cast<int>(ode_size + aux_size);
    names.push_back("log_density_state");
    // This matrix includes the boundary node as its last column, so the
    // Jacobian's trailing entry is used here (unlike r_log_densities()).
    const std::vector<double> jac = height_jacobian();
    for (int col = 0; col <= static_cast<int>(n_nodes); ++col) {
      const double nu = ret(ld, col);
      ret(st, col) = nu;
      ret(ld, col) = util::is_finite(jac[col]) ? nu - std::log(jac[col]) : NA_REAL;
    }
  }

  ret.attr("dimnames") = Rcpp::List::create(names, R_NilValue);
  // The top layer is no column of this matrix, so a reader integrating the
  // columns misses it; carried beside them so that reader can say so
  // (tidy_species). Only with the layer on, so any other state is unchanged.
  if (top_on()) {
    ret.attr("top_density") = odelia::util::to_passive(top_density_);
  }

  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->height()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights_rev() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->height()));
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void Species<T,E>::r_set_heights(std::vector<double> heights) {
  invalidate_height_scan();
  util::check_length(heights.size(), size());
  if (!util::is_decreasing(heights.begin(), heights.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (nodes_iterator it = nodes.begin(); it != nodes.end(); ++it, ++i) {
    it->individual.set_state("height", heights[i]);
  }
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_compute_competition_effect_by_nodes() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(odelia::util::to_passive(c.compute_competition(0.0)));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_compute_competition_effect_by_nodes_error(double scal) const {
  // Over the same abscissa the competition integral uses, so schedule
  // refinement measures the error of the quadrature actually being taken.
  // local_error_integration takes absolute differences, so the height branch's
  // sign flip leaves it unchanged.
  return util::local_error_integration(quadrature_abscissae(),
                                       r_compute_competition_effect_by_nodes(), scal);
}

// Central differences in the interior, one-sided at the two ends. The boundary
// node supplies the extra point at the young end, so the youngest real node --
// the one whose Jacobian a forward difference gets worst -- gets a genuine
// central difference. NA where the Jacobian vanishes: two cohorts at the same
// height is exactly where the height density is undefined (it is the multivalued
// case), and reporting +Inf there would be worse than reporting nothing.
template <typename T, typename E>
std::vector<double> Species<T,E>::height_jacobian() const {
  const size_t n = size();
  std::vector<double> h(n + 1), t(n + 1);
  for (size_t i = 0; i < n; ++i) {
    h[i] = nodes[i].height();
    t[i] = nodes[i].introduction_time();
  }
  h[n] = new_node.height();
  t[n] = new_node.introduction_time();

  std::vector<double> ret(n + 1, NA_REAL);

  // The boundary node needs no difference at all: dh/dtau = -g(H_0) at birth,
  // and compute_initial_conditions() re-evaluates that every step, so the value
  // it recorded is both exact and current. Only for *this* node -- an introduced
  // one has aged since, and its recorded rate is frozen at its own birth.
  const double g0 = new_node.growth_rate_at_birth();
  if (util::is_finite(g0) && g0 > 0.0) {
    ret[n] = g0;
  }

  if (n < 1) {
    return ret; // no interior to difference
  }
  for (size_t j = 0; j < n; ++j) {
    const size_t lo = (j == 0) ? 0 : j - 1;
    const size_t hi = j + 1;
    const double jac = std::abs((h[hi] - h[lo]) / (t[hi] - t[lo]));
    if (util::is_finite(jac) && jac > 0.0) {
      ret[j] = jac;
    }
  }
  if (!util::is_finite(ret[n])) {
    // Fall back to a one-sided difference if the birth rate is unavailable
    // (a node loaded from an exported state records zero).
    const double jac = std::abs((h[n] - h[n - 1]) / (t[n] - t[n - 1]));
    if (util::is_finite(jac) && jac > 0.0) {
      ret[n] = jac;
    }
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities_state() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->get_log_density()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities() const {
  std::vector<double> ret = r_log_densities_state();
  if (!density_in_birth_date()) {
    return ret;
  }
  // N = nu / |dh/dtau|. jac carries one extra trailing entry for the boundary
  // node, which this accessor does not report.
  const std::vector<double> jac = height_jacobian();
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i] = util::is_finite(jac[i]) ? ret[i] - std::log(jac[i]) : NA_REAL;
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_density_rates() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->get_log_density_rate()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_patch_densities() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(it->patch_density());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_pr_patch_survival_at_birth() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(it->get_pr_patch_survival_at_birth());
  }
  return ret;
}

template <typename T, typename E>
std::vector<bool> Species<T,E>::r_retired() const {
  std::vector<bool> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(it->retired());
  }
  return ret;
}

template <typename T, typename E>
void Species<T,E>::set_birth_state(const std::vector<double>& times,
                                   const std::vector<double>& patch_density,
                                   const std::vector<double>& pr_patch_survival) {
  util::check_length(times.size(), size());
  util::check_length(patch_density.size(), size());
  util::check_length(pr_patch_survival.size(), size());
  for (size_t i = 0; i < size(); ++i) {
    nodes[i].set_birth_state(times[i], patch_density[i], pr_patch_survival[i]);
  }
}

}

#endif /* SPECIES */
