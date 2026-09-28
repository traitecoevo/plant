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
  using base_type::ode_size;
  using base_type::ode_state;
  using base_type::ode_rates;
  using base_type::get_node_state;
  using base_type::get_node_aux;
  typename std::vector<node_type>::iterator node_begin() { return nodes.begin(); }
  typename std::vector<node_type>::iterator node_end() { return nodes.end(); }
  typename std::vector<node_type>::const_iterator node_begin() const { return nodes.begin(); }
  typename std::vector<node_type>::const_iterator node_end() const { return nodes.end(); }

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
                             double pr_patch_survival, double birth_rate) {
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
                              strategy);
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
    return base_type::set_ode_state(it);
  }
  void compute_rates(const environment_type& environment, double pr_patch_survival, double birth_rate);

  std::vector<double> net_reproduction_ratio_by_node() const;
  // Per-node lifetime offspring, weighted by patch-age density and S_D.
  std::vector<double> net_reproduction_ratio_by_node_weighted() const;
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
  // seeded or imported without per-node times.
  bool birth_dates_are_distinct() const;

  // Which coordinate this species' size distribution is carried in. Exposed so
  // Patch can check every species agrees before summing their contributions.
  bool density_in_birth_date() const {
    return control().node_density_in_birth_date;
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

  ExtrinsicDrivers extrinsic_drivers() const {return strategy->extrinsic_drivers;}

  // The prepared strategy this species and its nodes share.
  strategy_type_ptr strategy_ptr() const {return this->strategy;}

private:
  // The reduction, over the nodes in ascending abscissa. `order` names that order
  // where the node list is not already in it and is empty where it is; the early
  // exit is the decreasing heights', not this parameter's.
  competition_split reduce_competition(const value_type& height,
                                       const std::vector<std::size_t>& order) const;
  // The node positions in ascending abscissa, for the case where the heights are
  // no longer ordered and the node list cannot be the quadrature grid.
  // Positions only: nothing here evaluates a contribution, so no width the
  // reduction forms out of it can carry a derivative.
  std::vector<std::size_t> ascending_by_abscissa() const;
  // Whether the closing trapezium to the boundary node is taken. Below a node
  // that no longer contributes the density is zero and so is the interval; a
  // single node has no interval but that one, and a birth-date abscissa does not
  // order with the support at all.
  bool closes_on(const value_type& value_at_x1) const {
    return size() == 1 || control().node_density_in_birth_date ||
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
    const bool birth_date = control().node_density_in_birth_date;
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
    new_node(this->strategy) {
}

template <typename T, typename E>
Species<T,E>::Species(strategy_type_ptr s)
  : base_type(s),
    new_node(this->strategy) {
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
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node() {
  invalidate_height_scan();
  // new_node already holds the initial conditions computed against the current
  // environment by the most recent compute_rates() call (see compute_rates ->
  // new_node.compute_initial_conditions above), and the member is refreshed
  // again on the next compute_rates() ready for the following introduction.
  // Recomputing it here would be redundant, and would (wrongly) re-seed against
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
  if (nodes.empty()) {
    return new_node.height();
  }
  return scan_heights().h_max;
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
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const value_type h = it->height();
    if (h > h_prev) {
      ret.decreasing = false;
    }
    if (h > ret.h_max) {
      ret.h_max = h;
    }
    h_prev = h;
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
  return close_competition_and_slope(compute_competition_and_slope_split(height),
                                     height);
}

// One pass over the nodes for the whole height set. Where the prefix form does not
// hold -- a broken ordering, or a profile that is not a polynomial in w -- every
// height takes the walk, which is the same branch the walk's own early exit rests
// on.
template <typename T, typename E>
void Species<T,E>::field_splits(const std::vector<value_type>& heights,
                                std::vector<competition_split>& out) const {
  out.assign(heights.size(), competition_split());
  if (size() == 0) {
    return;
  }
  const HeightScan& scan = scan_heights();
  const bool birth_date = control().node_density_in_birth_date;
  const std::size_t n_moments = strategy->canopy_shape.n_moments();
  if (!scan.decreasing || n_moments == 0) {
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
  std::size_t crossing = n;
  moments weight, weight_slope;
  for (std::size_t k = 0; k < heights.size(); ++k) {
    const value_type& height = heights[k];
    if (scan.h_max < height) {
      continue;  // no node reaches it; the empty split stands
    }
    while (crossing > 0 && nodes[crossing - 1].height() < height) {
      --crossing;
    }
    // The walk stops AFTER the interval whose upper node is the first below the
    // height, so the last node it visited is that one -- or the last node of all,
    // where none is below.
    const std::size_t last = crossing < n ? crossing : n - 1;
    const std::size_t summed = last > 0 ? last - (crossing < n ? 1 : 0) : 0;

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
    if (crossing < n) {
      const with_slope<value_type> fs_above =
        nodes[crossing - 1].compute_competition_and_slope(height);
      const double width = abscissa[crossing] - abscissa[crossing - 1];
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
  const bool birth_date = control().node_density_in_birth_date;
  const HeightScan& scan = scan_heights();
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
    if (scan.decreasing && node.height() < height) {
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
  const bool birth_date = control().node_density_in_birth_date;
  std::vector<double> at(size());
  std::vector<std::size_t> order(size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    at[i] = abscissa_of(nodes[i], birth_date);
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
  if (control().node_density_in_birth_date || scan.decreasing) {
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
    abscissa_of(new_node, control().node_density_in_birth_date);
  return {(c.sum.value + (x0 - c.x1) * (c.at_x1.value + fs0.value)) / 2,
          (c.sum.slope + (x0 - c.x1) * (c.at_x1.slope + fs0.slope)) / 2};
}

template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment, double pr_patch_survival, double birth_rate) {
  for (auto& c : nodes) {
    c.compute_rates(environment, pr_patch_survival);
  }
  // The boundary condition, evaluated in the field the nodes above were just
  // rated in. This is not the evaluation the field itself reads -- that one is in
  // a field excluding the boundary interval, and Patch::compute_environment owns
  // it -- so the two are the same function at different arguments rather than one
  // computed twice. This value is the one an introduced node inherits.
  new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
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
  const bool birth_date = control().node_density_in_birth_date;
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
  return util::trapezium(x, rates);
}


template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::census_integral(const census_metric<T>& metric) const {
  if (size() == 0) {
    return value_type(0.0);
  }
  // The abscissa ascends as the nodes are walked in storage order and the
  // boundary node closes the grid: it is the youngest, and the shortest.
  const bool birth_date = control().node_density_in_birth_date;
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
  return util::trapezium(x, weighted);
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
  const size_t extra = control().node_density_in_birth_date ? 1 : 0;

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
  if (!control().node_density_in_birth_date) {
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
