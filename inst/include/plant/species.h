// -*-c++-*-
#ifndef SPECIES
#define SPECIES

#include <vector>
#include <algorithm>
#include <limits>
#include <utility>
#include <plant/util.h>
#include <plant/environment.h>
#include <odelia/ode_interface.hpp>
#include <plant/node.h>
#include <plant/species_base.h>
#include <odelia/drivers.hpp>

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
  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef Node<T,E> node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Species(strategy_type s);

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
  void introduce_new_node(double time, double patch_density);

  double height_max() const;
  double compute_competition(double height) const;
  // Whether the decreasing-height node ordering still holds (see height_max()).
  bool heights_are_decreasing() const;

  // The tallest node height and whether the heights are ordered, from a single
  // pass. compute_competition() needs both on every call, and walking the heights
  // twice was measurably slower on FF16 (~5% on the SCM benchmark) than the
  // O(1) nodes.front() it replaced.
  struct HeightScan { double h_max; bool decreasing; };
  // Cached: heights change only when the ODE state is set or a node is
  // introduced/cleared, whereas compute_competition() is called once per spline
  // knot, so this is hundreds of calls per change. Every mutator invalidates.
  HeightScan scan_heights() const;

  // Setting the ODE state rewrites every node's height, so the cached scan goes
  // with it. Shadows (rather than uses) the SpeciesBase version for that reason.
  odelia::ode::const_iterator set_ode_state(odelia::ode::const_iterator it) {
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
  double consumption_rate(int i) const;
  std::vector<double> consumption_rate_by_node_rev(int i) const;
  std::vector<double> consumption_rate_by_node(int i) const;

  odelia::ode::iterator       ode_aux(odelia::ode::iterator it) const;

  Rcpp::NumericMatrix r_get_state() const;

  // * R interface
  std::vector<double> r_heights() const;
  std::vector<double> r_heights_rev() const;
  void r_set_heights(std::vector<double> heights);
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

private:
  // compute_competition() for the case where the node heights are no longer
  // ordered, so the node list cannot be used directly as the quadrature grid.
  // Height coordinate only -- it integrates in height, and the birth-date
  // abscissa cannot invert (see compute_competition).
  double compute_competition_unordered(double height) const;

  // Cache for scan_heights(). Every path that can change a node height must call
  // invalidate_height_scan(); a stale cache here would silently reintroduce the
  // wrong competition profile of #571, so the coverage of these calls was checked
  // by asserting cache == freshly-computed on every call across the whole suite
  // and the scenario gateway.
  HeightScan compute_height_scan() const;
  void invalidate_height_scan() { height_scan_valid = false; }
  mutable HeightScan height_scan_cache{0.0, true};
  mutable bool height_scan_valid = false;

  // Storage (strategy, nodes) and control() live in SpeciesBase; the
  // using-declarations let the unqualified references below resolve through the
  // dependent base.
  using base_type::nodes;
  using base_type::strategy;
  using base_type::control;
  node_type new_node;

  // The abscissa both resource integrals are taken over, increasing as the node
  // list is walked from the tallest down. Heights are negated so that both
  // coordinates increase in the same direction; negation is exact, so the height
  // branch's trapezium widths are bit-identical to differencing the heights
  // themselves. Callers in hot loops read the coordinate once and pass it in.
  static double abscissa_of(const node_type& n, bool birth_date) {
    return birth_date ? n.introduction_time() : -n.height();
  }
  double quadrature_abscissa(const node_type& n) const {
    return abscissa_of(n, control().node_density_in_birth_date);
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
// This used to return nodes.front(), relying on the decreasing-height ordering
// asserted below. That ordering is guaranteed only while height growth is a
// function of height and the shared environment, which TF24 broke: its
// reserve-gated growth (#517) makes dh/dt depend on a cohort's own storage, so
// two cohorts born moments apart into a rapidly changing environment can cross
// in height. When they had, this returned a height 0.1 m *below* the tallest and
// only living cohort, truncating the light spline's domain (#571). Scanning is
// O(n) in heights only -- negligible against the crown integrals in
// compute_competition -- and returns exactly nodes.front() whenever the ordering
// does hold, so results are unchanged in that case.
template <typename T, typename E>
double Species<T,E>::height_max() const {
  if (nodes.empty()) {
    return new_node.height();
  }
  double ret = -std::numeric_limits<double>::infinity();
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret = std::max(ret, it->height());
  }
  return ret;
}

// Are the node heights still ordered largest to smallest? See height_max() above
// for why this can no longer be assumed. Heights only, so this is cheap relative
// to the per-node crown integrals it guards.
template <typename T, typename E>
bool Species<T,E>::heights_are_decreasing() const {
  return scan_heights().decreasing;
}

template <typename T, typename E>
typename Species<T,E>::HeightScan Species<T,E>::scan_heights() const {
  if (!height_scan_valid) {
    height_scan_cache = compute_height_scan();
    height_scan_valid = true;
  }
  return height_scan_cache;
}

// Tallest height and orderedness in one pass over the heights.
template <typename T, typename E>
typename Species<T,E>::HeightScan Species<T,E>::compute_height_scan() const {
  HeightScan ret{-std::numeric_limits<double>::infinity(), true};
  double h_prev = std::numeric_limits<double>::infinity();
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const double h = it->height();
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
double Species<T,E>::compute_competition(double height) const {
  if (size() == 0) {
    return 0.0;
  }
  const HeightScan scan = scan_heights();
  if (scan.h_max < height) {
    return 0.0;
  }
  // Read the coordinate once: this is the hottest loop in the solver (one pass
  // per spline knot per Runge-Kutta stage), so the control lookup does not
  // belong inside it.
  const bool birth_date = control().node_density_in_birth_date;
  // The loop below uses the node list itself as the quadrature grid, and the
  // early exit is valid only if that grid is monotone. When it is not, the exit
  // fires at the first node below `height` and silently drops every node beyond
  // it -- including, in #571, the only cohort with non-zero density, which put a
  // fictitious step in the competition profile. Take the ordered path instead.
  // Heights only, so the usual (ordered) case keeps this loop and its results
  // exactly.
  //
  // Only the height abscissa can invert. Introduction times are fixed at birth
  // and nodes are appended in that order, so the birth-date grid is monotone
  // whatever the heights do -- and compute_competition_unordered integrates in
  // height, so sending the birth-date coordinate down it would silently swap
  // coordinates mid-run.
  if (!birth_date && !scan.decreasing) {
    return compute_competition_unordered(height);
  }
  double tot = 0.0;
  nodes_const_iterator it = nodes.begin();
  double x1 = abscissa_of(*it, birth_date), f1 = it->compute_competition(height);

  // Loop over nodes
  for (++it; it != nodes.end(); ++it) {
    const double x0 = abscissa_of(*it, birth_date), h0 = it->height(),
                 f0 = it->compute_competition(height);
    if (!util::is_finite(f0)) {
      util::stop("Detected non-finite contribution");
    }
    // Integration
    tot += (x0 - x1) * (f1 + f0);
    // Upper point moves for next time:
    x1 = x0;
    f1 = f0;
    // It is the decreasing height ordering, not the abscissa, that licenses
    // stopping here: every later node is then shorter than `height` and
    // contributes nothing. On the birth-date axis that ordering can break while
    // the abscissa stays monotone, and a node below `height` may be followed by
    // a taller one, so walk the whole list instead. Always true on the height
    // path, which returned above otherwise.
    if (scan.decreasing && h0 < height) {
      break;
    }
  }

  // On the birth-date axis this segment is zero-width at the moment of
  // introduction and contributes nothing once the boundary node is below
  // `height`, so it is always safe to include; f1 can legitimately be zero here
  // when the walk ran to the end.
  if (size() == 1 || birth_date || f1 > 0) {
    const double x0 = abscissa_of(new_node, birth_date),
                 f0 = new_node.compute_competition(height);
    tot += (x0 - x1) * (f1 + f0);
  }

  return tot / 2;
}

// The same trapezium integral as compute_competition(), but over a height-sorted
// view of the nodes rather than the node list in place. Used only when the
// ordering has broken (#571): it agrees with the in-place version whenever the
// ordering holds, so this is a fallback rather than a change of method.
//
// Dropping the zero-density nodes instead would be wrong. A node whose density
// has collapsed to exactly zero contributes f = 0, and that zero is meaningful --
// it is the reconstruction saying density vanishes at that size. Removing those
// grid points would interpolate live density straight across the band and
// overestimate it, so they stay in and the grid gets sorted.
//
// No early exit here: it would need the same monotonicity that is missing. Nodes
// below `height` contribute f = 0 at both ends, so including them costs time but
// changes nothing. The scratch buffer is thread_local and reused, so the repeated
// calls that build one spline do not each allocate.
template <typename T, typename E>
double Species<T,E>::compute_competition_unordered(double height) const {
  thread_local std::vector<std::pair<double, double>> hf;
  hf.clear();
  hf.reserve(size());

  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const double f = it->compute_competition(height);
    if (!util::is_finite(f)) {
      util::stop("Detected non-finite contribution");
    }
    hf.push_back({it->height(), f});
  }
  std::sort(hf.begin(), hf.end(),
            [](std::pair<double, double> const& a,
               std::pair<double, double> const& b) {
              return a.first > b.first;
            });

  double tot = 0.0;
  double h1 = hf.front().first, f_h1 = hf.front().second;
  for (size_t j = 1; j < hf.size(); ++j) {
    const double h0 = hf[j].first, f_h0 = hf[j].second;
    tot += (h1 - h0) * (f_h1 + f_h0);
    h1   = h0;
    f_h1 = f_h0;
  }

  if (size() == 1 || f_h1 > 0) {
    const double h0 = new_node.height(), f_h0 = new_node.compute_competition(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
  }

  return tot / 2;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment, double pr_patch_survival, double birth_rate) {
  for (auto& c : nodes) {
    c.compute_rates(environment, pr_patch_survival);
  }
  new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node(double time, double patch_density) {
  invalidate_height_scan();
  // Stamp the pushed copy (not new_node) so the member stays pristine for
  // the no-arg introduction paths.
  nodes.push_back(new_node);
  nodes.back().set_introduction(time, patch_density);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.fecundity());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node_weighted() const {
  std::vector<double> ret;
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
}

template <typename T, typename E>
double Species<T,E>::consumption_rate(int i) const {
  if (size() == 0) {
    return 0.0;
  }
  if (control().node_density_in_birth_date) {
    // Introduction times are fixed at birth and nodes are appended in that
    // order, so this grid is ascending however the heights behave -- there is no
    // inverted case to sort. new_node's birth date is the current time, which is
    // the newest, so it goes on the end rather than the front.
    std::vector<double> times = node_times();
    times.push_back(new_node.introduction_time());
    std::vector<double> rates = consumption_rate_by_node(i);
    rates.push_back(new_node.consumption_rate(i));
    return util::trapezium(times, rates);
  }
  // node heights are in descending order - we need ascending for integration,
  // starting at new_node, which is where the size distribution starts.
  std::vector<double> heights = r_heights_rev();
  heights.insert(heights.begin(), new_node.height());
  std::vector<double> rates = consumption_rate_by_node_rev(i);

  // The node list is the quadrature grid here as it is in compute_competition,
  // so an inverted grid (#571) makes neighbouring trapezia cancel rather than
  // accumulate. Sort the pairs when the ordering has broken, as
  // compute_competition_unordered does; an already-ascending grid is untouched.
  if (!std::is_sorted(heights.begin(), heights.end())) {
    std::vector<std::pair<double, double>> hr;
    hr.reserve(heights.size());
    for (size_t j = 0; j < heights.size(); ++j) {
      hr.push_back({heights[j], rates[j]});
    }
    std::sort(hr.begin(), hr.end(),
              [](std::pair<double, double> const& a,
                 std::pair<double, double> const& b) {
                return a.first < b.first;
              });
    for (size_t j = 0; j < hr.size(); ++j) {
      heights[j] = hr[j].first;
      rates[j] = hr[j].second;
    }
  }
  return util::trapezium(heights, rates);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::consumption_rate_by_node_rev(int i) const {
  std::vector<double> ret;
  ret.reserve(size() + 1);
  ret.push_back(new_node.consumption_rate(i));
  for(auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
    ret.push_back(it->consumption_rate(i));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::consumption_rate_by_node(int i) const {
  std::vector<double> ret;
  ret.reserve(size());
  for(auto& c : nodes) {
    ret.push_back(c.consumption_rate(i));
  }
  return ret;
}

// bit clunky...
template <typename T, typename E>
size_t Species<T,E>::aux_size() const {
  return size() * strategy->aux_size();
}

template <typename T, typename E>
odelia::ode::iterator Species<T,E>::ode_aux(odelia::ode::iterator it) const {
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
    ret.push_back(it->height());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights_rev() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(it->height());
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
    ret.push_back(c.compute_competition(0.0));
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
    ret.push_back(it->get_log_density());
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
    ret.push_back(it->get_log_density_rate());
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
