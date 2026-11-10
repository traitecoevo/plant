// -*-c++-*-
// TF24_Strategy's allometry: the sizes and masses a plant of a given height
// and leaf area has, their derivatives, its seed, and how its roots are
// spread over soil depth.
// The class is declared in tf24_model.h, which includes this file.
#ifndef PLANT_PLANT_TF24_ALLOMETRY_H_
#define PLANT_PLANT_TF24_ALLOMETRY_H_

#include <plant/models/tf24_model.h>

namespace plant {

// area_leaf
template <typename S>
S TF24_Strategy<S>::area_leaf(const S& height) const {
  return pow(height / pars.a_l1, 1.0 / pars.a_l2);
}

// mass_leaf
template <typename S>
S TF24_Strategy<S>::mass_leaf(const S& area_leaf) const {
  return area_leaf * pars.lma;
}

// Area and mass of sapwood
template <typename S>
S TF24_Strategy<S>::area_sapwood(const S& area_leaf) const {
  return area_leaf * pars.theta;
}

// The conducting area a plant actually has, as against the one its leaf area
// prefers. At a resting departure sapwood_per_leaf_area() returns pars.theta
// unchanged, so this is bit-identical to the overload above.
template <typename S>
S TF24_Strategy<S>::area_sapwood(const S& area_leaf,
                                 const S& sapwood_departure) const {
  return area_leaf * sapwood_per_leaf_area(sapwood_departure);
}

template <typename S>
S TF24_Strategy<S>::mass_sapwood(const S& area_sapwood, const S& height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// Area and mass of bark
template <typename S>
S TF24_Strategy<S>::area_bark(const S& area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <typename S>
S TF24_Strategy<S>::mass_bark(const S& area_bark, const S& height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <typename S>
S TF24_Strategy<S>::area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <typename S>
S TF24_Strategy<S>::diameter_stem(const S& area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// Mass of (fine) roots
template <typename S>
S TF24_Strategy<S>::mass_root(const S& area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// Mass of coarse (structural) roots
template <typename S>
S TF24_Strategy<S>::mass_coarse_root(const S& mass_sapwood) const {
  return pars.a_cr1 * mass_sapwood;
}

// Total mass
// The coarse-root term is appended, not interleaved, so that with the default
// a_cr1 = 0 the sum is bit-identical to the pre-#349 expression (x + 0.0 == x
// exactly, and the association order of the other terms is untouched).
template <typename S>
S TF24_Strategy<S>::mass_live(const S& mass_leaf, const S& mass_bark,
                           const S& mass_sapwood, const S& mass_root,
                           const S& mass_coarse_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root + mass_coarse_root;
}

template <typename S>
S TF24_Strategy<S>::mass_total(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood,
                            const S& mass_root,
                            const S& mass_coarse_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root +
    mass_coarse_root;
}

template <typename S>
S TF24_Strategy<S>::mass_above_ground(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// for updating auxiliary state
template <typename S>
void TF24_Strategy<S>::update_dependent_aux(const int index, Internals<S>& vars) {
  // Leaf area depends on BOTH height and the departure state, so either write
  // has to refresh it. Individual::set_ode_state refreshes per slot in slot
  // order, so the height pass (slot 0) necessarily reads the departure still
  // sitting in vars from the previous evaluation; the departure pass recomputes
  // from both, which is why it must not be a partial update. Whichever fires
  // last leaves the aux consistent with the whole state vector.
  if (index == HEIGHT_INDEX || index == state_idx_log_area_leaf_departure) {
    const S& height = vars.state(HEIGHT_INDEX);
    // A = A*(h) * exp(phi). With plasticity off phi is exactly 0 and exp(0) is
    // exactly 1, so this is bit-identical to the plain area_leaf(height) it
    // replaces -- the multiplication is exact, and there is no addition for an
    // FMA to fuse.
    vars.set_aux(aux_idx_competition_effect, area_leaf_carried(vars));
    vars.set_aux(aux_idx_height_inverse, 1.0 / height);
  }
}

template <typename S>
S TF24_Strategy<S>::darea_leaf_dmass_live(const S& area_leaf,
                                          const S& sapwood_departure) const {
  // The coarse-root term is appended last so that at the default a_cr1 = 0 the
  // denominator is bit-identical to the pre-#349 sum.
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf, sapwood_departure)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf)
              + dmass_coarse_root_darea_leaf(area_leaf, sapwood_departure));
}

template <typename S>
S TF24_Strategy<S>::dheight_darea_leaf(const S& area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_leaf_darea_leaf(const S& /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_sapwood_darea_leaf(const S& area_leaf,
                                             const S& sapwood_departure) const {
  // Priced at the Huber value the plant ACTUALLY carries. At rest
  // sapwood_per_leaf_area() returns pars.theta unchanged, so this is bitwise
  // the expression it replaces.
  return pars.rho * eta_c * pars.a_l1 * sapwood_per_leaf_area(sapwood_departure) *
    (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_bark_darea_leaf(const S& area_leaf) const {
  // At a RESTING sapwood departure deliberately: bark follows leaf area, so it
  // must not pick up the Huber value's excursion. Expressed as a_b1 times the
  // sapwood form at zero departure rather than written out, because floating
  // point multiplication is not associative -- writing the product flat changes
  // the last bits and moved a whole SCM trajectory when first tried.
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf, S(0.0));
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_root_darea_leaf(const S& /* area_leaf */) const {
  return pars.a_r1;
}

// Mass of coarse root needed for new unit area leaf, d m_cr / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_coarse_root_darea_leaf(const S& area_leaf,
                                                 const S& sapwood_departure) const {
  return pars.a_cr1 * dmass_sapwood_darea_leaf(area_leaf, sapwood_departure);
}

// Heartwood grows from the sapwood already there, not from leaf-area growth.
template <typename S>
S TF24_Strategy<S>::area_heartwood_dt(const S& area_leaf,
                                      const S& sapwood_departure) const {
  return pars.k_s * area_sapwood(area_leaf, sapwood_departure);
}

template <typename S>
S TF24_Strategy<S>::mass_heartwood_dt(const S& mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <typename S>
S TF24_Strategy<S>::mass_live_given_height(const S& height) const {
  S area_leaf_ = area_leaf(height);
  const S mass_sapwood_ = mass_sapwood(area_sapwood(area_leaf_), height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood_ +
         mass_root(area_leaf_) +
         mass_coarse_root(mass_sapwood_);
}

template <typename S>
S TF24_Strategy<S>::height_given_mass_leaf(const S& mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

// Cumulative fraction of a quantity distributed over an extent with
//          shape exponent 'eta_x', above coordinate 'z' of a total 'height'.
//          Serves the root mass distribution over soil depth.
template <typename S>
S TF24_Strategy<S>::Q(const S& z, const S& rooting_depth, const S& eta_x) const {
  if (z > rooting_depth) {
    return S(0.0);
  }
  // u^eta_x. On double the plain pow; on an active scalar the recorded eta_x
  // derivative u^eta_x * log(u) is 0 * (-inf) -- a NaN -- at u = 0, where the
  // cumulative fraction is 1, so the guard supplies that value outright.
  S u_eta;
  if constexpr (std::is_same_v<S, double>) {
    u_eta = pow(z / rooting_depth, eta_x);
  } else {
    const S u = z / rooting_depth;
    u_eta = odelia::util::to_passive(u) <= 0.0 ? S(0.0) : pow(u, eta_x);
  }
  const S tmp = 1.0 - u_eta;
  return tmp * tmp;
}

// The aim is to find a plant height that gives the correct seed mass.
template <typename S>
double TF24_Strategy<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const S
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(pars.omega);

  const double tol = this->control.offspring_production_tol;
  const size_t max_iterations = this->control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> S {
    return mass_live_given_height(x) - pars.omega;
  };

  if constexpr (std::is_same_v<S, double>) {
    return util::uniroot(target, h0, h1, tol, max_iterations);
  } else {
    // Bisection is affine in its bracket and blind to the residual's values, so
    // recording the search would return d(h1)/d(trait) rather than the height at
    // which mass_live equals the seed mass. Declare it by the residual instead.
    static_assert(std::is_same_v<S, double>,
                  "height_seed() finds its root by iteration; an active scalar must "
                  "declare it through implicit_value on the residual.");
  }
}

}

#endif
