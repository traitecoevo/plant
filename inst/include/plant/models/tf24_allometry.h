// -*-c++-*-
// TF24_Model's allometry: the sizes and masses a plant of a given height
// and leaf area has, their derivatives, its seed, and how its roots are
// spread over soil depth.
// The class is declared in tf24_model.h, which includes this file.
#ifndef PLANT_PLANT_TF24_ALLOMETRY_H_
#define PLANT_PLANT_TF24_ALLOMETRY_H_

#include <plant/models/tf24_model.h>

namespace plant {

// Leaf area from height
template <typename S>
S TF24_Model<S>::area_leaf(const S& height) const {
  return pow(height / pars.a_l1, 1.0 / pars.a_l2);
}

// Leaf mass from leaf area
template <typename S>
S TF24_Model<S>::mass_leaf(const S& area_leaf) const {
  return area_leaf * pars.lma;
}

// Basal sapwood area: the conducting area a plant actually has (at a resting
// departure, the one its leaf area prefers), read at the base of the flow path.
// On a uniform stem at a resting departure the Huber value is theta_tip * 1 * 1,
// both multiplications exact, so this is area_leaf * theta_tip bit for bit.
template <typename S>
S TF24_Model<S>::area_sapwood(const S& area_leaf, const S& sapwood_departure,
                              const S& height) const {
  return area_leaf * sapwood_per_leaf_area_at_base(sapwood_departure, height);
}

template <typename S>
S TF24_Model<S>::mass_sapwood(const S& area_sapwood, const S& height) const {
  return stem_tissue_mass(area_sapwood, height);
}

// Basal bark area. The trailing factor is exactly 1 on a uniform stem.
template <typename S>
S TF24_Model<S>::area_bark(const S& area_leaf, const S& height) const {
  return pars.a_b1 * area_leaf * pars.theta_tip * huber_base_factor(height);
}

template <typename S>
S TF24_Model<S>::mass_bark(const S& area_bark, const S& height) const {
  return stem_tissue_mass(area_bark, height);
}

// A uniform stem's tissue is its basal area over the mean path length
// height * eta_c, in that order of operations. A profiled one is its basal area
// times V / theta_base-factor: the volume per unit tip area, over the basal
// area per unit tip area. The attached column is d(V/m)/dtheta_c at 0, which
// is dV/dtheta_c - L_top * log(L_top/L_tip) = L_tip - L_top above the segment.
template <typename S>
S TF24_Model<S>::stem_tissue_mass(const S& basal_area, const S& height) const {
  if (uniform_stem()) {
    S m = basal_area * height * eta_c * pars.rho;
    if constexpr (!std::is_same_v<S, double>) {
      const S L_top = height * eta_c;
      if (path_beyond_tip(L_top)) {
        m += basal_area * pars.theta_c * (pars.L_tip - L_top) * pars.rho;
      }
    }
    return m;
  }
  return basal_area *
    (stem_volume_per_tip_area(height) / huber_base_factor(height)) * pars.rho;
}

template <typename S>
S TF24_Model<S>::area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <typename S>
S TF24_Model<S>::diameter_stem(const S& area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// Mass of (fine) roots
template <typename S>
S TF24_Model<S>::mass_root(const S& area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// Mass of coarse (structural) roots
template <typename S>
S TF24_Model<S>::mass_coarse_root(const S& mass_sapwood) const {
  return pars.a_cr1 * mass_sapwood;
}

// Total mass
// The coarse-root term is appended, not interleaved, so that with the default
// a_cr1 = 0 the sum is bit-identical to the pre-#349 expression (x + 0.0 == x
// exactly, and the association order of the other terms is untouched).
template <typename S>
S TF24_Model<S>::mass_live(const S& mass_leaf, const S& mass_bark,
                           const S& mass_sapwood, const S& mass_root,
                           const S& mass_coarse_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root + mass_coarse_root;
}

template <typename S>
S TF24_Model<S>::mass_total(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood,
                            const S& mass_root,
                            const S& mass_coarse_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root +
    mass_coarse_root;
}

template <typename S>
S TF24_Model<S>::mass_above_ground(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// for updating auxiliary state
template <typename S>
void TF24_Model<S>::update_dependent_aux(const int index, Internals<S>& vars) {
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
S TF24_Model<S>::darea_leaf_dmass_live(const S& area_leaf,
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
S TF24_Model<S>::dheight_darea_leaf(const S& area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_l / d a_l
template <typename S>
S TF24_Model<S>::dmass_leaf_darea_leaf(const S& /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Model<S>::dmass_sapwood_darea_leaf(const S& area_leaf,
                                             const S& sapwood_departure) const {
  // Priced at the Huber value the plant ACTUALLY carries. At rest
  // sapwood_per_leaf_area() returns pars.theta_tip unchanged, so on a uniform
  // stem this is bitwise the expression it replaces.
  //
  // Derived by hand rather than by differentiating mass_sapwood, so it must be
  // kept in step with it by hand: test-tf24-stem-hydraulics.R checks it against
  // a tangent derivative of the mass. The mass is
  //   M_s(A) = rho * theta_tip * exp(psi) * A * V(H(A)),   H(A) = a_l1 * A^a_l2,
  // along the resting allometry through A, and since A * dH/dA = a_l2 * H and
  // dV/dH = eta_c * theta_base/theta_tip,
  //   dM_s/dA = rho * theta_tip * exp(psi) * (V + a_l2 * H * eta_c * m(H)).
  // On a uniform stem V = H * eta_c and m = 1, which is the line below.
  if (uniform_stem()) {
    S d = pars.rho * eta_c * pars.a_l1 * sapwood_per_leaf_area(sapwood_departure) *
      (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
    if constexpr (!std::is_same_v<S, double>) {
      // d/dtheta_c at 0 of V + a_l2 * H * eta_c * m: dV/dtheta_c plus
      // a_l2 * H * eta_c * log(L_top/L_tip). Both vanish below the segment.
      const S height = pars.a_l1 * pow(area_leaf, pars.a_l2);
      const S L_top = height * eta_c;
      if (path_beyond_tip(L_top)) {
        d += pars.rho * sapwood_per_leaf_area(sapwood_departure) * pars.theta_c *
          (dstem_volume_dtheta_c_uniform(height) +
           pars.a_l2 * L_top * log(L_top / pars.L_tip));
      }
    }
    return d;
  }
  const S height = pars.a_l1 * pow(area_leaf, pars.a_l2);
  return pars.rho * sapwood_per_leaf_area(sapwood_departure) *
    (stem_volume_per_tip_area(height) +
     pars.a_l2 * height * eta_c * huber_base_factor(height));
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <typename S>
S TF24_Model<S>::dmass_bark_darea_leaf(const S& area_leaf) const {
  // At a RESTING sapwood departure deliberately: bark follows leaf area, so it
  // must not pick up the Huber value's excursion. Expressed as a_b1 times the
  // sapwood form at zero departure rather than written out, because floating
  // point multiplication is not associative -- writing the product flat changes
  // the last bits and moved a whole stand trajectory when first tried.
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf, S(0.0));
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <typename S>
S TF24_Model<S>::dmass_root_darea_leaf(const S& /* area_leaf */) const {
  return pars.a_r1;
}

// Mass of coarse root needed for new unit area leaf, d m_cr / d a_l
template <typename S>
S TF24_Model<S>::dmass_coarse_root_darea_leaf(const S& area_leaf,
                                                 const S& sapwood_departure) const {
  return pars.a_cr1 * dmass_sapwood_darea_leaf(area_leaf, sapwood_departure);
}

// Heartwood grows from the sapwood already there, not from leaf-area growth,
// and from the basal area, which is the stem it thickens.
template <typename S>
S TF24_Model<S>::area_heartwood_dt(const S& area_leaf,
                                   const S& sapwood_departure,
                                   const S& height) const {
  return pars.k_s * area_sapwood(area_leaf, sapwood_departure, height);
}

template <typename S>
S TF24_Model<S>::mass_heartwood_dt(const S& mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <typename S>
S TF24_Model<S>::mass_live_given_height(const S& height) const {
  S area_leaf_ = area_leaf(height);
  const S mass_sapwood_ =
    mass_sapwood(area_sapwood(area_leaf_, S(0.0), height), height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_, height), height) +
         mass_sapwood_ +
         mass_root(area_leaf_) +
         mass_coarse_root(mass_sapwood_);
}

template <typename S>
S TF24_Model<S>::height_given_mass_leaf(const S& mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

// Cumulative fraction of a quantity distributed over an extent with
// shape exponent 'eta_x', above coordinate 'z' of a total 'height'.
// Serves the root mass distribution over soil depth.
template <typename S>
S TF24_Model<S>::Q(const S& z, const S& rooting_depth, const S& eta_x) const {
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

// Roots go as deep as the plant is tall, up to the soil's depth x_max. The cap
// is a smooth minimum,
//
//   D = x_max - w * softplus((x_max - h) / w),   w = 0.05 * x_max,
//
// which departs from min(h, x_max) only within a few w of the cap: by
// w * log(2) at the cap, and by w * exp(-|h - x_max| / w) away from it. Its
// slope in h is the logistic sigmoid((x_max - h) / w), so it has no kink and
// the root profile keeps responding to height above the cap, ever more weakly.
//
// D is increasing in h, and below the cap D = h - w * log1p(exp(-(x_max - h) / w)),
// which is positive for every h above about 2e-9 * w: far below any seed's height.
//
// softplus(y) = log(1 + exp(y)) is evaluated in the form that neither
// overflows for large y nor loses its tail for very negative y. The two
// branches agree in value and slope at y = 0, so the branch taken is no kink.
template <typename S>
S TF24_Model<S>::rooting_depth_at(const S& height) const {
  using std::exp;
  using std::log1p;
  const S w = 0.05 * pars.rooting_depth_max;
  const S y = (pars.rooting_depth_max - height) / w;
  const S softplus = odelia::util::to_passive(y) > 0.0
    ? S(y + log1p(exp(-y)))
    : S(log1p(exp(y)));
  return pars.rooting_depth_max - w * softplus;
}

// The aim is to find a plant height that gives the correct seed mass.
template <typename S>
double TF24_Model<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const S
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(pars.omega);

  const double tol = this->control.offspring_height_tol;

  auto target = [&] (double x) mutable -> S {
    return mass_live_given_height(x) - pars.omega;
  };

  if constexpr (std::is_same_v<S, double>) {
    return odelia::root::toms748(target, h0, h1, {tol, tol});
  } else {
    // Recording the search would differentiate the iterate, whose path through
    // the bracket is not the root's, rather than the height at which mass_live
    // equals the seed mass. Declare it by the residual instead.
    static_assert(std::is_same_v<S, double>,
                  "height_seed() finds its root by iteration; an active scalar must "
                  "declare it through implicit_value on the residual.");
  }
}

}

#endif
