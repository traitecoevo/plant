// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// plant's leaf gas-exchange and hydraulics model is phylloptim's, a header-only
// package reached through LinkingTo. This header aliases it into namespace plant,
// so plant's sources and the generated RcppR6/RcppExports glue, which name
// `plant::Leaf` throughout, compile against it.
//
// LinkingTo is not transitive in R: plant names BH and odelia in its own
// DESCRIPTION although it is phylloptim that includes them.

#include <phylloptim.hpp>

namespace plant {

// The leaf model itself.
using Leaf = ::phylloptim::Leaf;

// The leaf R's `Leaf()` builds: phylloptim's positional constructor, seated as
// this package's R code reads a leaf -- the layered supply, which also seats the
// collar route, and TF24. From phylloptim 0.29 a bare Leaf is TF24_floor on the
// leaf route over a single potential, which refuses the root network R hands it.
inline Leaf make_leaf(double vcmax_25, double stem_c, double stem_P50,
                      double root_c, double root_P50, double TF24_beta2,
                      double jmax_25, double phi_J, double curv_fact_elec_trans,
                      double curv_fact_colim, double GSS_tol_abs,
                      double vulnerability_curve_ncontrol, double ci_abs_tol,
                      double ci_niter, double TF24_cost_scale) {
  Leaf l(vcmax_25, stem_c, stem_P50, root_c, root_P50, TF24_beta2, jmax_25,
         phi_J, curv_fact_elec_trans, curv_fact_colim, GSS_tol_abs,
         vulnerability_curve_ncontrol, ci_abs_tol, ci_niter, TF24_cost_scale);
  l.set_supply_multilayer();
  l.set_model(Leaf::CostCurve::TF24, true);
  return l;
}

// The one leaf-model constant plant's own sources read: the molar mass of water,
// for TF24's soil-water consumption rate.
using ::phylloptim::kg_per_mol_h2o;

// The leaf converts a mole fraction to a partial pressure with the atmospheric
// pressure it is given, so TF24_Environment's `atm_kpa` driver (101.3 by
// default) reaches Gamma*, Kc, Ko, Km and the ci root-find's lower bound.
// Setting it for altitude moves all of those.

} // namespace plant

#endif
