// -*-c++-*-
#ifndef PLANT_PLANT_TF24_PRODUCTION_KERNEL_H_
#define PLANT_PLANT_TF24_PRODUCTION_KERNEL_H_

#include <cmath>   // std::pow / std::exp; XAD provides these for active types via ADL

// Scalar-templated core of the TF24 net-mass-production + demographic rate chain
// (#472 scope B / traitecoevo/plant#537, Phase F1-full). These pieces are the
// TF24 analogues of ff16_production_kernel.h. The MASS CASCADE, respiration,
// turnover and the whole demographic rate fill (growth / fecundity / heartwood /
// mortality downstream of net production) are ALGEBRAICALLY IDENTICAL to FF16;
// only the assimilation term differs -- TF24 forms it from the OPTIMISED leaf
// profit (assim = profit * area_leaf * conv) rather than the FF16 light-response
// hyperbola. They are the SINGLE SOURCE OF TRUTH: TF24_Strategy's double rate
// methods delegate to them (so the existing TF24 reference test validates
// faithfulness, bit-identical), and the AD calibration path instantiates them
// with an active scalar so net / the demographic rates differentiate w.r.t.
// traits by forward- or reverse-mode AD with no special handling.
//
// The leaf optimisation itself (profit*, a max over collar potential nesting a
// psi_stem->ci root-find) is NOT in this header: it stays in the leaf submodel,
// and its trait sensitivities enter the AD path as the Leaf::dprofit_d* numbers
// injected first-order into an active `profit` (the #539 IFT / FF16 height_0
// injection pattern). This kernel takes `profit` (or `assimilation`) as a given
// and carries it -- and every mass-cascade trait -- through to the five ODE rates.

namespace plant {

// Unit conversion folding leaf-area assimilation (umol CO2 m^-2 s^-1) to
// canopy-level yearly assimilation (mol yr^-1): 60*60 s/h * 12 h/day daylight *
// 365 d/yr / 1e6 umol/mol. Recurs in TF24_Strategy::net_mass_production_dt.
constexpr double tf24_assimilation_conv = 60.0 * 60.0 * 12.0 * 365.0 / 1e6;

// [eqn 2] Leaf area as a function of height (inverse of [eqn 3]). Carries the
// allometric traits a_l1, a_l2 -- the entry point through which an allometric
// trait reaches both the mass cascade and area_leaf itself.
template <typename S>
S tf24_area_leaf(S a_l1, S a_l2, S height) {
  using std::pow;
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 13] Total maintenance respiration (linear in the mass cascade).
template <typename S>
S tf24_respiration(S mass_leaf, S mass_sapwood, S mass_bark, S mass_root,
                   S r_l, S r_s, S r_b, S r_r) {
  return r_l * mass_leaf + r_b * mass_bark + r_s * mass_sapwood + r_r * mass_root;
}

// [eqn 14] Total turnover.
template <typename S>
S tf24_turnover(S mass_leaf, S mass_bark, S mass_sapwood, S mass_root,
                S k_l, S k_b, S k_s, S k_r) {
  return k_l * mass_leaf + k_b * mass_bark + k_s * mass_sapwood + k_r * mass_root;
}

// [eqn 15] Net production from assimilation/respiration/turnover.
template <typename S>
S tf24_net_production_A(S a_bio, S a_y, S assimilation, S respiration, S turnover) {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

// The TF24 parameters the net-production chain + demographic rate fill read, plus
// the prepare_strategy()-derived eta_c. The hot double path uses the per-piece
// functions directly with pars.* members; the all-in-one functions below (the AD
// calibration entry points) read these. NOTE: unlike FF16ProdPars there are no
// a_p1/a_p2 fields -- TF24 assimilation comes from the optimised leaf profit, not
// the light-response hyperbola.
template <typename S>
struct TF24ProdPars {
  S lma, rho, theta, a_b1, a_r1, eta_c;
  S r_l, r_s, r_b, r_r;
  S k_l, k_b, k_s, k_r;
  S a_bio, a_y;
  // Allometry + reproduction-allocation parameters for the growth/fecundity rates.
  S a_l1, a_l2;       // height <-> leaf-area allometry [eqn 2/3]
  S a_f1, a_f2, hmat; // reproduction-allocation logistic [eqn 16]
  S omega, a_f3;      // seed mass + accessory reproduction cost (fecundity_dt) [eqn 17]
  S d_I, a_dG1, a_dG2; // mortality: growth-independent + growth-dependent [eqn 21]
};

// Mass cascade -> respiration/turnover -> net production, GIVEN the assimilation
// rate. Shared by every assimilation variant (crown-centre, mean-light,
// deep-crown), which differ only in how `assimilation` is formed. Mirrors
// TF24_Strategy::net_mass_production_dt's tail exactly (FF16-identical algebra).
template <typename S>
S tf24_net_from_components(const TF24ProdPars<S>& p, S height, S area_leaf,
                          S assimilation) {
  const S mass_leaf    = area_leaf * p.lma;
  const S area_sapwood = area_leaf * p.theta;
  const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
  const S area_bark    = p.a_b1 * area_leaf * p.theta;
  const S mass_bark    = area_bark * height * p.eta_c * p.rho;
  const S mass_root    = p.a_r1 * area_leaf;

  const S respiration = tf24_respiration(mass_leaf, mass_sapwood, mass_bark, mass_root,
                                         p.r_l, p.r_s, p.r_b, p.r_r);
  const S turnover    = tf24_turnover(mass_leaf, mass_bark, mass_sapwood, mass_root,
                                      p.k_l, p.k_b, p.k_s, p.k_r);
  return tf24_net_production_A(p.a_bio, p.a_y, assimilation, respiration, turnover);
}

// Net production from the OPTIMISED leaf profit: assim = profit * area_leaf * conv,
// then the shared mass cascade. `profit` is the leaf submodel's profit_ (a max over
// collar potential); in the AD path it is an active scalar holding the injected
// Leaf::dprofit_d* sensitivities. Mirrors net_mass_production_dt's last lines.
template <typename S>
S tf24_net_mass_production(const TF24ProdPars<S>& p, S height, S area_leaf, S profit) {
  const S assimilation = profit * area_leaf * S(tf24_assimilation_conv);
  return tf24_net_from_components(p, height, area_leaf, assimilation);
}

// ---------------------------------------------------------------------------
// Demographic rate pieces (the rates downstream of net production). Mirror
// TF24_Strategy::{fraction_allocation_reproduction, fraction_allocation_growth,
// fecundity_dt, dheight_darea_leaf, dmass_*_darea_leaf, darea_leaf_dmass_live,
// mortality_growth_independent_dt, mortality_growth_dependent_dt} and the
// dheight/dt + heartwood assembly in compute_rates. Elementary arithmetic, so the
// whole demographic rate fill differentiates w.r.t. a trait.
// ---------------------------------------------------------------------------

// [eqn 16] Fraction of production allocated to reproduction (logistic in height).
template <typename S>
S tf24_fraction_allocation_reproduction(S a_f1, S a_f2, S hmat, S height) {
  using std::exp;
  return a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// [eqn 16] Fraction of production allocated to growth = 1 - reproduction.
template <typename S>
S tf24_fraction_allocation_growth(S a_f1, S a_f2, S hmat, S height) {
  return 1.0 - tf24_fraction_allocation_reproduction(a_f1, a_f2, hmat, height);
}

// [eqn 17] Rate of offspring production.
template <typename S>
S tf24_fecundity_dt(S net, S fraction_allocation_reproduction, S omega, S a_f3) {
  return net * fraction_allocation_reproduction / (omega + a_f3);
}

// d(height)/d(area_leaf): derivative of the [eqn 2] allometry.
template <typename S>
S tf24_dheight_darea_leaf(S a_l1, S a_l2, S area_leaf) {
  using std::pow;
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1.0);
}

// d(area_leaf)/d(mass_live): reciprocal of the summed per-component mass
// derivatives (leaf + sapwood + bark + root) w.r.t. area_leaf. Matches the
// TF24_Strategy::dmass_*_darea_leaf sum exactly.
template <typename S>
S tf24_darea_leaf_dmass_live(const TF24ProdPars<S>& p, S area_leaf) {
  using std::pow;
  const S dmass_leaf    = p.lma;                                   // d(area_leaf*lma)
  const S dmass_sapwood = p.rho * p.eta_c * p.a_l1 * p.theta *
                          (p.a_l2 + 1.0) * pow(area_leaf, p.a_l2);
  const S dmass_bark    = p.a_b1 * dmass_sapwood;
  const S dmass_root    = p.a_r1;
  return 1.0 / (dmass_leaf + dmass_sapwood + dmass_bark + dmass_root);
}

// dheight/dt given net production (the compute_rates growth assembly): returns 0
// when net is non-positive (the growth clamp), else dheight_darea_leaf *
// area_leaf_dt with area_leaf_dt = net * frac_growth * darea_leaf_dmass_live.
template <typename S>
S tf24_height_dt_from_net(const TF24ProdPars<S>& p, S height, S area_leaf, S net) {
  if (net <= 0.0) return S(0.0);
  const S frac_growth  = tf24_fraction_allocation_growth(p.a_f1, p.a_f2, p.hmat, height);
  const S darea_dmass  = tf24_darea_leaf_dmass_live(p, area_leaf);
  const S area_leaf_dt = net * frac_growth * darea_dmass;
  return tf24_dheight_darea_leaf(p.a_l1, p.a_l2, area_leaf) * area_leaf_dt;
}

// [eqn 21] growth-independent mortality (intrinsic baseline).
template <typename S>
S tf24_mortality_growth_independent_dt(S d_I) { return d_I; }

// [eqn 21] growth-dependent mortality; productivity_area = net / area_leaf.
template <typename S>
S tf24_mortality_growth_dependent_dt(S a_dG1, S a_dG2, S productivity_area) {
  using std::exp;
  return a_dG1 * exp(-a_dG2 * productivity_area);
}

// The five ODE state rates TF24_Strategy::compute_rates writes, plus the net
// production aux. Scalar-templated (#472 scope B) so the whole demographic rate
// fill differentiates w.r.t. a trait.
template <typename S>
struct TF24Rates {
  S net_mass_production_dt;
  S height_dt;
  S fecundity_dt;
  S area_heartwood_dt;
  S mass_heartwood_dt;
  S mortality_dt;
};

// Full TF24 compute_rates fill GIVEN the net production (however `net` was
// formed -- crown-centre / mean-light / deep-crown). Mirrors
// TF24_Strategy::compute_rates downstream of net EXACTLY: the net>0 growth clamp
// gates the growth/fecundity/heartwood rates, and mortality is the [eqn 21]
// growth-independent + growth-dependent sum (productivity = net/area_leaf).
// `mortality_finite` is the frozen util::is_finite(cumulative mortality) branch --
// a pass-1 (double) control-flow decision, passed in so the taped replay is
// branch-free (it never differentiates the is_finite test). This is the part of
// compute_rates the demographic kernel owns; the leaf optimisation upstream
// supplies `net`.
template <typename S>
TF24Rates<S> tf24_compute_rates_from_net(const TF24ProdPars<S>& p, S height,
                                         S area_leaf, S net,
                                         bool mortality_finite) {
  TF24Rates<S> r;
  r.net_mass_production_dt = net;
  if (net > 0.0) {
    const S frac_repro = tf24_fraction_allocation_reproduction(p.a_f1, p.a_f2,
                                                               p.hmat, height);
    r.height_dt          = tf24_height_dt_from_net(p, height, area_leaf, net);
    r.fecundity_dt       = tf24_fecundity_dt(net, frac_repro, p.omega, p.a_f3);
    const S area_sapwood = area_leaf * p.theta;            // [eqn 4]
    r.area_heartwood_dt  = p.k_s * area_sapwood;           // turnover of sapwood area
    const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
    r.mass_heartwood_dt  = p.k_s * mass_sapwood;           // turnover_sapwood(mass)
  } else {
    r.height_dt = S(0.0); r.fecundity_dt = S(0.0);
    r.area_heartwood_dt = S(0.0); r.mass_heartwood_dt = S(0.0);
  }
  // [eqn 21] instantaneous mortality rate; productivity_area = net / area_leaf.
  if (mortality_finite) {
    const S productivity_area = net / area_leaf;
    r.mortality_dt = tf24_mortality_growth_independent_dt(p.d_I) +
                     tf24_mortality_growth_dependent_dt(p.a_dG1, p.a_dG2, productivity_area);
  } else {
    r.mortality_dt = S(0.0);
  }
  return r;
}

}  // namespace plant

#endif
