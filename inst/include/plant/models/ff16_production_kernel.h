// -*-c++-*-
#ifndef PLANT_PLANT_FF16_PRODUCTION_KERNEL_H_
#define PLANT_PLANT_FF16_PRODUCTION_KERNEL_H_

#include <vector>
#include <cstddef>
#include <cmath>   // std::pow; XAD provides pow for active types via ADL

// Scalar-templated core of the FF16 single-plant net-mass-production chain
// (#472 scope B / traitecoevo/plant#537, Milestone A). The pieces below [eqn
// 12-15] are elementary arithmetic, so templating them on the scalar S makes net
// production differentiable w.r.t. traits by reverse-mode AD with no special
// handling. They are the SINGLE SOURCE OF TRUTH: FF16_Strategy's double methods
// delegate to them (so the existing FF16 test suite validates faithfulness), and
// the AD calibration path instantiates them with an active scalar.

namespace plant {

// [eqn 2] Leaf area as a function of height (inverse of [eqn 3]). Carries the
// allometric traits a_l1, a_l2, so it is the entry point through which a trait
// reaches both the mass cascade and the resident competition / light spline.
template <typename S>
S ff16_area_leaf(S a_l1, S a_l2, S height) {
  using std::pow;
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 12] Photosynthetic rate per leaf area; x is openness in [0,1].
template <typename S>
S ff16_assimilation_leaf(S a_p1, S a_p2, S x) {
  return a_p1 * x / (x + a_p2);
}

// [eqn 13] Total maintenance respiration (linear in the mass cascade).
template <typename S>
S ff16_respiration(S mass_leaf, S mass_sapwood, S mass_bark, S mass_root,
                   S r_l, S r_s, S r_b, S r_r) {
  return r_l * mass_leaf + r_b * mass_bark + r_s * mass_sapwood + r_r * mass_root;
}

// [eqn 14] Total turnover.
template <typename S>
S ff16_turnover(S mass_leaf, S mass_bark, S mass_sapwood, S mass_root,
                S k_l, S k_b, S k_s, S k_r) {
  return k_l * mass_leaf + k_b * mass_bark + k_s * mass_sapwood + k_r * mass_root;
}

// [eqn 15] Net production from assimilation/respiration/turnover.
template <typename S>
S ff16_net_production_A(S a_bio, S a_y, S assimilation, S respiration, S turnover) {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

// The FF16 parameters the crown-top net-production chain reads, plus the
// prepare_strategy()-derived eta_c. Used by the all-in-one convenience below
// (the AD calibration entry point); the hot double path uses the per-piece
// functions directly with pars.* members.
template <typename S>
struct FF16ProdPars {
  S lma, rho, theta, a_b1, a_r1, eta_c;
  S a_p1, a_p2;
  S r_l, r_s, r_b, r_r;
  S k_l, k_b, k_s, k_r;
  S a_bio, a_y;
  // Allometry + allocation parameters for the height-growth rate (Milestone C).
  S a_l1, a_l2;       // height <-> leaf-area allometry [eqn 2/3]
  S a_f1, a_f2, hmat; // reproduction-allocation logistic [eqn 16]
  // Demographic rate parameters for the full compute_rates fill (Milestone C):
  // fecundity [eqn 17] and mortality [eqn 21].
  S omega, a_f3;       // seed mass + accessory reproduction cost (fecundity_dt)
  S d_I, a_dG1, a_dG2; // mortality: growth-independent + growth-dependent
};

// Whole single-plant net production under the CROWN-TOP assimilation variant (a
// single light evaluation light_E at the canopy top; no quadrature -- the
// adaptive crown integral is a separate frozen-replay concern, Milestone B).
// Mirrors FF16_Strategy's mass cascade + assimilation_crown_top + the pieces
// above exactly. light_E is the (fixed, double on the resident path) light
// fraction at the crown.
// Mass cascade -> respiration/turnover -> net production, GIVEN the assimilation
// rate. Shared by every assimilation variant (crown-top, deep-crown, ...), so
// the variants differ only in how they compute `assimilation`.
template <typename S>
S ff16_net_from_components(const FF16ProdPars<S>& p, S height, S area_leaf,
                           S assimilation) {
  const S mass_leaf    = area_leaf * p.lma;
  const S area_sapwood = area_leaf * p.theta;
  const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
  const S area_bark    = p.a_b1 * area_leaf * p.theta;
  const S mass_bark    = area_bark * height * p.eta_c * p.rho;
  const S mass_root    = p.a_r1 * area_leaf;

  const S respiration = ff16_respiration(mass_leaf, mass_sapwood, mass_bark, mass_root,
                                         p.r_l, p.r_s, p.r_b, p.r_r);
  const S turnover    = ff16_turnover(mass_leaf, mass_bark, mass_sapwood, mass_root,
                                      p.k_l, p.k_b, p.k_s, p.k_r);
  return ff16_net_production_A(p.a_bio, p.a_y, assimilation, respiration, turnover);
}

template <typename S>
S ff16_net_mass_production_crown_top(const FF16ProdPars<S>& p,
                                     S height, S area_leaf, S light_E) {
  const S assimilation = area_leaf * ff16_assimilation_leaf(p.a_p1, p.a_p2, light_E);
  return ff16_net_from_components(p, height, area_leaf, assimilation);
}

// Deep-crown assimilation as a FROZEN-REPLAY weighted sum (#472 scope B,
// Milestone B). The production model's default assimilation integrates
// assimilation_leaf(light(z)) * q(z/height, z) over crown depth with adaptive
// Gauss-Kronrod. Per the two-pass plan we do NOT differentiate the adaptive
// controller: pass 1 (double) discovers the nodes z_j and the COMBINED weights
// wq_j = w_j * q(z_j/height, z_j) (the leaf-area density q is constant in the
// physiology traits, so it folds into the frozen weight); pass 2 replays
//   A = area_leaf * sum_j wq_j * assimilation_leaf(a_p1, a_p2, light(z_j)).
// `light` is any callable z -> S (e.g. an AD-capable resident light spline,
// odelia::interpolator::basic_interpolator<S>), so A is differentiable w.r.t.
// a_p1/a_p2 and w.r.t. light's knot values (the resident self-shading coupling).
template <typename S, typename LightFn>
S ff16_assimilation_deep_crown_replay(S a_p1, S a_p2, S area_leaf,
                                      const std::vector<double>& z,
                                      const std::vector<double>& wq,
                                      LightFn&& light) {
  S A = S(0.0);
  for (std::size_t j = 0; j < z.size(); ++j) {
    A += wq[j] * ff16_assimilation_leaf(a_p1, a_p2, S(light(z[j])));
  }
  return area_leaf * A;
}

// Total live mass given height (#472 scope B): the same leaf+sapwood+bark+root mass
// cascade as ff16_net_from_components, as a function of height. Mirrors
// FF16_Strategy::mass_live_given_height. The seedling height_0 solves
// mass_live_given_height(h0) = omega, so this is the residual whose root is height_0;
// differentiating it gives d(height_0)/d(trait) by the implicit function theorem
//   d(h0)/d(theta) = -(d mass_live/d theta - [theta==omega]) / (d mass_live/d height)
// at h0 -- the seedling-size response of the emergent gradient (the #539 IFT pattern).
template <typename S>
S ff16_mass_live_given_height(const FF16ProdPars<S>& p, S height) {
  const S area_leaf    = ff16_area_leaf(p.a_l1, p.a_l2, height);
  const S area_sapwood = area_leaf * p.theta;
  const S mass_leaf    = area_leaf * p.lma;
  const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
  const S mass_bark    = p.a_b1 * area_sapwood * height * p.eta_c * p.rho;
  const S mass_root    = p.a_r1 * area_leaf;
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

// Establishment probability (the recruitment filter), scalar-templated as a
// function of the SEEDLING net production (#472 scope B). Mirrors
// FF16_Strategy::establishment_probability:
//   pr_estab = decay_over_time / ((a_d0 * area_leaf_0 / net0)^2 + 1)   (net0 > 0),
// where net0 is the seedling's net production in the birth environment and
// decay_over_time = exp(-recruitment_decay * birth_time). recruitment_decay, the
// birth time and a_d0 are not physiology traits, so they fold to doubles; the trait
// dependence enters through net0 (and area_leaf_0). A node's initial mortality state
// is -log(pr_estab), so seeding the taped replay with -log(ff16_establishment_
// probability(...)) -- rather than a frozen constant -- makes the recruitment filter
// part of the emergent gradient. net0 > 0 is a frozen pass-1 sign.
template <typename S>
S ff16_establishment_probability(S area_leaf_0, S net0, double a_d0,
                                 double decay_over_time) {
  const S tmp = a_d0 * area_leaf_0 / net0;
  return decay_over_time / (tmp * tmp + 1.0);
}

// ---------------------------------------------------------------------------
// Height-growth rate pieces (#472 scope B, Milestone C). Mirror
// FF16_Strategy::{fraction_allocation_growth, dheight_darea_leaf,
// dmass_*_darea_leaf, darea_leaf_dmass_live} and the dheight/dt assembly in
// compute_rates. Elementary, so dheight/dt is differentiable w.r.t. height
// (A1 -- the exact growth-rate gradient now done by finite difference in
// Node::growth_rate_gradient) and w.r.t. traits.
// ---------------------------------------------------------------------------

// [eqn 16] Fraction of production allocated to reproduction (logistic in height).
template <typename S>
S ff16_fraction_allocation_reproduction(S a_f1, S a_f2, S hmat, S height) {
  using std::exp;
  return a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// [eqn 16] Fraction of production allocated to growth = 1 - reproduction.
template <typename S>
S ff16_fraction_allocation_growth(S a_f1, S a_f2, S hmat, S height) {
  return 1.0 - ff16_fraction_allocation_reproduction(a_f1, a_f2, hmat, height);
}

// d(height)/d(area_leaf): derivative of the [eqn 2] allometry.
template <typename S>
S ff16_dheight_darea_leaf(S a_l1, S a_l2, S area_leaf) {
  using std::pow;
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1.0);
}

// d(area_leaf)/d(mass_live): reciprocal of the summed per-component mass
// derivatives (leaf + sapwood + bark + root) w.r.t. area_leaf.
template <typename S>
S ff16_darea_leaf_dmass_live(const FF16ProdPars<S>& p, S area_leaf) {
  using std::pow;
  const S dmass_leaf    = p.lma;                                   // d(area_leaf*lma)
  const S dmass_sapwood = p.rho * p.eta_c * p.a_l1 * p.theta *
                          (p.a_l2 + 1.0) * pow(area_leaf, p.a_l2);
  const S dmass_bark    = p.a_b1 * dmass_sapwood;
  const S dmass_root    = p.a_r1;
  return 1.0 / (dmass_leaf + dmass_sapwood + dmass_bark + dmass_root);
}

// dheight/dt given net production (the compute_rates growth assembly): returns
// 0 when net is non-positive (the growth clamp), else dheight_darea_leaf *
// area_leaf_dt. Shared by every assimilation variant.
template <typename S>
S ff16_height_dt_from_net(const FF16ProdPars<S>& p, S height, S area_leaf, S net) {
  if (net <= 0.0) return S(0.0);
  const S frac_growth = ff16_fraction_allocation_growth(p.a_f1, p.a_f2, p.hmat, height);
  const S darea_dmass = ff16_darea_leaf_dmass_live(p, area_leaf);
  const S area_leaf_dt = net * frac_growth * darea_dmass;
  return ff16_dheight_darea_leaf(p.a_l1, p.a_l2, area_leaf) * area_leaf_dt;
}

// dheight/dt for a plant of the given height under the CROWN-TOP assimilation
// variant, in light light_E. area_leaf is derived from height so the gradient
// w.r.t. height flows through the whole chain.
template <typename S>
S ff16_height_dt_crown_top(const FF16ProdPars<S>& p, S height, S light_E) {
  const S area_leaf = ff16_area_leaf(p.a_l1, p.a_l2, height);
  const S net = ff16_net_mass_production_crown_top(p, height, area_leaf, light_E);
  return ff16_height_dt_from_net(p, height, area_leaf, net);
}

// The five ODE state rates FF16_Strategy::compute_rates writes, plus the net
// production aux. Scalar-templated (#472 scope B, Milestone C) so the whole
// demographic rate fill differentiates w.r.t. a trait by reverse-mode AD.
template <typename S>
struct FF16Rates {
  S net_mass_production_dt;
  S height_dt;
  S fecundity_dt;
  S area_heartwood_dt;
  S mass_heartwood_dt;
  S mortality_dt;
};

// Full FF16 compute_rates fill for the CROWN-TOP assimilation variant (single
// light evaluation light_E). Mirrors FF16_Strategy::compute_rates EXACTLY: the
// net>0 growth clamp gates the growth/fecundity/heartwood rates, and mortality
// is the [eqn 21] growth-independent + growth-dependent sum (productivity =
// net/area_leaf). `mortality_finite` is the frozen util::is_finite(cumulative
// mortality) branch -- a pass-1 (double) control-flow decision, passed in so the
// taped replay is branch-free (it never differentiates the is_finite test).
// Deep-crown differs only in how `net` is formed (the frozen-replay crown
// integral, ff16_assimilation_deep_crown_replay -> ff16_net_from_components),
// so a deep-crown fill reuses everything below by substituting that `net`.
// The rate fill SHARED by every assimilation variant: given the net production
// (however it was formed -- single-light crown-top, or the deep-crown crown
// integral) and the area_leaf, write the five ODE rates. This is the part of
// FF16_Strategy::compute_rates downstream of `net`; the net>0 growth clamp gates
// growth/fecundity/heartwood, mortality is the [eqn 21] sum (productivity =
// net/area_leaf). `mortality_finite` is the frozen is_finite branch (a pass-1
// control-flow decision, so the taped replay never differentiates the test).
template <typename S>
FF16Rates<S> ff16_compute_rates_from_net(const FF16ProdPars<S>& p, S height,
                                         S area_leaf, S net,
                                         bool mortality_finite) {
  FF16Rates<S> r;
  r.net_mass_production_dt = net;
  if (net > 0.0) {
    const S frac_repro = ff16_fraction_allocation_reproduction(p.a_f1, p.a_f2,
                                                               p.hmat, height);
    r.height_dt          = ff16_height_dt_from_net(p, height, area_leaf, net);
    r.fecundity_dt       = net * frac_repro / (p.omega + p.a_f3);
    const S area_sapwood = area_leaf * p.theta;            // [eqn 4]
    r.area_heartwood_dt  = p.k_s * area_sapwood;           // turnover of sapwood area
    const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
    r.mass_heartwood_dt  = p.k_s * mass_sapwood;           // turnover_sapwood(mass)
  } else {
    r.height_dt = S(0.0); r.fecundity_dt = S(0.0);
    r.area_heartwood_dt = S(0.0); r.mass_heartwood_dt = S(0.0);
  }
  // [eqn 21] instantaneous mortality rate; productivity_area = net / area_leaf.
  using std::exp;
  if (mortality_finite) {
    const S productivity_area = net / area_leaf;
    r.mortality_dt = p.d_I + p.a_dG1 * exp(-p.a_dG2 * productivity_area);
  } else {
    r.mortality_dt = S(0.0);
  }
  return r;
}

template <typename S>
FF16Rates<S> ff16_compute_rates_crown_top(const FF16ProdPars<S>& p, S height,
                                          S light_E, bool mortality_finite) {
  const S area_leaf = ff16_area_leaf(p.a_l1, p.a_l2, height);
  const S net = ff16_net_mass_production_crown_top(p, height, area_leaf, light_E);
  return ff16_compute_rates_from_net(p, height, area_leaf, net, mortality_finite);
}

// [eqn] Yokozawa leaf-area density q(z,H) = 2 eta (1 - u^eta) u^eta / z,
// u = z/H. Mirrors CanopyShape::q exactly (eta is a fixed double; the gradient
// flows through the active u and z). The deep-crown crown integral weights the
// per-depth assimilation by this.
template <typename S>
S ff16_canopy_q(double eta, S u, S z) {
  using std::pow;
  const S u_eta = pow(u, eta);
  return 2.0 * eta * (1.0 - u_eta) * u_eta / z;
}

// Single-plant height TRAJECTORY: integrate dheight/dt from h0 over n fixed RK4
// steps to age t_end, in a fixed crown-top light light_E (#472 scope B). This is
// the bridge from instantaneous-rate gradients to time-integrated emergent
// outputs: the whole trajectory is scalar-templated, so reverse/forward AD gives
// d(height at age t_end)/d(trait) -- a calibration gradient through the growth
// ODE. A fixed step schedule is exactly the frozen-schedule formulation
// end-to-end AD needs (the adaptive stepper is a pass-1 discovery, replayed).
template <typename S>
S ff16_grow_height(const FF16ProdPars<S>& p, S h0, S light_E,
                   double t_end, int n_steps) {
  S h = h0;
  const double dt = t_end / n_steps;
  for (int i = 0; i < n_steps; ++i) {
    // Materialise each stage as S: with XAD expression templates `h + c*k`
    // is an expression type, not S, which would break template deduction of
    // ff16_height_dt_crown_top's scalar.
    const S k1 = ff16_height_dt_crown_top(p, h, light_E);
    const S h2 = h + S(0.5 * dt) * k1;
    const S k2 = ff16_height_dt_crown_top(p, h2, light_E);
    const S h3 = h + S(0.5 * dt) * k2;
    const S k3 = ff16_height_dt_crown_top(p, h3, light_E);
    const S h4 = h + S(dt) * k3;
    const S k4 = ff16_height_dt_crown_top(p, h4, light_E);
    h = h + S(dt / 6.0) * (k1 + S(2.0) * k2 + S(2.0) * k3 + k4);
  }
  return h;
}

// The five FF16 ODE states, as a value the trajectory integrator carries.
template <typename S>
struct FF16State {
  S height, mortality, fecundity, area_heartwood, mass_heartwood;
};

// Single-plant FULL-STATE demographic TRAJECTORY (#472 scope B, Milestone C):
// integrate all five FF16 states from y0 over n_steps fixed RK4 steps to age
// t_end, in a fixed crown-top light light_E, using ff16_compute_rates_crown_top.
// Generalises ff16_grow_height (height only) to the demographic vector, so
// reverse AD gives d(any emergent state at age t_end)/d(trait) -- e.g. lifetime
// fecundity (cumulative offspring) sensitivity, a calibration target through the
// whole demographic ODE. mortality_finite is the frozen pass-1 branch (see the
// rate kernel). Each RK4 stage is materialised as S so XAD expression templates
// don't break scalar deduction (same caveat as ff16_grow_height).
template <typename S>
FF16State<S> ff16_grow_demography(const FF16ProdPars<S>& p, FF16State<S> y,
                                  S light_E, double t_end, int n_steps,
                                  bool mortality_finite) {
  const double dt = t_end / n_steps;
  auto deriv = [&](const FF16State<S>& s) -> FF16State<S> {
    const FF16Rates<S> r =
        ff16_compute_rates_crown_top(p, s.height, light_E, mortality_finite);
    return FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
                        r.area_heartwood_dt, r.mass_heartwood_dt};
  };
  auto axpy = [](const FF16State<S>& a, S c, const FF16State<S>& k) -> FF16State<S> {
    return FF16State<S>{a.height + c * k.height, a.mortality + c * k.mortality,
                        a.fecundity + c * k.fecundity,
                        a.area_heartwood + c * k.area_heartwood,
                        a.mass_heartwood + c * k.mass_heartwood};
  };
  for (int i = 0; i < n_steps; ++i) {
    const FF16State<S> k1 = deriv(y);
    const FF16State<S> k2 = deriv(axpy(y, S(0.5 * dt), k1));
    const FF16State<S> k3 = deriv(axpy(y, S(0.5 * dt), k2));
    const FF16State<S> k4 = deriv(axpy(y, S(dt), k3));
    const S c = S(dt / 6.0);
    y.height = y.height + c * (k1.height + S(2.0) * k2.height + S(2.0) * k3.height + k4.height);
    y.mortality = y.mortality + c * (k1.mortality + S(2.0) * k2.mortality + S(2.0) * k3.mortality + k4.mortality);
    y.fecundity = y.fecundity + c * (k1.fecundity + S(2.0) * k2.fecundity + S(2.0) * k3.fecundity + k4.fecundity);
    y.area_heartwood = y.area_heartwood + c * (k1.area_heartwood + S(2.0) * k2.area_heartwood + S(2.0) * k3.area_heartwood + k4.area_heartwood);
    y.mass_heartwood = y.mass_heartwood + c * (k1.mass_heartwood + S(2.0) * k2.mass_heartwood + S(2.0) * k3.mass_heartwood + k4.mass_heartwood);
  }
  return y;
}

// Frozen-schedule forward-Euler replay of ONE cohort's full demographic state
// over a per-step crown-light schedule (#472 scope B, Milestone C -- the two-pass
// SCM replay primitive). `light[k]` is the (frozen, pass-1 double) light the
// cohort's crown reads at global replay step k; the cohort is born at step0 and
// integrated to the end of the schedule with forward Euler at fixed dt. Euler
// (not RK4) is deliberate: it reproduces the SCM's control.fixed_time_step
// integration EXACTLY, so replaying a fixed_time_step resident run is faithful.
// Templated on S, so reverse AD over a weighted sum of cohort outcomes
//   J(theta) = sum_i w_i * f(replay_i)   (w_i, light frozen from pass 1)
// gives d(emergent stand output)/d(trait), holding the resident light schedule
// fixed -- the legitimate "resident-light-frozen" gradient. The full resident
// self-shading gradient additionally makes light active (odelia #32 active-query
// spline / ff16_assimilation_deep_crown_replay), deferred.
template <typename S>
FF16State<S> ff16_replay_cohort(const FF16ProdPars<S>& p, FF16State<S> y,
                                double dt, const std::vector<double>& light,
                                std::size_t step0, bool mortality_finite) {
  for (std::size_t k = step0; k < light.size(); ++k) {
    const FF16Rates<S> r =
        ff16_compute_rates_crown_top(p, y.height, S(light[k]), mortality_finite);
    y.height         = y.height + S(dt) * r.height_dt;
    y.mortality      = y.mortality + S(dt) * r.mortality_dt;
    y.fecundity      = y.fecundity + S(dt) * r.fecundity_dt;
    y.area_heartwood = y.area_heartwood + S(dt) * r.area_heartwood_dt;
    y.mass_heartwood = y.mass_heartwood + S(dt) * r.mass_heartwood_dt;
  }
  return y;
}

// As ff16_replay_cohort, but the cohort reads its crown light ACTIVELY from a
// frozen resident profile at each step (#472 scope B, Milestone C). `crown_light`
// is a caller-supplied callable S -> S returning the light at the cohort's crown
// for its current (active) height; in the AD context the caller seeds it from the
// resident profile's value + slope (FF16_Environment::get_environment_at_height /
// get_environment_deriv_at_height) so d(light)/d(height) flows -- the within-cohort
// self-shading feedback (a taller cohort reads higher in the canopy -> more light
// -> faster growth) that the frozen per-step light of the plain overload omits.
// The profile KNOTS stay frozen double (resident held fixed); making them active
// is the full self-shading gradient (odelia #32 active-knot spline). LightFn is a
// template so XAD never enters this header (mirrors ff16_assimilation_deep_crown_replay).
template <typename S, typename LightFn>
FF16State<S> ff16_replay_cohort_active_light(const FF16ProdPars<S>& p, FF16State<S> y,
                                             double dt, LightFn&& crown_light,
                                             int n_steps, bool mortality_finite) {
  for (int k = 0; k < n_steps; ++k) {
    const S light_E = crown_light(y.height);
    const FF16Rates<S> r =
        ff16_compute_rates_crown_top(p, y.height, light_E, mortality_finite);
    y.height         = y.height + S(dt) * r.height_dt;
    y.mortality      = y.mortality + S(dt) * r.mortality_dt;
    y.fecundity      = y.fecundity + S(dt) * r.fecundity_dt;
    y.area_heartwood = y.area_heartwood + S(dt) * r.area_heartwood_dt;
    y.mass_heartwood = y.mass_heartwood + S(dt) * r.mass_heartwood_dt;
  }
  return y;
}

// PRODUCTION two-pass replay primitive: a single cohort's full demographic state
// integrated with the SAME adaptive Cash-Karp RKCK scheme the live SCM used
// (#472 scope B, Milestone C -- the FAITHFUL replacement for the Euler
// ff16_replay_cohort). Forward Euler mirrors only the non-default
// control.fixed_time_step path; the real SCM integrates with odelia's embedded
// 4/5 RKCK (ode_step.hpp), and the mutant-fitness replay (run_mutant ->
// advance_fixed -> step_to) re-uses that SAME stepper over the resident's pinned
// step times, swapping in a FROZEN per-RK-stage environment (environment_history
// [step][stage], 6 stages/step). This kernel is that path lifted to the scalar S.
//
// step_h[n] = the resident's actual adaptive step size for global step n (=
// step_history[n+1]-step_history[n]); the cohort is integrated from global step
// `step0` (its birth step) to the end of the schedule. The per-stage crown light
// is supplied by `crown_light(n, stage, height) -> S` so XAD/odelia stay out of
// this header (mirrors ff16_replay_cohort_active_light). The caller wires it to
// the FROZEN resident env: in the AD context it seeds value + slope from
// FF16_Environment::get_environment_at_height / get_environment_deriv_at_height
// so d(light)/d(height) flows -- the mutant-through-frozen-canopy feedback (a
// taller focal cohort reads higher in the resident profile). Stage codes:
//   stage 0 -> k1 env  = the env at the step START (environment_history[n-1][5],
//                        or the birth env for n==step0); recomputing k1 against
//                        it is numerically identical to the solver's FSAL reuse;
//   stage 1..5 -> the envs for the k2..k6 derivs (environment_history[n][0..4]).
// The 6th cached env (environment_history[n][5], the solver's dydt_out stage) is
// re-used as the next step's stage-0 env, exactly as first_same_as_last does. The
// y-update uses k1,k3,k4,k6 (c2==c5==0), matching ode_step.hpp::step line-for-line.
// Generic Cash-Karp RKCK driver over the FROZEN resident schedule, shared by the
// demographic replay (ff16_replay_cohort_rkck) and the lifetime-offspring replay
// (ff16_replay_cohort_offspring_rkck). The integrator logic -- the GSL Cash-Karp
// tableau, the FSAL stage-0 reuse, the c2==c5==0 final sum -- lives here ONCE.
// Callers supply the state type and its two operations:
//   deriv(state, n, stage) -> State : the state-derivative at RK `stage` (0..5) of
//        global step n. stage 0 is the FSAL k1 (evaluated at the step START);
//        stages 1..5 are the k2..k6 derivs. Callers map (n, stage) to the frozen
//        per-RK-stage resident environment (see ff16_replay_cohort_rkck).
//   axpy(a, c, k) -> State : a + c*k with c a double RK coefficient; each component
//        MUST be materialised into the scalar type in the returned State (brace-
//        init), so XAD expression templates never escape with dangling references.
// The y-update is the 5th-order sum y += h*(c1 k1 + c3 k3 + c4 k4 + c6 k6) written
// as an axpy chain (c2==c5==0), matching odelia ode_step.hpp::step.
template <typename State, typename DerivFn, typename AxpyFn>
State ff16_cashkarp_replay(State y, const std::vector<double>& step_h,
                           std::size_t step0, DerivFn&& deriv, AxpyFn&& axpy) {
  // Cash-Karp coefficients, identical to odelia::ode::Step (from GSL).
  const double b21   = 1.0 / 5.0;
  const double b3[2] = {3.0 / 40.0, 9.0 / 40.0};
  const double b4[3] = {0.3, -0.9, 1.2};
  const double b5[4] = {-11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0};
  const double b6[5] = {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0,
                        44275.0 / 110592.0, 253.0 / 4096.0};
  const double c1 = 37.0 / 378.0, c3 = 250.0 / 621.0,
               c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0;

  for (std::size_t n = step0; n < step_h.size(); ++n) {
    const double h = step_h[n];
    const State k1 = deriv(y, n, 0);
    const State k2 = deriv(axpy(y, b21 * h, k1), n, 1);
    State y3 = axpy(y, h * b3[0], k1); y3 = axpy(y3, h * b3[1], k2);
    const State k3 = deriv(y3, n, 2);
    State y4 = axpy(y, h * b4[0], k1);
    y4 = axpy(y4, h * b4[1], k2); y4 = axpy(y4, h * b4[2], k3);
    const State k4 = deriv(y4, n, 3);
    State y5 = axpy(y, h * b5[0], k1);
    y5 = axpy(y5, h * b5[1], k2); y5 = axpy(y5, h * b5[2], k3); y5 = axpy(y5, h * b5[3], k4);
    const State k5 = deriv(y5, n, 4);
    State y6 = axpy(y, h * b6[0], k1);
    y6 = axpy(y6, h * b6[1], k2); y6 = axpy(y6, h * b6[2], k3);
    y6 = axpy(y6, h * b6[3], k4); y6 = axpy(y6, h * b6[4], k5);
    const State k6 = deriv(y6, n, 5);
    y = axpy(axpy(axpy(axpy(y, h * c1, k1), h * c3, k3), h * c4, k4), h * c6, k6);
  }
  return y;
}

template <typename S, typename StageLightFn>
FF16State<S> ff16_replay_cohort_rkck(const FF16ProdPars<S>& p, FF16State<S> y,
                                     const std::vector<double>& step_h,
                                     std::size_t step0,
                                     StageLightFn&& crown_light,
                                     bool mortality_finite) {
  // 5-state FF16 derivative at a trial state, reading the frozen stage env.
  auto deriv = [&](const FF16State<S>& s, std::size_t n, int stage) -> FF16State<S> {
    const S light_E = crown_light(n, stage, s.height);
    const FF16Rates<S> r =
        ff16_compute_rates_crown_top(p, s.height, light_E, mortality_finite);
    return FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
                        r.area_heartwood_dt, r.mass_heartwood_dt};
  };
  // a + c*k, each component materialised as S in the brace-init.
  auto axpy = [](const FF16State<S>& a, double c, const FF16State<S>& k) -> FF16State<S> {
    return FF16State<S>{a.height + c * k.height, a.mortality + c * k.mortality,
                        a.fecundity + c * k.fecundity,
                        a.area_heartwood + c * k.area_heartwood,
                        a.mass_heartwood + c * k.mass_heartwood};
  };
  return ff16_cashkarp_replay(y, step_h, step0, deriv, axpy);
}

// State for the lifetime-offspring replay: the 5 FF16 states + the cumulative
// survival-weighted offspring (the SCM's offspring_produced_survival_weighted).
template <typename S>
struct FF16LifeState {
  FF16State<S> demog;
  S offspring;
};

// Lifetime survival-weighted offspring replay (#472 scope B): augments the
// demographic replay with a 6th accumulator mirroring Node::compute_rates,
//   d(offspring)/dt = fecundity_dt * exp(-mortality) * surv_weight(n, stage),
// integrated with the SAME Cash-Karp driver. `surv_weight(n, stage) -> double`
// supplies the FROZEN pr_patch_survival(t_stage)/pr_patch_survival_at_birth at the
// step's RK stage time; set y.demog.mortality = -log(establishment_probability) at
// birth (the node's initial condition) before calling. The stand's emergent
// offspring_production is then the node-spacing trapezium of
//   offspring * patch_density_at_birth * S_D * birth_rate
// over the cohorts (a frozen, linear post-weighting), so one reverse sweep of that
// sum gives d(offspring_production)/d(trait). crown_light/surv_weight are callables
// so XAD/odelia stay out of this header.
template <typename S, typename StageLightFn, typename SurvFn>
FF16LifeState<S> ff16_replay_cohort_offspring_rkck(
    const FF16ProdPars<S>& p, FF16LifeState<S> y, const std::vector<double>& step_h,
    std::size_t step0, StageLightFn&& crown_light, SurvFn&& surv_weight,
    bool mortality_finite) {
  using std::exp;  // XAD provides exp for active types via ADL
  auto deriv = [&](const FF16LifeState<S>& s, std::size_t n, int stage) -> FF16LifeState<S> {
    const S light_E = crown_light(n, stage, s.demog.height);
    const FF16Rates<S> r =
        ff16_compute_rates_crown_top(p, s.demog.height, light_E, mortality_finite);
    const S off_dt = r.fecundity_dt * exp(-s.demog.mortality) * S(surv_weight(n, stage));
    return FF16LifeState<S>{FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
                                         r.area_heartwood_dt, r.mass_heartwood_dt},
                            off_dt};
  };
  auto axpy = [](const FF16LifeState<S>& a, double c, const FF16LifeState<S>& k) -> FF16LifeState<S> {
    return FF16LifeState<S>{
        FF16State<S>{a.demog.height + c * k.demog.height,
                     a.demog.mortality + c * k.demog.mortality,
                     a.demog.fecundity + c * k.demog.fecundity,
                     a.demog.area_heartwood + c * k.demog.area_heartwood,
                     a.demog.mass_heartwood + c * k.demog.mass_heartwood},
        a.offspring + c * k.offspring};
  };
  return ff16_cashkarp_replay(y, step_h, step0, deriv, axpy);
}

// Resident light availability E(z) = exp( - sum_i density_i * k_I * area_leaf_i *
// Q(z/h_i) ) at height z from a FROZEN stand (heights/densities are pass-1
// doubles), ACTIVE in the traits through each cohort's area_leaf [eqn 2]
// (#472 scope B, Milestone C -- the resident self-shading coupling). Q is the
// deep/Yokozawa leaf-area-above (1 - u^eta)^2 with eta a fixed double (other
// shading variants swap Q); contributions vanish above each plant's top. Beer's
// law E = exp(-projected leaf area), matching FF16_Environment::compute_environment.
// Evaluated at FROZEN knot positions z_k to fill an active-VALUE light spline
// (odelia basic_interpolator<S>), this is what makes a resident emergent output
// differentiable w.r.t. a trait THROUGH the self-shaded light profile -- the
// full self-shading gradient (vs the frozen-knot active-query of
// ff16_replay_cohort_active_light). Heights frozen here (a fixed stand census);
// coupling growth back in is the live two-pass replay.
template <typename S>
S ff16_resident_light_at(double z, S a_l1, S a_l2, double k_I, double eta,
                         const std::vector<double>& height,
                         const std::vector<double>& density) {
  using std::pow; using std::exp;
  S L = S(0.0);
  for (std::size_t i = 0; i < height.size(); ++i) {
    if (z >= height[i]) continue;            // no leaf area above the plant's crown
    const double u = z / height[i];
    const double one_minus = 1.0 - pow(u, eta);
    const double Q = one_minus * one_minus;  // Yokozawa leaf-area-above
    L += S(density[i] * k_I) * ff16_area_leaf(a_l1, a_l2, S(height[i])) * S(Q);
  }
  return exp(-L);
}

}  // namespace plant

#endif
