// Reverse-mode AD emergent-gradient routines for FF16, compiled into plant.so
// (#472 scope B, Phase C). Unlike growth_rate_gradient_height_ad (forward mode,
// header-only), these use the XAD adjoint TAPE (xad::adj). plant.so has no tape of
// its own: the tape symbols are odelia's single compiled copy (src/Tape.cpp),
// resolved at load against odelia's globally-loaded DLL -- the same mechanism the
// odelia ODE Solver already relies on (see src/Makevars; odelia is imported first
// via importFrom(odelia, odelia_load_dll), so its DLL loads before plant's).
//
// The headline routine differentiates the SCM's emergent offspring_production
// w.r.t. a set of FF16 traits in ONE reverse sweep, over the frozen resident
// schedule + per-RK-stage resident light harvested by a save_RK45_cache run
// (deep-crown / default shading). Establishment is frozen (a separable partial).
// It takes the harvested data as plain arrays, so it carries no templating into the
// SCM class; the R-facing offspring_production_gradient() gathers these from a run
// SCM and calls it.
#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <functional>
#include <XAD/XAD.hpp>
#include <plant.h>                 // RcppR6 as<>/wrap for FF16_Environment etc.
#include <plant/models/ff16_production_kernel.h>
#include <plant/gradient/scm_harvest.h> // shared birth_steps / recover_birth_rate / census_trapezium
#include <plant/gradient/coupled_canopy.h> // shared Yokozawa canopy_comp_at
#include <odelia/interpolator.hpp> // basic_interpolator<S> (active-value light spline)

using ad   = xad::adj<double>;
using ad_t = ad::active_type;
using fwd_t = xad::fwd<double>::active_type;   // forward (tangent) mode, header-only

namespace {

double as_double(double v)      { return v; }
double as_double(const ad_t& v) { return xad::value(v); }
double as_double(const fwd_t& v){ return xad::value(v); }

plant::FF16_Strategy make_strategy(const Rcpp::NumericVector& pp) {
  plant::FF16_Strategy s; auto& q = s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];q.a_p1=pp["a_p1"];q.a_p2=pp["a_p2"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.S_D=pp["S_D"];q.a_d0=pp["a_d0"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.recruitment_decay=pp["recruitment_decay"];
  s.prepare_strategy(); return s;
}
template <typename S> plant::FF16ProdPars<S> lift(const plant::FF16ProdPars<double>& d) {
  plant::FF16ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.a_p1=d.a_p1;p.a_p2=d.a_p2;p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;
  p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;p.a_bio=d.a_bio;p.a_y=d.a_y;
  p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2; return p;
}
template <typename S> std::vector<S*> field_ptrs(plant::FF16ProdPars<S>& p) {
  return {&p.lma,&p.rho,&p.theta,&p.a_b1,&p.a_r1,&p.eta_c,&p.a_p1,&p.a_p2,
          &p.r_l,&p.r_s,&p.r_b,&p.r_r,&p.k_l,&p.k_b,&p.k_s,&p.k_r,&p.a_bio,&p.a_y,
          &p.a_l1,&p.a_l2,&p.a_f1,&p.a_f2,&p.hmat,&p.omega,&p.a_f3,&p.d_I,&p.a_dG1,&p.a_dG2};
}
std::vector<std::string> field_names() {
  return {"lma","rho","theta","a_b1","a_r1","eta_c","a_p1","a_p2","r_l","r_s","r_b","r_r",
          "k_l","k_b","k_s","k_r","a_bio","a_y","a_l1","a_l2","a_f1","a_f2","hmat","omega",
          "a_f3","d_I","a_dG1","a_dG2"};
}

struct Frozen {
  std::vector<std::vector<plant::FF16_Environment>> eh;  // [step][0..5]
  std::vector<double> step_h, ppsab, tw, decay;          // decay = exp(-recr_decay*t_birth)
  Rcpp::NumericMatrix ppsurv;                            // [step][0..5] stage survival
  std::vector<int> birth;
  double eta, h0, a_d0;
  // Harvest birth_rate, ONLY for the net_reproduction_ratio unit weights w_i =
  // tw_i / birth_rate0 (tw is harvested birth-scaled; R0 uses per-seed unit weights). It
  // must stay the FIXED harvest value so that perturbing/activating `birth_rate` (which
  // scales the resident density / canopy) does NOT spuriously rescale the weight -- the
  // mutant reading: fixed per-seed weight, varying resident density. <=0 => fall back to
  // the value of the birth_rate passed to assemble (correct for the native AD, whose
  // active base IS the harvest value; the R-list FD reference passes it explicitly).
  double birth_rate0 = -1.0;
  const plant::quadrature::QK* integ;

  // RESIDENT total-gradient harvest (#472 scope B, R0). When `resident` is true the
  // cohorts read a value-ANCHORED active light: the frozen env VALUE plus the
  // theta-derivative of a trapezium reconstruction over the per-RK-stage frozen
  // stand. st_h[n][s] / st_C[n][s] are the species-0 stand at step n, RK stage s,
  // sorted by DESCENDING height; st_C is the frozen per-node weight
  // C_i = ce_i / area_leaf(theta0, h_i) so that C_i * area_leaf(theta, h_i) flows
  // the allometric trait while reproducing ce_i at theta0. `area` is the patch area
  // (the Beer's-law competition is trapezium/area, matching Patch::compute_competition).
  bool resident = false;
  // anchor=true (the shipped resident path): cohorts read the frozen env VALUE with
  // the recon theta-DERIVATIVE added (zero value) -> baseline metrics bit-identical
  // to the frozen engine. anchor=false: cohorts read the genuine recon VALUE
  // (theta-dependent, ~1e-4 off the frozen env) -- used only to FD-validate R1
  // ("re-run the reconstruction, re-reduce"); AD and FD then differentiate the SAME
  // function and agree tightly, and AD(anchored) ~ AD(noanchor) confirms anchoring
  // is gradient-neutral while it nails the baseline value.
  bool anchor = true;
  double area = 1.0;
  std::vector<std::vector<std::vector<double>>> st_h, st_C;

  // COUPLED resident replay (#472 scope B, R0-R1, the course-corrected build). The
  // whole stand is re-evolved together over the frozen schedule and the canopy light
  // is reconstructed per RK stage from the ACTIVE stand (heights AND densities respond
  // to theta), so EVERY trait feeds back -- not just the leaf-area channel of the
  // frozen-geometry graft above. coupled=true selects it.
  //   kI            : light extinction coefficient (competition = kI*leaf-area*Q).
  //   knot_x[n][s]  : the env light-spline's frozen knot x-positions at step n,
  //                   cached stage s (the "freeze knot positions" of odelia #32).
  //   knot_y0[n][s] : the SCM's stored knot light VALUES there (the R0 ground truth
  //                   the active reconstruction must reproduce).
  //   nn_h/nn_c[n][s]: the boundary new_node height + competition effect at (n,s)
  //                   (FROZEN; the seedling tail term, a tiny ground-level channel).
  bool coupled = false;
  double kI = 0.0;
  std::vector<std::vector<std::vector<double>>> knot_x, knot_y0;
  std::vector<std::vector<double>> nn_h, nn_c;
  // ACTIVE birth-env establishment (#472 scope B, R1 step 2). When true the cohort
  // initial conditions (net0/pr_estab/mort0/g0/logd0) read the canopy RECONSTRUCTED
  // from the already-alive cohorts at the start of the birth step (active in theta),
  // not the FROZEN harvested birth env -- so a trait feeds back through establishment
  // too. Exact at theta0 (the alive set there == the harvested stand). false restores
  // the frozen-birth-env behaviour (used to A/B-measure the birth-env channel).
  bool coupled_active_birthenv = true;
};

// Deep-crown net at `height` reading the frozen env `e` (moving-node GK integral).
template <typename S>
S deep_net(const plant::FF16ProdPars<S>& pd, const plant::quadrature::QK* integ,
           double eta, const plant::FF16_Environment* e, S height) {
  const double canopy_top = e->max_environment_height();
  auto integrand = [&](S z) -> S {
    double zv = as_double(z);
    double lv = e->get_environment_at_height(zv, canopy_top);
    double ld = e->get_environment_deriv_at_height(zv);
    S light = S(lv) + S(std::isfinite(ld)?ld:0.0) * (z - S(zv));
    return plant::ff16_assimilation_leaf<S>(pd.a_p1, pd.a_p2, light) *
           plant::ff16_canopy_q<S>(eta, z / height, z);
  };
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S assim = area_leaf * integ->integrate_ad<S>(integrand, S(0.0), height);
  return plant::ff16_net_from_components(pd, height, area_leaf, assim);
}

// Resident competition at z reconstructed from a FROZEN per-RK-stage stand
// (descending heights h, frozen weights C_i = ce_i/area_leaf(theta0,h_i)), ACTIVE
// in the allometric trait through area_leaf(theta, h_i). This is the C-27/C-28
// trapezium matching Species::compute_competition: comp(z) = (1/2) sum_adjacent
// (h_i - h_{i+1})(g_i + g_{i+1}), g_i = C_i*area_leaf(theta,h_i)*Q(z/h_i), Q the
// Yokozawa leaf-area-above (1-u^eta)^2. Returned UN-divided by area (caller divides).
template <typename S>
S recon_competition(S z, const plant::FF16ProdPars<S>& pd, double eta,
                    const std::vector<double>& h, const std::vector<double>& C) {
  using std::pow;
  if (h.size() < 2) return S(0.0);               // <2 nodes: SCM trapezium is 0
  auto g = [&](std::size_t i) -> S {
    if (as_double(z) >= h[i]) return S(0.0);     // no leaf area above the crown
    const S u = z / S(h[i]);                     // z may be active (focal-height channel)
    const S om = S(1.0) - pow(u, S(eta));
    return S(C[i]) * plant::ff16_area_leaf(pd.a_l1, pd.a_l2, S(h[i])) * (om * om);
  };
  S comp = S(0.0);
  S gp = g(0); double hp = h[0];
  for (std::size_t i = 1; i < h.size(); ++i) {
    S gi = g(i);
    comp = comp + S(hp - h[i]) * (gp + gi);
    hp = h[i]; gp = gi;
  }
  return S(0.5) * comp;
}

// Value-ANCHORED resident light at zv: the exact frozen-env VALUE plus the
// theta-derivative of the trapezium reconstruction (zero value). The baseline light
// is therefore bit-identical to the frozen engine; only the resident feedback
// channel theta->canopy->light is added. zv is a double (the GK integration point);
// the focal-height->light channel is carried by deep_net's frozen z-linearisation.
template <typename S>
S resident_light_anchored(double zv, double lv, const plant::FF16ProdPars<S>& pd,
                          double eta, const std::vector<double>& h,
                          const std::vector<double>& C, double area) {
  using std::exp;
  // recon at the FROZEN point zv (constant): the focal-height->light channel is
  // carried separately by deep_net's analytic frozen-env z-derivative, so here only
  // the trait (a_l1/a_l2) feedback flows. Value-anchored to the frozen env light lv.
  S comp  = recon_competition<S>(S(zv), pd, eta, h, C) * S(1.0 / area);
  S recon = exp(-comp);
  return S(lv) + (recon - S(as_double(recon)));
}

// Deep-crown net at `height` reading the value-anchored RESIDENT light at stage
// (n,stage): the frozen env e + the active reconstruction over the matching frozen
// per-RK-stage stand (sh heights, sC weights). Mirrors deep_net but with the light
// active in the allometric trait through the canopy (the resident feedback term).
template <typename S>
S deep_net_resident(const plant::FF16ProdPars<S>& pd, const plant::quadrature::QK* integ,
                    double eta, const plant::FF16_Environment* e, S height,
                    const std::vector<double>& sh, const std::vector<double>& sC,
                    double area, bool anchor) {
  using std::exp;
  const double canopy_top = e->max_environment_height();
  auto integrand = [&](S z) -> S {
    double zv = as_double(z);
    double lv = e->get_environment_at_height(zv, canopy_top);
    double ld = e->get_environment_deriv_at_height(zv);
    S light;
    if (anchor) {
      light = S(lv) + S(std::isfinite(ld)?ld:0.0) * (z - S(zv));   // frozen value + z-lin
      light = light + (resident_light_anchored<S>(zv, lv, pd, eta, sh, sC, area) - S(lv));
    } else {
      // genuine recon value with z ACTIVE -> carries BOTH the trait feedback and the
      // focal-height->light channel (recon's own z-derivative). Self-consistent target
      // for the R1 finite-difference (forward-AD == FD over this exact function).
      S comp = recon_competition<S>(z, pd, eta, sh, sC) * S(1.0 / area);
      light = exp(-comp);
    }
    return plant::ff16_assimilation_leaf<S>(pd.a_p1, pd.a_p2, light) *
           plant::ff16_canopy_q<S>(eta, z / height, z);
  };
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S assim = area_leaf * integ->integrate_ad<S>(integrand, S(0.0), height);
  return plant::ff16_net_from_components(pd, height, area_leaf, assim);
}

// Replay ONE cohort (index i) over the frozen schedule, returning its FINAL
// FF16LifeState<S> (the 5 demographic states + survival-weighted offspring
// accumulator at the end of the run, i.e. at the final census). Establishment is
// differentiated via the seedling net production -> mortality_0 = -log(pr_estab);
// h0 carries its own d/d(trait) via the IFT injection done by the caller. This is
// the single per-cohort replay that the generic reduction engine sums over: every
// emergent metric is a weighted reduction over these final cohort states, so the
// replay is recorded ONCE and reused for every metric (one tape, one reverse sweep
// per metric).
template <typename S>
plant::FF16LifeState<S> replay_cohort_final(const plant::FF16ProdPars<S>& pd,
                                            const Frozen& F, std::size_t i, S h0) {
  using std::exp; using std::log;
  const std::size_t b = (std::size_t)F.birth[i];
  const double ppsab = F.ppsab[i];
  const plant::FF16_Environment* eb = (b > 0) ? &F.eh[b - 1][5] : &F.eh[0][0];
  S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
  S net0 = deep_net<S>(pd, F.integ, F.eta, eb, h0);
  S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay[i]);
  S mort0 = -log(pr_estab);
  auto deriv = [&](const plant::FF16LifeState<S>& s, std::size_t n, int stage)
      -> plant::FF16LifeState<S> {
    const plant::FF16_Environment* e =
      (stage==0)?((n>0)?&F.eh[n-1][5]:&F.eh[0][0]):&F.eh[n][stage-1];
    S net = deep_net<S>(pd, F.integ, F.eta, e, s.demog.height);
    S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, s.demog.height);
    plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, s.demog.height, area_leaf, net, true);
    S off_dt = r.fecundity_dt * exp(-s.demog.mortality) * S(F.ppsurv(n, stage) / ppsab);
    return plant::FF16LifeState<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt,
      r.fecundity_dt, r.area_heartwood_dt, r.mass_heartwood_dt}, off_dt};
  };
  auto axpy = [](const plant::FF16LifeState<S>& a, double c, const plant::FF16LifeState<S>& k)
      -> plant::FF16LifeState<S> {
    return plant::FF16LifeState<S>{plant::FF16State<S>{
      a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
      a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
      a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.offspring+c*k.offspring};
  };
  plant::FF16LifeState<S> y{plant::FF16State<S>{h0, mort0, S(0), S(0), S(0)}, S(0)};
  return plant::ff16_cashkarp_replay(y, F.step_h, b, deriv, axpy);
}

// Deep-crown dheight/dt at `height` reading the frozen env `e` (= the resident
// growth rate g of a plant of that height). area_leaf is derived from height.
template <typename S>
S deep_height_dt(const plant::FF16ProdPars<S>& pd, const plant::quadrature::QK* integ,
                 double eta, const plant::FF16_Environment* e, S height) {
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S net = deep_net<S>(pd, integ, eta, e, height);
  return plant::ff16_height_dt_from_net(pd, height, area_leaf, net);
}

// Resident counterpart of deep_height_dt: dheight/dt reading the value-anchored
// resident light (used by the census log_density rate / its backward-FD gradient).
template <typename S>
S deep_height_dt_resident(const plant::FF16ProdPars<S>& pd,
                          const plant::quadrature::QK* integ, double eta,
                          const plant::FF16_Environment* e, S height,
                          const std::vector<double>& sh, const std::vector<double>& sC,
                          double area, bool anchor) {
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S net = deep_net_resident<S>(pd, integ, eta, e, height, sh, sC, area, anchor);
  return plant::ff16_height_dt_from_net(pd, height, area_leaf, net);
}

// Dispatch helpers: net production / height rate at replay (step n, RK stage), via
// either the frozen env (mutant path) or the value-anchored resident reconstruction.
// The (n, stage) -> env mapping mirrors the replay deriv rule; the matching stand is
// looked up at the SAME (n, stage). Bit-identical to the direct deep_net call when
// F.resident is false.
template <typename S>
const plant::FF16_Environment* stage_env(const Frozen& F, std::size_t n, int stage) {
  return (stage == 0) ? ((n > 0) ? &F.eh[n-1][5] : &F.eh[0][0]) : &F.eh[n][stage-1];
}
template <typename S>
S net_at(const plant::FF16ProdPars<S>& pd, const Frozen& F, std::size_t n, int stage, S h) {
  const plant::FF16_Environment* e = stage_env<S>(F, n, stage);
  if (!F.resident) return deep_net<S>(pd, F.integ, F.eta, e, h);
  const auto& sh = (stage == 0) ? ((n > 0) ? F.st_h[n-1][5] : F.st_h[0][0]) : F.st_h[n][stage-1];
  const auto& sC = (stage == 0) ? ((n > 0) ? F.st_C[n-1][5] : F.st_C[0][0]) : F.st_C[n][stage-1];
  return deep_net_resident<S>(pd, F.integ, F.eta, e, h, sh, sC, F.area, F.anchor);
}
template <typename S>
S height_dt_at(const plant::FF16ProdPars<S>& pd, const Frozen& F, std::size_t n, int stage, S h) {
  const plant::FF16_Environment* e = stage_env<S>(F, n, stage);
  if (!F.resident) return deep_height_dt<S>(pd, F.integ, F.eta, e, h);
  const auto& sh = (stage == 0) ? ((n > 0) ? F.st_h[n-1][5] : F.st_h[0][0]) : F.st_h[n][stage-1];
  const auto& sC = (stage == 0) ? ((n > 0) ? F.st_C[n-1][5] : F.st_C[0][0]) : F.st_C[n][stage-1];
  return deep_height_dt_resident<S>(pd, F.integ, F.eta, e, h, sh, sC, F.area, F.anchor);
}

// Census-state replay: the 5 demographic states + log_density (the cohort's
// size-distribution number density along its characteristic). log_density evolves
// by the method-of-characteristics rate the SCM uses (Node::compute_rates):
//   d(log_density)/dt = - growth_rate_gradient - mortality_dt,
// growth_rate_gradient = d(height_dt)/d(height). We reproduce the SCM's EXACT
// scheme for that gradient (Node::growth_rate_gradient default): a BACKWARD finite
// difference with ABSOLUTE step node_gradient_eps = 1e-6, reusing the already-known
// height_dt as fx -- g' = (height_dt(h) - height_dt(h - eps)) / eps. Replicating
// the SCM's own g' (rather than an exact derivative) is what makes the replayed
// census density reproduce the SCM's stored density; the trait derivative of that
// FD expression is taken exactly by AD, consistent with how the SCM formed it.
template <typename S> struct CensusState { plant::FF16State<S> demog; S log_density; };

template <typename S>
CensusState<S> replay_cohort_census(const plant::FF16ProdPars<S>& pd, const Frozen& F,
                                    std::size_t i, S h0, double birth_rate) {
  using std::log;
  const std::size_t b = (std::size_t)F.birth[i];
  const plant::FF16_Environment* eb = (b > 0) ? &F.eh[b - 1][5] : &F.eh[0][0];
  const double GEPS = 1e-6;                          // Control::node_gradient_eps
  // Establishment + seedling growth rate g0 in the birth env -> log_density_0.
  S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
  S net0 = deep_net<S>(pd, F.integ, F.eta, eb, h0);
  S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay[i]);
  S mort0 = -log(pr_estab);
  S g0 = deep_height_dt<S>(pd, F.integ, F.eta, eb, h0);
  S logd0 = log(S(birth_rate) * pr_estab / g0);      // = log(birth_rate*pr_estab/g0)

  auto deriv = [&](const CensusState<S>& s, std::size_t n, int stage) -> CensusState<S> {
    const plant::FF16_Environment* e =
      (stage==0)?((n>0)?&F.eh[n-1][5]:&F.eh[0][0]):&F.eh[n][stage-1];
    const S h = s.demog.height;
    S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h);
    S net = deep_net<S>(pd, F.integ, F.eta, e, h);
    plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, h, area_leaf, net, true);
    // g' = backward FD (abs step GEPS) of height_dt, reusing r.height_dt as fx.
    S g_back = deep_height_dt<S>(pd, F.integ, F.eta, e, h - S(GEPS));
    S gprime = (r.height_dt - g_back) / S(GEPS);
    S log_density_dt = -gprime - r.mortality_dt;
    return CensusState<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
      r.area_heartwood_dt, r.mass_heartwood_dt}, log_density_dt};
  };
  auto axpy = [](const CensusState<S>& a, double c, const CensusState<S>& k) -> CensusState<S> {
    return CensusState<S>{plant::FF16State<S>{
      a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
      a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
      a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.log_density+c*k.log_density};
  };
  CensusState<S> y{plant::FF16State<S>{h0, mort0, S(0), S(0), S(0)}, logd0};
  return plant::ff16_cashkarp_replay(y, F.step_h, b, deriv, axpy);
}


// Unified cohort state for the generic engine: the 5 demographic states + the
// survival-weighted offspring accumulator (Lagrangian metrics like
// offspring_production) + log_density (census metrics like LAI / biomass / size
// moments). One replay carries BOTH so every metric shares ONE recorded tape. The
// offspring and log_density accumulators are independent fields, so this is
// bit-identical to the dedicated replays for each (the offspring gate is unchanged).
template <typename S>
struct FullState { plant::FF16State<S> demog; S offspring; S log_density; };

template <typename S>
FullState<S> replay_cohort_full(const plant::FF16ProdPars<S>& pd, const Frozen& F,
                                std::size_t i, S h0, double birth_rate) {
  using std::exp; using std::log;
  const std::size_t b = (std::size_t)F.birth[i];
  const double ppsab = F.ppsab[i];
  const double GEPS = 1e-6;                              // Control::node_gradient_eps
  S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
  S net0 = net_at<S>(pd, F, b, 0, h0);
  S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay[i]);
  S mort0 = -log(pr_estab);
  S g0 = height_dt_at<S>(pd, F, b, 0, h0);
  S logd0 = log(S(birth_rate) * pr_estab / g0);

  auto deriv = [&](const FullState<S>& s, std::size_t n, int stage) -> FullState<S> {
    const S h = s.demog.height;
    S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h);
    S net = net_at<S>(pd, F, n, stage, h);
    plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, h, area_leaf, net, true);
    S off_dt = r.fecundity_dt * exp(-s.demog.mortality) * S(F.ppsurv(n, stage) / ppsab);
    // log_density rate: -(growth_rate_gradient + mortality_dt); g' replicates the
    // SCM's backward FD (abs step GEPS), reusing r.height_dt as fx.
    S g_back = height_dt_at<S>(pd, F, n, stage, h - S(GEPS));
    S gprime = (r.height_dt - g_back) / S(GEPS);
    S log_density_dt = -gprime - r.mortality_dt;
    return FullState<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
      r.area_heartwood_dt, r.mass_heartwood_dt}, off_dt, log_density_dt};
  };
  auto axpy = [](const FullState<S>& a, double c, const FullState<S>& k) -> FullState<S> {
    return FullState<S>{plant::FF16State<S>{
      a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
      a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
      a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.offspring+c*k.offspring,
      a.log_density+c*k.log_density};
  };
  FullState<S> y{plant::FF16State<S>{h0, mort0, S(0), S(0), S(0)}, S(0), logd0};
  return plant::ff16_cashkarp_replay(y, F.step_h, b, deriv, axpy);
}

// The new_node (pending seed) census contribution, born into the FINAL-time env
// (the SCM's compute_competition tail term). Establishment + seedling growth rate
// are evaluated in the last cached stage env; returns (density_new, h0).
template <typename S>
void new_node_census(const plant::FF16ProdPars<S>& pd, const Frozen& F, S h0,
                     double birth_rate, S& density_new) {
  using std::log; using std::exp;
  // Final-time env (eh.back()[5]) + matching final-stage stand: net_at(N, 0, .) maps
  // stage 0 of the would-be next step N=eh.size() to eh[N-1][5] = eh.back()[5].
  const std::size_t N = F.eh.size();
  S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
  S net0 = net_at<S>(pd, F, N, 0, h0);
  // The seed's establishment uses the run-end decay factor (last birth time slot).
  S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0,
                                                        F.decay.back());
  S g0 = height_dt_at<S>(pd, F, N, 0, h0);
  density_new = S(birth_rate) * pr_estab / g0;           // = exp(log_density_new)
}

// ---------------------------------------------------------------------------
// Generic (w, f) reduction engine (#472 scope B, build-order step 1). A stand
// metric is a weighted reduction over the replayed cohort final states,
//   metric = sum_i w_i * f(state_i),
// where w_i is a frozen per-cohort weight (from pass 1) and f maps a cohort's
// final FF16LifeState to a scalar contribution. offspring_production, LAI, biomass
// and size-distribution moments are all symmetric instantiations -- NONE privileged.
// The engine records ONE forward replay of every cohort onto a single tape, then
// takes one cheap reverse (adjoint) sweep PER metric, so M metrics cost
// replay + M*sweep, not M replays. This is the calibration-facing core: a
// metrics x traits Jacobian + the metric values, out of one resident baseline.
// ---------------------------------------------------------------------------

// Assemble the Frozen harvest (resident schedule + per-RK-stage env + per-cohort
// survival/weights/decay) from the R-side arrays. Shared by every entry point.
Frozen build_frozen(const plant::FF16_Strategy& s, Rcpp::List eh_list,
                    const std::vector<double>& sh, std::vector<int> birth,
                    Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
                    std::vector<double> tw) {
  Frozen F; F.eta = s.pars.eta; F.h0 = s.initial_height(); F.birth = birth;
  F.ppsab = ppsab; F.tw = tw; F.ppsurv = ppsurv; F.integ = &s.function_integrator;
  F.a_d0 = s.pars.a_d0;
  const std::size_t N = eh_list.size(); F.eh.resize(N); F.step_h.resize(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n]; for(R_xlen_t k=0;k<st.size();++k) F.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));}
  for (std::size_t n=0;n<N;++n) F.step_h[n]=sh[n+1]-sh[n];
  // Per-cohort establishment decay exp(-recruitment_decay * birth_time).
  F.decay.resize(birth.size());
  for (std::size_t i=0;i<birth.size();++i)
    F.decay[i] = std::exp(-s.pars.recruitment_decay * sh[(std::size_t)birth[i]]);
  return F;
}

// FULL native harvest (#472 scope B, refactor+optimize phase): build the ENTIRE
// Frozen harvest from the live SCM's Patch -- the per-RK-stage env (faithful C++
// copy, no Rcpp::as<>), the step schedule, and the per-cohort birth steps /
// trapezoid weights / per-RK-stage patch survival / survival-at-birth -- replacing
// the R-side ff16_harvest (whose pr_survival loop is O(stand size) because each
// scm$patch access rebuilds the whole RcppR6 patch; see R/emergent_gradient.R:70).
// Bit-identical to ff16_harvest + build_frozen (the same arithmetic, native data).
Frozen build_frozen_scm(
    const plant::FF16_Strategy& s,
    const plant::Patch<plant::FF16_Strategy, plant::FF16_Environment>& patch,
    std::size_t species, double birth_rate) {
  Frozen F; F.eta = s.pars.eta; F.h0 = s.initial_height();
  F.integ = &s.function_integrator; F.a_d0 = s.pars.a_d0;
  const auto& EH = patch.environment_history;
  const auto& sh = patch.step_history;
  const std::size_t N = EH.size();
  F.eh = EH;                                   // faithful copy (no Rcpp::as<>)
  F.step_h.resize(N);
  for (std::size_t n = 0; n < N; ++n) F.step_h[n] = sh[n + 1] - sh[n];

  // Birth steps / trapezoid weights / per-RK-stage + at-birth survival: the shared
  // native offspring harvest (identical computation across FF16 / TF24 / TF24f).
  auto w = plant::gradient::offspring_weights(patch, species, s.pars.S_D, birth_rate);
  F.birth = w.birth; F.ppsab = w.ppsab; F.tw = w.tw;
  F.birth_rate0 = birth_rate;                  // fixed harvest br for the R0 unit weights
  const std::size_t nC = F.birth.size();
  F.decay.resize(nC);
  for (std::size_t i = 0; i < nC; ++i)
    F.decay[i] = std::exp(-s.pars.recruitment_decay * sh[(std::size_t)F.birth[i]]);
  Rcpp::NumericMatrix ppsurv((int)N, 6);
  for (std::size_t k = 0; k < N; ++k)
    for (int st = 0; st < 6; ++st) ppsurv((int)k, st) = w.ppsurv[k * 6 + (std::size_t)st];
  F.ppsurv = ppsurv;
  return F;
}

// Populate the RESIDENT per-RK-stage stand harvest on F from the R-side nested lists
// stand_height_stage_history / stand_competition_stage_history ([step][stage][cohort]).
// Each (step, stage) cohort vector is sorted DESCENDING by height and the frozen
// weight C_i = ce_i / area_leaf(theta0, h_i) is precomputed (so C_i*area_leaf(theta)
// reproduces ce_i at theta0 and flows the allometric trait). a_l1_0/a_l2_0 are the
// base (resident) allometry parameters.
void attach_resident_stand(Frozen& F, Rcpp::List sh_h_list, Rcpp::List sh_c_list,
                           double a_l1_0, double a_l2_0, double area, bool anchor) {
  F.resident = true; F.area = area; F.anchor = anchor;
  const std::size_t N = sh_h_list.size();
  F.st_h.resize(N); F.st_C.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List hn = sh_h_list[n], cn = sh_c_list[n];
    const std::size_t ns = hn.size();
    F.st_h[n].resize(ns); F.st_C[n].resize(ns);
    for (std::size_t st = 0; st < ns; ++st) {
      std::vector<double> hv = Rcpp::as<std::vector<double>>(hn[st]);
      std::vector<double> cv = Rcpp::as<std::vector<double>>(cn[st]);
      // descending-height order (matches Species::compute_competition's trapezium walk)
      std::vector<std::size_t> ord(hv.size());
      for (std::size_t i = 0; i < hv.size(); ++i) ord[i] = i;
      std::sort(ord.begin(), ord.end(),
                [&](std::size_t a, std::size_t b){ return hv[a] > hv[b]; });
      F.st_h[n][st].resize(hv.size()); F.st_C[n][st].resize(hv.size());
      for (std::size_t k = 0; k < ord.size(); ++k) {
        const double hh = hv[ord[k]];
        const double al = std::pow(hh / a_l1_0, 1.0 / a_l2_0);   // area_leaf(theta0, hh)
        F.st_h[n][st][k] = hh;
        F.st_C[n][st][k] = (al > 0.0) ? cv[ord[k]] / al : 0.0;
      }
    }
  }
}

// Map requested trait names to FF16ProdPars field indices (errors on unknown).
std::vector<std::size_t> resolve_traits(const std::vector<std::string>& traits) {
  auto names = field_names();
  std::vector<std::size_t> idx;
  for (auto& t : traits) {
    auto it = std::find(names.begin(), names.end(), t);
    if (it == names.end()) Rcpp::stop("unknown FF16 trait: " + t);
    idx.push_back(std::distance(names.begin(), it));
  }
  return idx;
}

// Seed an active h0 carrying the IFT injection d(h0)/d(theta_k) for each registered
// trait (zero for traits not in mass_live), shared by every reverse-sweep entry point.
ad_t inject_h0(double h0v, const std::vector<ad_t*>& fp,
               const std::vector<std::size_t>& idx, const std::vector<double>& dh0) {
  ad_t h0 = h0v;
  for (std::size_t k = 0; k < idx.size(); ++k)
    h0 = h0 + ad_t(dh0[k]) * (*fp[idx[k]] - ad_t(xad::value(*fp[idx[k]])));
  return h0;
}

// d(offspring_production)/d(trait) accumulated PER COHORT. offspring_production is the
// weighted sum sum_i tw_i * offspring_i, and under the frozen resident the cohorts are
// independent (each reads the harvested env as a double constant), so the gradient of
// the sum is the sum of the gradients. Taping one cohort at a time keeps each tape tiny
// (~1/nC the nodes) -- bit-identical to one monolithic tape but ~25x faster, because the
// single giant tape (millions of nodes) thrashed cache. Fills grad += and value +=.
void accumulate_offspring_gradient(const plant::FF16ProdPars<double>& pd, const Frozen& F,
    double h0v, const std::vector<double>& dh0, const std::vector<std::size_t>& idx,
    std::vector<double>& grad, double& value) {
  const std::size_t nC = F.birth.size(), nT = idx.size();
  grad.assign(nT, 0.0); value = 0.0;
  for (std::size_t i = 0; i < nC; ++i) {
    ad::tape_type tape;
    auto pa = lift<ad_t>(pd);
    auto fp = field_ptrs<ad_t>(pa);
    for (auto j : idx) tape.registerInput(*fp[j]);
    tape.newRecording();
    ad_t h0 = inject_h0(h0v, fp, idx, dh0);
    ad_t Ji = ad_t(F.tw[i]) * replay_cohort_final<ad_t>(pa, F, i, h0).offspring;
    tape.registerOutput(Ji); xad::derivative(Ji) = 1.0; tape.computeAdjoints();
    for (std::size_t k = 0; k < nT; ++k) grad[k] += xad::derivative(*fp[idx[k]]);
    value += as_double(Ji);
  }
}

// d(height_0)/d(trait_k) by the implicit function theorem at the height_seed root
// (mass_live_given_height(h0) == omega), in a SCOPED reverse sweep of mass_live at
// h0 (the #539 seedling-size pattern). Returns the per-requested-trait dh0.
std::vector<double> compute_dh0(const plant::FF16ProdPars<double>& pd, double h0v,
                                const std::vector<std::size_t>& idx) {
  std::vector<double> dh0(idx.size(), 0.0);
  ad::tape_type tape0;
  auto pm = lift<ad_t>(pd);
  auto fm = field_ptrs<ad_t>(pm);
  ad_t hin = h0v;
  for (auto i : idx) tape0.registerInput(*fm[i]);
  tape0.registerInput(hin);
  tape0.newRecording();
  ad_t m = plant::ff16_mass_live_given_height<ad_t>(pm, hin);
  tape0.registerOutput(m); xad::derivative(m) = 1.0; tape0.computeAdjoints();
  const double dm_dh = xad::derivative(hin);
  auto names_o = field_names();
  for (std::size_t k = 0; k < idx.size(); ++k) {
    double dm_dtheta = xad::derivative(*fm[idx[k]]);
    double dF_dtheta = dm_dtheta - (names_o[idx[k]] == "omega" ? 1.0 : 0.0);
    dh0[k] = (dm_dh != 0.0) ? -dF_dtheta / dm_dh : 0.0;
  }
  return dh0;
}

// Templated stand-metric assembly: replay every cohort, then reduce each requested
// metric to a scalar. Shared by (a) the FROZEN reverse sweep (S = ad_t, F.resident
// false), (b) the RESIDENT forward sweep (S = fwd_t, F.resident true -- one trait
// direction, NO tape, so the O(stand) recon per evaluation never builds a giant
// adjoint graph), and (c) plain double value reconstruction. Census metrics are the
// height-trapezium of density_i * psi(state_i) over the cohorts (descending height)
// + the pending-seed tail; offspring_production is the frozen-weighted offspring sum.
template <typename S>
std::vector<S> assemble_metrics(const plant::FF16ProdPars<S>& pa, const Frozen& F,
    const std::vector<std::string>& metrics, double birth_rate, double kI, S h0) {
  using std::exp;
  const std::size_t nC = F.birth.size(), M = metrics.size();
  std::vector<FullState<S>> finals; finals.reserve(nC);
  for (std::size_t i = 0; i < nC; ++i)
    finals.push_back(replay_cohort_full<S>(pa, F, i, h0, birth_rate));
  S dens_new; new_node_census<S>(pa, F, h0, birth_rate, dens_new);

  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
    return as_double(finals[a].demog.height) > as_double(finals[b].demog.height); });

  auto census_reduce = [&](auto psi) -> S {
    return plant::gradient::census_trapezium<S>(nC, ord, h0, dens_new,
      [&](std::size_t i) -> S { return finals[i].demog.height; },
      [&](std::size_t i) -> S { return exp(finals[i].log_density); },
      [&](std::size_t i) -> S { return finals[i].demog.mass_heartwood; }, psi);
  };

  std::vector<S> J(M);
  for (std::size_t m = 0; m < M; ++m) {
    const std::string& nm = metrics[m];
    if (nm == "offspring_production") {
      S acc = S(0.0);
      for (std::size_t i = 0; i < nC; ++i) acc = acc + S(F.tw[i]) * finals[i].offspring;
      J[m] = acc;
    } else if (nm == "LAI") {
      J[m] = census_reduce([&](S h, S dens, S) -> S {
        return dens * S(kI) * plant::ff16_area_leaf(pa.a_l1, pa.a_l2, h); });
    } else if (nm == "biomass") {
      J[m] = census_reduce([&](S h, S dens, S mhw) -> S {
        return dens * (plant::ff16_mass_live_given_height(pa, h) + mhw); });
    } else { // size_moment
      J[m] = census_reduce([&](S h, S dens, S) -> S { return dens * h; });
    }
  }
  return J;
}

// ===========================================================================
// COUPLED resident replay (#472 scope B, R0-R1 -- the course-corrected build).
// All alive cohorts are stepped TOGETHER over the frozen schedule; at each RK
// stage the canopy light is reconstructed from the CURRENT active stand (active
// heights h_i, active densities exp(log_density_i), active area_leaf), filling an
// odelia basic_interpolator<S> at the env spline's FROZEN knot x-positions. So
// every trait that moves a height or a density re-shades the whole stand -- the
// genuine resident feedback, not the leaf-area-only frozen-geometry graft. Built
// ONCE per stage and shared across cohorts (O(N) per knot, not O(N^2)).
// ===========================================================================

// Competition at a frozen knot z over the active stand: see
// plant::gradient::canopy_comp_at (shared Yokozawa trapezium, used by both the
// single- and multi-species coupled reconstructions here and by TF24f).

// Deep-crown net at `height` reading the ACTIVE light interpolator (the coupled
// resident canopy). Mirrors deep_net but light(z) = interp(z) [+ analytic z-slope
// for the focal-height channel], active in the interpolator's knot values. Above the
// frozen cap light is the open value 1; the #253 floor (max(0,.)) and its zero
// derivative are replicated so it matches FF16_Environment::get_environment_at_height
// /get_environment_deriv_at_height at theta0 (where interp == the frozen env spline).
template <typename S>
S deep_net_coupled(const plant::FF16ProdPars<S>& pd, const plant::quadrature::QK* integ,
                   double eta, const odelia::interpolator::basic_interpolator<S>& interp,
                   double cap, S height) {
  auto integrand = [&](S z) -> S {
    double zv = as_double(z);
    S lv, ld;
    if (zv > cap) { lv = S(1.0); ld = S(0.0); }
    else {
      lv = interp(zv);
      if (as_double(lv) > 0.0) { ld = interp.deriv(zv); }
      else { lv = S(0.0); ld = S(0.0); }
    }
    S light = lv + ld * (z - S(zv));
    return plant::ff16_assimilation_leaf<S>(pd.a_p1, pd.a_p2, light) *
           plant::ff16_canopy_q<S>(eta, z / height, z);
  };
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S assim = area_leaf * integ->integrate_ad<S>(integrand, S(0.0), height);
  return plant::ff16_net_from_components(pd, height, area_leaf, assim);
}

template <typename S>
S deep_height_dt_coupled(const plant::FF16ProdPars<S>& pd,
                         const plant::quadrature::QK* integ, double eta,
                         const odelia::interpolator::basic_interpolator<S>& interp,
                         double cap, S height) {
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S net = deep_net_coupled<S>(pd, integ, eta, interp, cap, height);
  return plant::ff16_height_dt_from_net(pd, height, area_leaf, net);
}

// Single Cash-Karp RKCK step over an arbitrary state (the whole-stand vector), one
// frozen step of size h at replay step rn. Same tableau as ff16_cashkarp_replay,
// peeled to a single step so the coupled driver can introduce cohorts between steps.
template <typename State, typename DerivFn, typename AxpyFn>
State rkck_one_step(State y, double h, std::size_t rn, DerivFn&& deriv, AxpyFn&& axpy) {
  const double b21 = 1.0 / 5.0;
  const double b3[2] = {3.0 / 40.0, 9.0 / 40.0};
  const double b4[3] = {0.3, -0.9, 1.2};
  const double b5[4] = {-11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0};
  const double b6[5] = {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0,
                        44275.0 / 110592.0, 253.0 / 4096.0};
  const double c1 = 37.0 / 378.0, c3 = 250.0 / 621.0,
               c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0;
  const State k1 = deriv(y, rn, 0);
  const State k2 = deriv(axpy(y, b21 * h, k1), rn, 1);
  State y3 = axpy(y, h * b3[0], k1); y3 = axpy(y3, h * b3[1], k2);
  const State k3 = deriv(y3, rn, 2);
  State y4 = axpy(y, h * b4[0], k1);
  y4 = axpy(y4, h * b4[1], k2); y4 = axpy(y4, h * b4[2], k3);
  const State k4 = deriv(y4, rn, 3);
  State y5 = axpy(y, h * b5[0], k1);
  y5 = axpy(y5, h * b5[1], k2); y5 = axpy(y5, h * b5[2], k3); y5 = axpy(y5, h * b5[3], k4);
  const State k5 = deriv(y5, rn, 4);
  State y6 = axpy(y, h * b6[0], k1);
  y6 = axpy(y6, h * b6[1], k2); y6 = axpy(y6, h * b6[2], k3);
  y6 = axpy(y6, h * b6[3], k4); y6 = axpy(y6, h * b6[4], k5);
  const State k6 = deriv(y6, rn, 5);
  return axpy(axpy(axpy(axpy(y, h * c1, k1), h * c3, k3), h * c4, k4), h * c6, k6);
}

// The coupled whole-stand replay + emergent-metric reduction. Re-evolves all
// cohorts together over the frozen schedule, reconstructing the active canopy each
// RK stage, and reduces the final census to the requested metrics. Cohort initial
// conditions use the FROZEN harvested birth env (exact at theta0; the birth-env
// feedback is a small deferred channel) -- only the GROWTH-phase env is active.
// If env_err != nullptr (double only), accumulates max|reconstructed knot light -
// SCM knot light| across all stages (the R0 coupled-drift gauge).
template <typename S>
std::vector<S> assemble_metrics_coupled(const plant::FF16ProdPars<S>& pd, const Frozen& F,
    const std::vector<std::string>& metrics, S birth_rate, S h0,
    double* env_err = nullptr, double* env_err_z = nullptr) {
  using std::exp; using std::log;
  using interp_t = odelia::interpolator::basic_interpolator<S>;
  const std::size_t nC = F.birth.size(), N = F.eh.size();
  const double GEPS = 1e-6;                          // Control::node_gradient_eps

  // FullState carries a per-cohort survival-weighted OFFSPRING accumulator alongside the
  // {5 demog, log_density} census state. Offspring is fully decoupled from the demog /
  // log_density rates, so carrying it leaves the census + density trajectories bit-
  // identical; it is reduced only for the net_reproduction_ratio metric (the demographic-
  // equilibrium / mutant-fitness axis -- birth_rate enters offspring ONLY through the
  // active canopy, the pure density feedback).
  std::vector<FullState<S>> stand(nC);               // all cohorts, zero until birth
  std::vector<char> alive(nC, 0);

  // (rn, stage) -> the cached env (en, es) it reads, mirroring stage_env().
  auto env_idx = [&](std::size_t rn, int rs, std::size_t& en, int& es) {
    if (rs == 0) { if (rn > 0) { en = rn - 1; es = 5; } else { en = 0; es = 0; } }
    else { en = rn; es = rs - 1; }
  };

  // Build the active light interpolator at (en, es)'s frozen knots from the current
  // alive stand + the frozen boundary node. Shared by all cohorts in this stage.
  auto build_interp = [&](const std::vector<FullState<S>>& y,
                          std::size_t en, int es) -> interp_t {
    const std::vector<double>& kx = F.knot_x[en][es];
    std::vector<S> hv, gv;
    hv.reserve(nC + 1); gv.reserve(nC + 1);
    for (std::size_t i = 0; i < nC; ++i) if (alive[i]) {
      S hi = y[i].demog.height;
      S dens = exp(y[i].log_density);
      S al = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, hi);
      hv.push_back(hi); gv.push_back(dens * S(F.kI) * al);
    }
    hv.push_back(S(F.nn_h[en][es])); gv.push_back(S(F.nn_c[en][es]));  // boundary tail
    std::vector<std::size_t> ord(hv.size());
    for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(),
              [&](std::size_t a, std::size_t b){ return as_double(hv[a]) > as_double(hv[b]); });
    std::vector<S> hs(hv.size()), gs(hv.size());
    for (std::size_t k = 0; k < ord.size(); ++k) { hs[k] = hv[ord[k]]; gs[k] = gv[ord[k]]; }
    std::vector<S> ly(kx.size());
    for (std::size_t k = 0; k < kx.size(); ++k)
      ly[k] = exp(-plant::gradient::canopy_comp_at<S>(kx[k], hs, gs, F.eta) * S(1.0 / F.area));
    if (env_err) {
      const std::vector<double>& y0 = F.knot_y0[en][es];
      for (std::size_t k = 0; k < ly.size() && k < y0.size(); ++k) {
        double e = std::abs(as_double(ly[k]) - y0[k]);
        if (e > *env_err) { *env_err = e; if (env_err_z) *env_err_z = kx[k]; }
      }
    }
    interp_t interp; interp.init(kx, ly);
    return interp;
  };

  // Whole-stand derivative at replay step rn, RK stage rs (alive cohorts only).
  auto deriv = [&](const std::vector<FullState<S>>& y, std::size_t rn, int rs)
      -> std::vector<FullState<S>> {
    std::size_t en; int es; env_idx(rn, rs, en, es);
    interp_t interp = build_interp(y, en, es);
    const double cap = F.knot_x[en][es].back();
    std::vector<FullState<S>> dy(nC);
    for (std::size_t i = 0; i < nC; ++i) {
      if (!alive[i]) {
        dy[i] = FullState<S>{plant::FF16State<S>{S(0),S(0),S(0),S(0),S(0)}, S(0), S(0)};
        continue;
      }
      const S h = y[i].demog.height;
      S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h);
      S net = deep_net_coupled<S>(pd, F.integ, F.eta, interp, cap, h);
      plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, h, area_leaf, net, true);
      // log_density rate: -(g' + mortality_dt); g' = backward FD (abs step GEPS) of
      // height_dt, reusing r.height_dt as fx (the SCM's own scheme).
      S g_back = deep_height_dt_coupled<S>(pd, F.integ, F.eta, interp, cap, h - S(GEPS));
      S gprime = (r.height_dt - g_back) / S(GEPS);
      S log_density_dt = -gprime - r.mortality_dt;
      // Survival-weighted offspring rate (same as the frozen offspring tape): the focal
      // individual's fecundity discounted by its own survival exp(-mortality) and the
      // per-stage patch survival ppsurv/ppsab. fecundity reads the ACTIVE canopy, so
      // offspring carries the birth_rate density feedback (and nothing else).
      S off_dt = r.fecundity_dt * exp(-y[i].demog.mortality)
                 * S(F.ppsurv(rn, rs) / F.ppsab[i]);
      dy[i] = FullState<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, off_dt, log_density_dt};
    }
    return dy;
  };
  auto axpy = [&](const std::vector<FullState<S>>& a, double c,
                  const std::vector<FullState<S>>& k) -> std::vector<FullState<S>> {
    std::vector<FullState<S>> r(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
      r[i] = FullState<S>{plant::FF16State<S>{
        a[i].demog.height + c * k[i].demog.height,
        a[i].demog.mortality + c * k[i].demog.mortality,
        a[i].demog.fecundity + c * k[i].demog.fecundity,
        a[i].demog.area_heartwood + c * k[i].demog.area_heartwood,
        a[i].demog.mass_heartwood + c * k[i].demog.mass_heartwood},
        a[i].offspring + c * k[i].offspring,
        a[i].log_density + c * k[i].log_density};
    return r;
  };

  // March the frozen schedule: introduce cohorts born at step rn, then step all.
  for (std::size_t rn = 0; rn < N; ++rn) {
    // Birth canopy at the start of step rn (stage 0 env env_idx(rn,0)). ACTIVE path:
    // reconstruct it from the already-alive cohorts (born < rn) + the boundary node
    // BEFORE adding the cohorts born at rn -- non-circular (cohort i is not yet
    // alive), exact at theta0 (the alive set here reproduces the harvested stand at
    // (rn-1,5)). FROZEN path: read the harvested env F.eh as a constant.
    bool any_born = false;
    for (std::size_t i = 0; i < nC; ++i)
      if ((std::size_t)F.birth[i] == rn) { any_born = true; break; }
    if (any_born) {
      std::size_t en; int es; env_idx(rn, 0, en, es);
      const plant::FF16_Environment* eb = stage_env<S>(F, rn, 0);
      interp_t binterp;
      double bcap = 0.0;
      const bool active_be = F.coupled_active_birthenv;
      if (active_be) { binterp = build_interp(stand, en, es); bcap = F.knot_x[en][es].back(); }
      for (std::size_t i = 0; i < nC; ++i) {
        if ((std::size_t)F.birth[i] != rn) continue;
        S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
        S net0, g0;
        if (active_be) {
          net0 = deep_net_coupled<S>(pd, F.integ, F.eta, binterp, bcap, h0);
          g0   = deep_height_dt_coupled<S>(pd, F.integ, F.eta, binterp, bcap, h0);
        } else {
          net0 = deep_net<S>(pd, F.integ, F.eta, eb, h0);
          g0   = deep_height_dt<S>(pd, F.integ, F.eta, eb, h0);
        }
        S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay[i]);
        S mort0 = -log(pr_estab);
        S logd0 = log(birth_rate * pr_estab / g0);
        stand[i] = FullState<S>{plant::FF16State<S>{h0, mort0, S(0), S(0), S(0)}, S(0), logd0};
        alive[i] = 1;
      }
    }
    stand = rkck_one_step(stand, F.step_h[rn], rn, deriv, axpy);
  }

  // Cohorts whose birth step lands on the final boundary (birth >= N) are introduced
  // at the last node time but never stepped: mirror replay_cohort_full (b >= N runs
  // zero Cash-Karp steps, so the cohort sits at h0). Establish them at h0 in the frozen
  // final birth env -- otherwise they keep the zero-initialised height and the a_l2
  // channel of area_leaf = (h/a_l1)^(1/a_l2) is differentiated at h=0 (0*log(0)=NaN).
  for (std::size_t i = 0; i < nC; ++i) {
    if (alive[i]) continue;
    S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
    S net0 = deep_net<S>(pd, F.integ, F.eta, &F.eh.back()[5], h0);
    S g0   = deep_height_dt<S>(pd, F.integ, F.eta, &F.eh.back()[5], h0);
    S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay[i]);
    stand[i] = FullState<S>{plant::FF16State<S>{h0, -log(pr_estab), S(0), S(0), S(0)},
                            S(0), log(birth_rate * pr_estab / g0)};
    alive[i] = 1;
  }

  // Final-census reduction: the size-distribution integral = the height-trapezium of
  // density_i * psi(state_i) over the cohorts (descending height) + the pending-seed
  // tail. Mirrors assemble_metrics' census_reduce (the boundary seed via the FINAL
  // active env reconstructed from the final stand).
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){
    return as_double(stand[a].demog.height) > as_double(stand[b].demog.height); });

  // Pending-seed density via establishment in the FINAL env (frozen harvested tail).
  S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
  S net0 = deep_net<S>(pd, F.integ, F.eta, &F.eh.back()[5], h0);
  S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay.back());
  S g0 = deep_height_dt<S>(pd, F.integ, F.eta, &F.eh.back()[5], h0);
  S dens_new = birth_rate * pr_estab / g0;

  auto census_reduce = [&](auto psi) -> S {
    std::vector<S> phi(nC);
    for (std::size_t i = 0; i < nC; ++i)
      phi[i] = psi(stand[i].demog.height, exp(stand[i].log_density), stand[i].demog.mass_heartwood);
    S J = S(0.0);
    for (std::size_t j = 0; j + 1 < nC; ++j) {
      const std::size_t a = ord[j], b = ord[j+1];
      J = J + S(0.5) * (stand[a].demog.height - stand[b].demog.height) * (phi[a] + phi[b]);
    }
    if (nC > 0) {
      const std::size_t last = ord[nC-1];
      S phi_new = psi(h0, dens_new, S(0.0));
      J = J + S(0.5) * (stand[last].demog.height - h0) * (phi[last] + phi_new);
    }
    return J;
  };

  // net_reproduction_ratio (the demographic-equilibrium / mutant-fitness axis): the
  // per-seed expected offspring = sum_i w_i * offspring_i with UNIT birth weights
  // w_i = tw_i / birth_rate0 (the CONSTANT harvested birth_rate; tw was harvested
  // birth-scaled). birth_rate enters offspring_i ONLY through the active canopy, so the
  // reverse sweep w.r.t. birth_rate is the pure density feedback dR0/d(birth_rate) -- the
  // change in a (trait-identical) mutant's fitness as the resident density moves. No
  // pending-seed / new_node term (the pending seed has not reproduced).
  const double br0 = (F.birth_rate0 > 0.0) ? F.birth_rate0 : as_double(birth_rate);
  std::vector<S> J(metrics.size());
  for (std::size_t m = 0; m < metrics.size(); ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      J[m] = census_reduce([&](S h, S dens, S) -> S {
        return dens * S(F.kI) * plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h); });
    } else if (nm == "biomass") {
      J[m] = census_reduce([&](S h, S dens, S mhw) -> S {
        return dens * (plant::ff16_mass_live_given_height(pd, h) + mhw); });
    } else if (nm == "net_reproduction_ratio") {
      S R0 = S(0.0);
      for (std::size_t i = 0; i < nC; ++i)
        R0 = R0 + S(F.tw[i] / br0) * stand[i].offspring;
      J[m] = R0;
    } else { // size_moment
      J[m] = census_reduce([&](S h, S dens, S) -> S { return dens * h; });
    }
  }
  return J;
}

// Populate the COUPLED harvest on F: kI, the frozen knot x-positions + SCM knot light
// values per cached stage (from F.eh's light spline), and the boundary new_node
// (height + competition effect) per stage from the R-side harvest.
void attach_coupled(Frozen& F, double kI, double area,
                    Rcpp::List nn_h_list, Rcpp::List nn_c_list,
                    bool active_birthenv = true) {
  F.coupled = true; F.kI = kI; F.area = area;
  F.coupled_active_birthenv = active_birthenv;
  const std::size_t N = F.eh.size();
  F.knot_x.resize(N); F.knot_y0.resize(N);
  F.nn_h.resize(N); F.nn_c.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    const std::size_t ns = F.eh[n].size();
    F.knot_x[n].resize(ns); F.knot_y0[n].resize(ns);
    for (std::size_t s = 0; s < ns; ++s) {
      F.knot_x[n][s]  = F.eh[n][s].light_availability.spline.get_x();
      F.knot_y0[n][s] = F.eh[n][s].light_availability.spline.get_y();
    }
    F.nn_h[n] = Rcpp::as<std::vector<double>>(nn_h_list[n]);
    F.nn_c[n] = Rcpp::as<std::vector<double>>(nn_c_list[n]);
  }
}

// Native variant: the boundary new_node height/effect per RK stage come directly from
// the live Patch's stand_newnode_*_stage_history (native [step][stage] doubles), so the
// coupled path needs no R harvest at all. knot_x/knot_y0 are read from F.eh as above.
void attach_coupled_native(Frozen& F, double kI, double area,
                           const std::vector<std::vector<double>>& nn_h,
                           const std::vector<std::vector<double>>& nn_c,
                           bool active_birthenv = true) {
  F.coupled = true; F.kI = kI; F.area = area;
  F.coupled_active_birthenv = active_birthenv;
  const std::size_t N = F.eh.size();
  F.knot_x.resize(N); F.knot_y0.resize(N);
  F.nn_h.resize(N); F.nn_c.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    const std::size_t ns = F.eh[n].size();
    F.knot_x[n].resize(ns); F.knot_y0[n].resize(ns);
    for (std::size_t s = 0; s < ns; ++s) {
      F.knot_x[n][s]  = F.eh[n][s].light_availability.spline.get_x();
      F.knot_y0[n][s] = F.eh[n][s].light_availability.spline.get_y();
    }
    F.nn_h[n] = nn_h[n];
    F.nn_c[n] = nn_c[n];
  }
}

// ===========================================================================
// MULTI-SPECIES COUPLED replay (#472 scope B, R2 -- the cross-species resident
// Jacobian). All species' cohorts are re-evolved TOGETHER over the frozen schedule;
// the canopy light each RK stage is the JOINT reconstruction
//   competition(z) = (1/area) * sum_species [ trapezium over species s's active
//                    cohorts + its boundary node, with species s's eta ]
// (matching Patch::compute_competition = sum_s Species_s::compute_competition/area --
// per-species trapeziums summed, NOT one merged trapezium). Every census metric is
// the sum over species of that species' size-distribution trapezium. Differentiating
// w.r.t. ONE species' traits gives the cross-species total: the target species'
// traits move its cohorts, which re-shade the joint canopy that EVERY species reads,
// so the other species' contributions respond too (the cross term the frozen/mutant
// gradient sets to zero).
// ===========================================================================

// Per-species frozen harvest for the multi-species coupled replay. The schedule, the
// joint env spline knots (knot_x/knot_y0 from the joint eh) and the patch area are
// shared; everything else is per species.
struct FrozenMS {
  std::vector<std::vector<plant::FF16_Environment>> eh;          // joint env [step][0..5]
  std::vector<double> step_h;
  const plant::quadrature::QK* integ = nullptr;
  double area = 1.0;
  std::vector<std::vector<std::vector<double>>> knot_x, knot_y0; // [step][stage][knot]
  bool active_birthenv = true;

  std::size_t nS = 0;
  std::vector<plant::FF16ProdPars<double>> pd;                   // [species]
  std::vector<double> eta, kI, a_d0, h0, birth_rate;             // [species]
  std::vector<std::vector<int>> birth;                           // [species][cohort]
  std::vector<std::vector<double>> decay;                        // [species][cohort]
  std::vector<std::vector<std::vector<double>>> nn_h, nn_c;      // [species][step][stage]
};

// The coupled whole-stand replay + emergent-metric reduction, ALL SPECIES. Re-evolves
// every species' cohorts together, reconstructing the joint canopy each RK stage, and
// reduces the final census to the requested TOTAL-stand metrics (summed over species).
// pds carries the per-species prod_pars lifted to S; h0s the per-species seedling
// height (active for the differentiated species via the IFT injection). If env_err !=
// nullptr (double only) accumulates max|reconstructed knot light - SCM knot light|.
template <typename S>
std::vector<S> assemble_metrics_coupled_ms(const std::vector<plant::FF16ProdPars<S>>& pds,
    const FrozenMS& F, const std::vector<std::string>& metrics,
    const std::vector<S>& h0s, double birth_rate_unused = 0.0,
    double* env_err = nullptr, double* env_err_z = nullptr,
    const std::vector<S>* br_active = nullptr) {
  using std::exp; using std::log;
  using interp_t = odelia::interpolator::basic_interpolator<S>;
  const std::size_t nS = F.nS, N = F.eh.size();
  const double GEPS = 1e-6;                          // Control::node_gradient_eps
  (void)birth_rate_unused;
  // Per-species birth-rate driver: the harvested constant, or -- for the birth-rate
  // gradient -- an active scalar (one species registered as the tape input).
  auto BR = [&](std::size_t s) -> S {
    return br_active ? (*br_active)[s] : S(F.birth_rate[s]);
  };

  // Flatten cohorts across species into a global list (species, local, birth step).
  std::vector<std::size_t> gspec, glocal;
  std::vector<int> gbirth;
  for (std::size_t s = 0; s < nS; ++s)
    for (std::size_t i = 0; i < F.birth[s].size(); ++i) {
      gspec.push_back(s); glocal.push_back(i); gbirth.push_back(F.birth[s][i]);
    }
  const std::size_t nG = gspec.size();
  std::vector<CensusState<S>> stand(nG);
  std::vector<char> alive(nG, 0);

  auto env_idx = [&](std::size_t rn, int rs, std::size_t& en, int& es) {
    if (rs == 0) { if (rn > 0) { en = rn - 1; es = 5; } else { en = 0; es = 0; } }
    else { en = rn; es = rs - 1; }
  };

  // Joint canopy interpolator at (en,es)'s frozen knots: per-species trapezium summed.
  auto build_interp = [&](const std::vector<CensusState<S>>& y,
                          std::size_t en, int es) -> interp_t {
    const std::vector<double>& kx = F.knot_x[en][es];
    std::vector<std::vector<S>> hv(nS), gv(nS);
    for (std::size_t g = 0; g < nG; ++g) if (alive[g]) {
      const std::size_t s = gspec[g];
      S hi = y[g].demog.height;
      S dens = exp(y[g].log_density);
      S al = plant::ff16_area_leaf(pds[s].a_l1, pds[s].a_l2, hi);
      hv[s].push_back(hi); gv[s].push_back(dens * S(F.kI[s]) * al);
    }
    std::vector<S> comp(kx.size(), S(0.0));
    for (std::size_t s = 0; s < nS; ++s) {
      hv[s].push_back(S(F.nn_h[s][en][es]));            // boundary tail (per species)
      gv[s].push_back(S(F.nn_c[s][en][es]));
      std::vector<std::size_t> ord(hv[s].size());
      for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
      std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){
        return as_double(hv[s][a]) > as_double(hv[s][b]); });
      std::vector<S> hs(hv[s].size()), gs(hv[s].size());
      for (std::size_t k = 0; k < ord.size(); ++k) { hs[k] = hv[s][ord[k]]; gs[k] = gv[s][ord[k]]; }
      for (std::size_t k = 0; k < kx.size(); ++k)
        comp[k] = comp[k] + plant::gradient::canopy_comp_at<S>(kx[k], hs, gs, F.eta[s]);
    }
    std::vector<S> ly(kx.size());
    for (std::size_t k = 0; k < kx.size(); ++k) ly[k] = exp(-comp[k] * S(1.0 / F.area));
    if (env_err) {
      const std::vector<double>& y0 = F.knot_y0[en][es];
      for (std::size_t k = 0; k < ly.size() && k < y0.size(); ++k) {
        double e = std::abs(as_double(ly[k]) - y0[k]);
        if (e > *env_err) { *env_err = e; if (env_err_z) *env_err_z = kx[k]; }
      }
    }
    interp_t interp; interp.init(kx, ly);
    return interp;
  };

  auto deriv = [&](const std::vector<CensusState<S>>& y, std::size_t rn, int rs)
      -> std::vector<CensusState<S>> {
    std::size_t en; int es; env_idx(rn, rs, en, es);
    interp_t interp = build_interp(y, en, es);
    const double cap = F.knot_x[en][es].back();
    std::vector<CensusState<S>> dy(nG);
    for (std::size_t g = 0; g < nG; ++g) {
      if (!alive[g]) {
        dy[g] = CensusState<S>{plant::FF16State<S>{S(0),S(0),S(0),S(0),S(0)}, S(0)};
        continue;
      }
      const std::size_t s = gspec[g];
      const S h = y[g].demog.height;
      S area_leaf = plant::ff16_area_leaf(pds[s].a_l1, pds[s].a_l2, h);
      S net = deep_net_coupled<S>(pds[s], F.integ, F.eta[s], interp, cap, h);
      plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pds[s], h, area_leaf, net, true);
      S g_back = deep_height_dt_coupled<S>(pds[s], F.integ, F.eta[s], interp, cap, h - S(GEPS));
      S gprime = (r.height_dt - g_back) / S(GEPS);
      S log_density_dt = -gprime - r.mortality_dt;
      dy[g] = CensusState<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, log_density_dt};
    }
    return dy;
  };
  auto axpy = [&](const std::vector<CensusState<S>>& a, double c,
                  const std::vector<CensusState<S>>& k) -> std::vector<CensusState<S>> {
    std::vector<CensusState<S>> r(a.size());
    for (std::size_t g = 0; g < a.size(); ++g)
      r[g] = CensusState<S>{plant::FF16State<S>{
        a[g].demog.height + c * k[g].demog.height,
        a[g].demog.mortality + c * k[g].demog.mortality,
        a[g].demog.fecundity + c * k[g].demog.fecundity,
        a[g].demog.area_heartwood + c * k[g].demog.area_heartwood,
        a[g].demog.mass_heartwood + c * k[g].demog.mass_heartwood},
        a[g].log_density + c * k[g].log_density};
    return r;
  };

  for (std::size_t rn = 0; rn < N; ++rn) {
    bool any_born = false;
    for (std::size_t g = 0; g < nG; ++g)
      if ((std::size_t)gbirth[g] == rn) { any_born = true; break; }
    if (any_born) {
      std::size_t en; int es; env_idx(rn, 0, en, es);
      interp_t binterp;
      double bcap = 0.0;
      const plant::FF16_Environment* eb = (rn > 0) ? &F.eh[rn-1][5] : &F.eh[0][0];
      if (F.active_birthenv) { binterp = build_interp(stand, en, es); bcap = F.knot_x[en][es].back(); }
      for (std::size_t g = 0; g < nG; ++g) {
        if ((std::size_t)gbirth[g] != rn) continue;
        const std::size_t s = gspec[g];
        S h0 = h0s[s];
        S area_leaf_0 = plant::ff16_area_leaf(pds[s].a_l1, pds[s].a_l2, h0);
        S net0, g0;
        if (F.active_birthenv) {
          net0 = deep_net_coupled<S>(pds[s], F.integ, F.eta[s], binterp, bcap, h0);
          g0   = deep_height_dt_coupled<S>(pds[s], F.integ, F.eta[s], binterp, bcap, h0);
        } else {
          net0 = deep_net<S>(pds[s], F.integ, F.eta[s], eb, h0);
          g0   = deep_height_dt<S>(pds[s], F.integ, F.eta[s], eb, h0);
        }
        S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0[s],
                       F.decay[s][glocal[g]]);
        S mort0 = -log(pr_estab);
        S logd0 = log(BR(s) * pr_estab / g0);
        stand[g] = CensusState<S>{plant::FF16State<S>{h0, mort0, S(0), S(0), S(0)}, logd0};
        alive[g] = 1;
      }
    }
    stand = rkck_one_step(stand, F.step_h[rn], rn, deriv, axpy);
  }

  // Establish final-boundary cohorts (birth >= N) at h0 (see the single-species note in
  // assemble_metrics_coupled): never stepped, mirrors the frozen per-cohort replay;
  // prevents the a_l2 0*log(0)=NaN from a zero-initialised census height.
  for (std::size_t g = 0; g < nG; ++g) {
    if (alive[g]) continue;
    const std::size_t s = gspec[g];
    S h0 = h0s[s];
    S area_leaf_0 = plant::ff16_area_leaf(pds[s].a_l1, pds[s].a_l2, h0);
    S net0 = deep_net<S>(pds[s], F.integ, F.eta[s], &F.eh.back()[5], h0);
    S g0   = deep_height_dt<S>(pds[s], F.integ, F.eta[s], &F.eh.back()[5], h0);
    S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0[s],
                   F.decay[s][glocal[g]]);
    stand[g] = CensusState<S>{plant::FF16State<S>{h0, -log(pr_estab), S(0), S(0), S(0)},
                              log(BR(s) * pr_estab / g0)};
    alive[g] = 1;
  }

  // Total-stand reduction: sum over species of that species' size-distribution
  // trapezium (descending height within species + the species' pending-seed tail).
  std::vector<S> J(metrics.size(), S(0.0));
  for (std::size_t s = 0; s < nS; ++s) {
    std::vector<std::size_t> gs;                       // global indices of species s
    for (std::size_t g = 0; g < nG; ++g) if (gspec[g] == s) gs.push_back(g);
    const std::size_t nc = gs.size();
    std::vector<std::size_t> ord(nc);
    for (std::size_t i = 0; i < nc; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){
      return as_double(stand[gs[a]].demog.height) > as_double(stand[gs[b]].demog.height); });
    // pending seed for species s, established in the frozen final joint env.
    S h0 = h0s[s];
    S area_leaf_0 = plant::ff16_area_leaf(pds[s].a_l1, pds[s].a_l2, h0);
    S net0 = deep_net<S>(pds[s], F.integ, F.eta[s], &F.eh.back()[5], h0);
    S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0[s],
                   F.decay[s].back());
    S g0 = deep_height_dt<S>(pds[s], F.integ, F.eta[s], &F.eh.back()[5], h0);
    S dens_new = BR(s) * pr_estab / g0;

    auto census_reduce = [&](auto psi) -> S {
      std::vector<S> phi(nc);
      for (std::size_t i = 0; i < nc; ++i) {
        const CensusState<S>& y = stand[gs[i]];
        phi[i] = psi(y.demog.height, exp(y.log_density), y.demog.mass_heartwood);
      }
      S Js = S(0.0);
      for (std::size_t j = 0; j + 1 < nc; ++j) {
        const std::size_t a = ord[j], b = ord[j+1];
        Js = Js + S(0.5) * (stand[gs[a]].demog.height - stand[gs[b]].demog.height) *
                  (phi[a] + phi[b]);
      }
      if (nc > 0) {
        const std::size_t last = ord[nc-1];
        S phi_new = psi(h0, dens_new, S(0.0));
        Js = Js + S(0.5) * (stand[gs[last]].demog.height - h0) * (phi[last] + phi_new);
      }
      return Js;
    };

    for (std::size_t m = 0; m < metrics.size(); ++m) {
      const std::string& nm = metrics[m];
      if (nm == "LAI") {
        J[m] = J[m] + census_reduce([&](S h, S dens, S) -> S {
          return dens * S(F.kI[s]) * plant::ff16_area_leaf(pds[s].a_l1, pds[s].a_l2, h); });
      } else if (nm == "biomass") {
        J[m] = J[m] + census_reduce([&](S h, S dens, S mhw) -> S {
          return dens * (plant::ff16_mass_live_given_height(pds[s], h) + mhw); });
      } else { // size_moment
        J[m] = J[m] + census_reduce([&](S h, S dens, S) -> S { return dens * h; });
      }
    }
  }
  return J;
}

// Build the multi-species frozen harvest from per-species R-side arrays. pp_list:
// per-species parameter vectors; birth_list: per-species cohort birth steps; nn_h/c:
// the all-species boundary harvest [step][stage][species]; sh: shared step times.
FrozenMS build_frozen_ms(Rcpp::List pp_list, Rcpp::List eh_list,
    const std::vector<double>& sh, Rcpp::List birth_list,
    std::vector<double> birth_rate, Rcpp::List nn_h_list, Rcpp::List nn_c_list,
    double patch_area, bool active_birthenv) {
  FrozenMS F;
  F.area = patch_area; F.active_birthenv = active_birthenv;
  const std::size_t nS = pp_list.size(); F.nS = nS;
  const std::size_t N = eh_list.size();
  F.eh.resize(N); F.step_h.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      F.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));
  }
  for (std::size_t n = 0; n < N; ++n) F.step_h[n] = sh[n+1] - sh[n];
  // Joint env spline knots (shared canopy).
  F.knot_x.resize(N); F.knot_y0.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    const std::size_t ns = F.eh[n].size();
    F.knot_x[n].resize(ns); F.knot_y0[n].resize(ns);
    for (std::size_t s = 0; s < ns; ++s) {
      F.knot_x[n][s]  = F.eh[n][s].light_availability.spline.get_x();
      F.knot_y0[n][s] = F.eh[n][s].light_availability.spline.get_y();
    }
  }
  // Per-species scalars, prod_pars, birth steps, decay, boundary.
  F.pd.resize(nS); F.eta.resize(nS); F.kI.resize(nS); F.a_d0.resize(nS);
  F.h0.resize(nS); F.birth_rate.resize(nS); F.birth.resize(nS); F.decay.resize(nS);
  F.nn_h.resize(nS); F.nn_c.resize(nS);
  // F.integ is the Gauss-Kronrod quadrature rule -- strategy-independent (the same
  // nodes/weights regardless of biology), so the caller sets it from its own strategy
  // (whose lifetime spans the assemble); we only copy out per-species values here.
  F.integ = nullptr;
  for (std::size_t s = 0; s < nS; ++s) {
    Rcpp::NumericVector pp = pp_list[s];
    plant::FF16_Strategy st = make_strategy(pp);
    F.pd[s] = st.prod_pars();
    F.eta[s] = st.pars.eta; F.kI[s] = st.pars.k_I; F.a_d0[s] = st.pars.a_d0;
    F.h0[s] = st.initial_height(); F.birth_rate[s] = birth_rate[s];
    F.birth[s] = Rcpp::as<std::vector<int>>(birth_list[s]);
    F.decay[s].resize(F.birth[s].size());
    for (std::size_t i = 0; i < F.birth[s].size(); ++i)
      F.decay[s][i] = std::exp(-st.pars.recruitment_decay * sh[(std::size_t)F.birth[s][i]]);
    // Reshape the all-species boundary harvest [step][stage][species] -> [step][stage].
    F.nn_h[s].resize(N); F.nn_c[s].resize(N);
    for (std::size_t n = 0; n < N; ++n) {
      Rcpp::List hns = nn_h_list[n], cns = nn_c_list[n];
      const std::size_t ns = hns.size();
      F.nn_h[s][n].resize(ns); F.nn_c[s][n].resize(ns);
      for (std::size_t k = 0; k < ns; ++k) {
        std::vector<double> hv = Rcpp::as<std::vector<double>>(hns[k]);
        std::vector<double> cv = Rcpp::as<std::vector<double>>(cns[k]);
        F.nn_h[s][n][k] = (s < hv.size()) ? hv[s] : 0.0;
        F.nn_c[s][n][k] = (s < cv.size()) ? cv[s] : 0.0;
      }
    }
  }
  return F;
}

// Native-SCM multi-species frozen harvest: the joint env, step schedule, per-species
// birth steps and the all-species boundary harvest are read DIRECTLY from the live
// Patch (no R ff16_harvest_ms, no Rcpp::as<> env round-trip) -- the multi-species
// counterpart of build_frozen_scm. pp_list (cheap, from $parameters$strategies) gives
// each species' prod_pars / scalars; birth_rate is the per-species constant driver
// (recovered cheaply in R from offspring_production / net_reproduction_ratios).
// Bit-identical to build_frozen_ms (same arithmetic, faithful native data).
FrozenMS build_frozen_ms_scm(
    Rcpp::List pp_list,
    const plant::Patch<plant::FF16_Strategy, plant::FF16_Environment>& patch,
    std::vector<double> birth_rate, double patch_area, bool active_birthenv) {
  FrozenMS F;
  F.area = patch_area; F.active_birthenv = active_birthenv;
  const std::size_t nS = pp_list.size(); F.nS = nS;
  const auto& EH = patch.environment_history;
  const auto& sh = patch.step_history;
  const std::size_t N = EH.size();
  F.eh = EH;                                          // faithful copy (no Rcpp::as<>)
  F.step_h.resize(N);
  for (std::size_t n = 0; n < N; ++n) F.step_h[n] = sh[n + 1] - sh[n];
  // Joint env spline knots (shared canopy), straight off the native env splines.
  F.knot_x.resize(N); F.knot_y0.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    const std::size_t ns = F.eh[n].size();
    F.knot_x[n].resize(ns); F.knot_y0[n].resize(ns);
    for (std::size_t s = 0; s < ns; ++s) {
      F.knot_x[n][s]  = F.eh[n][s].light_availability.spline.get_x();
      F.knot_y0[n][s] = F.eh[n][s].light_availability.spline.get_y();
    }
  }
  F.pd.resize(nS); F.eta.resize(nS); F.kI.resize(nS); F.a_d0.resize(nS);
  F.h0.resize(nS); F.birth_rate.resize(nS); F.birth.resize(nS); F.decay.resize(nS);
  F.nn_h.resize(nS); F.nn_c.resize(nS);
  F.integ = nullptr;
  const auto& NNH = patch.stand_newnode_height_stage_history_all;       // [step][stage][species]
  const auto& NNC = patch.stand_newnode_competition_stage_history_all;
  for (std::size_t s = 0; s < nS; ++s) {
    Rcpp::NumericVector pp = pp_list[s];
    plant::FF16_Strategy st = make_strategy(pp);
    F.pd[s] = st.prod_pars();
    F.eta[s] = st.pars.eta; F.kI[s] = st.pars.k_I; F.a_d0[s] = st.pars.a_d0;
    F.h0[s] = st.initial_height(); F.birth_rate[s] = birth_rate[s];
    F.birth[s] = plant::gradient::birth_steps(patch, s);                // native birth steps
    F.decay[s].resize(F.birth[s].size());
    for (std::size_t i = 0; i < F.birth[s].size(); ++i)
      F.decay[s][i] = std::exp(-st.pars.recruitment_decay * sh[(std::size_t)F.birth[s][i]]);
    // All-species boundary harvest [step][stage][species] -> [species][step][stage].
    F.nn_h[s].resize(N); F.nn_c[s].resize(N);
    for (std::size_t n = 0; n < N; ++n) {
      const std::size_t ns = NNH[n].size();
      F.nn_h[s][n].resize(ns); F.nn_c[s][n].resize(ns);
      for (std::size_t k = 0; k < ns; ++k) {
        F.nn_h[s][n][k] = (s < NNH[n][k].size()) ? NNH[n][k][s] : 0.0;
        F.nn_c[s][n][k] = (s < NNC[n][k].size()) ? NNC[n][k][s] : 0.0;
      }
    }
  }
  return F;
}

// ===========================================================================
// grow_individual_to_size trait gradient (#472 scope B, the last FF16 surface).
// A single plant grown in a FIXED environment to a target size, differentiated
// w.r.t. traits. No resident feedback (the env is given), so the only machinery
// beyond the SCM replay is the stopping-time implicit-function step.
//
// Pass 1 (double, R): grow_individual_to_size finds t* (height==target) and hands
//   us the adaptive step schedule + initial state via grow_individual_bracket.
// Pass 2 (here): replay the demographic ODE over the FROZEN schedule reading the
//   FIXED env (deep-crown) to a single partial final step landing on t*. The
//   schedule interval containing t* is integrated as exactly ONE RKCK step whether
//   it is the partial step or (post node-crossing) a full step, so the function is
//   C1-smooth across node boundaries and its FD reference matches AD to the floor.
//   One reverse sweep per state component gives the PARTIAL d(state at t*)/d(theta)
//   holding t* fixed; the stopping time responds to the trait via the IFT on
//   size(t*, theta) = target:  d(t*)/d(theta) = -(d size/d theta|t*) / size_dt(t*),
//   and the TOTAL derivative of each returned component is
//      d y_c/d theta = d y_c/d theta|t* + y_dot_c(t*) * d(t*)/d(theta).
//   For the size component itself the two terms cancel (it is pinned to target).
// ===========================================================================

// 5-state FF16 demographic derivative at a state, reading the FIXED env (deep-crown,
// the FF16 default). mortality_finite frozen true (matches the live grow path).
template <typename S>
plant::FF16State<S> grow_deriv_fixed(const plant::FF16ProdPars<S>& pd,
    const plant::quadrature::QK* integ, double eta, const plant::FF16_Environment* e,
    const plant::FF16State<S>& st) {
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, st.height);
  S net = deep_net<S>(pd, integ, eta, e, st.height);
  plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, st.height, area_leaf, net, true);
  return plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
                             r.area_heartwood_dt, r.mass_heartwood_dt};
}
template <typename S>
plant::FF16State<S> grow_axpy(const plant::FF16State<S>& a, double c, const plant::FF16State<S>& k) {
  return plant::FF16State<S>{a.height+c*k.height, a.mortality+c*k.mortality,
    a.fecundity+c*k.fecundity, a.area_heartwood+c*k.area_heartwood, a.mass_heartwood+c*k.mass_heartwood};
}
// Replay over a frozen step schedule reading the FIXED env (RKCK), from step 0.
template <typename S>
plant::FF16State<S> grow_replay_fixed(const plant::FF16ProdPars<S>& pd,
    const plant::quadrature::QK* integ, double eta, const plant::FF16_Environment* e,
    plant::FF16State<S> y, const std::vector<double>& sched) {
  auto deriv=[&](const plant::FF16State<S>& st, std::size_t, int){ return grow_deriv_fixed<S>(pd,integ,eta,e,st); };
  auto axpy=[](const plant::FF16State<S>& a, double c, const plant::FF16State<S>& k){ return grow_axpy<S>(a,c,k); };
  return plant::ff16_cashkarp_replay(y, sched, 0, deriv, axpy);
}
// One FF16State component by index (0..4): height, mortality, fecundity, area_hw, mass_hw.
double grow_state_at(const plant::FF16State<double>& s, int c) {
  switch (c) { case 0: return s.height; case 1: return s.mortality; case 2: return s.fecundity;
               case 3: return s.area_heartwood; default: return s.mass_heartwood; }
}
// Double discovery: integrate y0 over step_h (fixed env) until the size component
// `sidx` crosses `target`; return the full-step count before crossing (nfull) and the
// bisected partial final step dt_final so the trajectory lands on size==target.
void grow_discover(const plant::FF16ProdPars<double>& pd, const plant::quadrature::QK* integ,
    double eta, const plant::FF16_Environment* e, const plant::FF16State<double>& y0,
    const std::vector<double>& step_h, double target, int sidx,
    int& nfull, double& dt_final, plant::FF16State<double>& final_state) {
  plant::FF16State<double> y = y0;
  for (std::size_t n = 0; n < step_h.size(); ++n) {
    plant::FF16State<double> ynext = grow_replay_fixed<double>(pd, integ, eta, e, y, {step_h[n]});
    if (grow_state_at(ynext, sidx) >= target) {
      nfull = (int)n;
      double lo = 0.0, hi = step_h[n];
      for (int it = 0; it < 80; ++it) {
        double mid = 0.5 * (lo + hi);
        double sm = grow_state_at(grow_replay_fixed<double>(pd, integ, eta, e, y, {mid}), sidx);
        if (sm < target) lo = mid; else hi = mid;
      }
      dt_final = 0.5 * (lo + hi);
      final_state = grow_replay_fixed<double>(pd, integ, eta, e, y, {dt_final});
      return;
    }
    y = ynext;
  }
  nfull = -1;  // target not reached within the schedule
}

} // namespace

// Compiled core of offspring_production_gradient(). Takes the harvested resident
// schedule (env per RK stage, step sizes), the per-cohort birth steps / weights /
// survival, and the trait names to differentiate. Returns d(offspring_production)/
// d(trait), establishment frozen. eh_list is steps x 6 FF16_Environment objects.
// [[Rcpp::export]]
Rcpp::NumericVector ff16_offspring_production_gradient_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits) {
  auto s = make_strategy(pp);
  auto pd = s.prod_pars();
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);

  // Per-cohort reverse sweep (cohorts independent under the frozen resident).
  std::vector<double> g; double J;
  accumulate_offspring_gradient(pd, F, h0v, dh0, idx, g, J);

  Rcpp::NumericVector grad(idx.size());
  for (std::size_t k=0;k<idx.size();++k) grad[k] = g[k];
  grad.attr("names") = Rcpp::wrap(traits);
  grad.attr("offspring_production") = J;
  return grad;
}

// Compiled core of the generic stand-gradient engine (#472 scope B, build-order
// step 1). Records ONE forward replay of every cohort (final demographic state +
// offspring accumulator + log_density) onto a single adjoint tape, then for EACH
// requested metric takes one reverse sweep, giving the metrics x traits Jacobian
// (+ the metric values) out of one resident baseline. No metric is privileged.
// Two symmetric reduction kinds, both weighted reductions over the replayed cohorts:
//   - Lagrangian (offspring_production): J = sum_i w_i * f(state_i), frozen w_i;
//   - census (LAI / biomass / size_moment): the size-distribution integral at the
//     final census = the height-trapezium of density_i * psi(state_i) over the
//     cohorts (descending height) + the pending-seed tail term, matching the SCM's
//     compute_competition arithmetic. density_i = exp(log_density_i) is the replayed
//     census number density (log_density evolved with the SCM's own g' scheme).
// Shared core of the FF16 stand-gradient entry points: given an already-built
// Frozen harvest, compute the metrics x traits Jacobian + values. Split out so the
// R-list (`_impl`) and native-SCM (`_native`) entries differ ONLY in how F's env is
// sourced (lossy Rcpp::as<> vs faithful native copy), not in the replay/sweep.
Rcpp::List ff16_stand_gradient_core(
    Frozen& F, const plant::FF16_Strategy& s, std::vector<std::string> traits,
    std::vector<std::string> metrics, double birth_rate, std::string feedback,
    Rcpp::List sh_h_list, Rcpp::List sh_c_list, double patch_area,
    double al1_base, double al2_base) {
  auto pd = s.prod_pars();
  // Weight basis: the per-node weight C_i = ce_i / area_leaf(al1_base, al2_base, h_i)
  // is frozen at the RESIDENT's base allometry, independent of the differentiated
  // a_l1/a_l2 in `pp` -- so that a FD that perturbs pp's allometry (the R1
  // reconstruction-FD) actually moves the recon (C_i fixed, area_leaf active) instead
  // of cancelling. Defaults (al1_base<0) to pp's own a_l1/a_l2 (the normal call).
  const double alb1 = (al1_base > 0) ? al1_base : pd.a_l1;
  const double alb2 = (al2_base > 0) ? al2_base : pd.a_l2;
  if (feedback == "resident") {
    attach_resident_stand(F, sh_h_list, sh_c_list, alb1, alb2, patch_area, true);
  } else if (feedback == "resident_noanchor") {
    attach_resident_stand(F, sh_h_list, sh_c_list, alb1, alb2, patch_area, false);
  } else if (feedback != "frozen") {
    Rcpp::stop("unknown feedback: " + feedback +
               " (expected 'frozen', 'resident' or 'resident_noanchor')");
  }
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);
  const double kI = s.pars.k_I;
  const std::size_t M = metrics.size(), nT = idx.size();
  bool need_census = false;
  for (auto& nm : metrics) {
    if (nm != "offspring_production" && nm != "LAI" && nm != "biomass" &&
        nm != "size_moment")
      Rcpp::stop("unknown stand metric: " + nm);
    if (nm != "offspring_production") need_census = true;
  }

  // Fast path: with no census metric requested every metric is offspring_production,
  // a per-cohort-independent reduction under the frozen resident -> tape one cohort at
  // a time (tiny tapes) rather than the whole stand at once. Bit-identical, ~25x faster.
  if (!need_census) {
    std::vector<double> g; double val;
    accumulate_offspring_gradient(pd, F, h0v, dh0, idx, g, val);
    Rcpp::NumericMatrix jac(M, nT);
    Rcpp::NumericVector values(M);
    for (std::size_t m = 0; m < M; ++m) {
      for (std::size_t k = 0; k < nT; ++k) jac(m, k) = g[k];
      values[m] = val;
    }
    jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
    values.attr("names") = Rcpp::wrap(metrics);
    return Rcpp::List::create(Rcpp::Named("jacobian") = jac,
                              Rcpp::Named("values")   = values);
  }

  // --- FROZEN baseline (reverse mode): all traits, one tape, one sweep per metric.
  // The resident feedback enters only through area_leaf (a_l1/a_l2), so the reverse
  // baseline runs the FROZEN env (no O(stand) recon on the tape -> no blow-up); the
  // a_l1/a_l2 columns get the resident TOTAL grafted on by forward mode below.
  const bool want_resident = F.resident;
  F.resident = false;
  Rcpp::NumericMatrix jac(M, nT);
  Rcpp::NumericVector values(M);
  {
    // Scope the reverse tape so it is DESTROYED before the forward pass runs below
    // (keep the adj global-tape state strictly separate from the fwd pass).
    ad::tape_type tape;
    auto pa = lift<ad_t>(pd);
    auto fp = field_ptrs<ad_t>(pa);
    for (auto i : idx) tape.registerInput(*fp[i]);
    tape.newRecording();
    ad_t h0 = inject_h0(h0v, fp, idx, dh0);   // h0 active via the IFT first-order injection
    std::vector<ad_t> J = assemble_metrics<ad_t>(pa, F, metrics, birth_rate, kI, h0);
    for (std::size_t m = 0; m < M; ++m) tape.registerOutput(J[m]);

    // One cheap reverse sweep PER metric over the single recording (clear adjoints
    // between, the XAD multi-output Jacobian pattern).
    for (std::size_t m = 0; m < M; ++m) {
      tape.clearDerivatives();
      xad::derivative(J[m]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < nT; ++k) jac(m, k) = xad::derivative(*fp[idx[k]]);
      values[m] = as_double(J[m]);
    }
  }
  F.resident = want_resident;

  // --- RESIDENT feedback (forward mode): graft the resident TOTAL d(metric)/d(a_l1)
  // and /d(a_l2) over the FROZEN baseline. At frozen geometry the canopy-light formula
  // Sum density*k_I*area_leaf(a_l1,a_l2,h)*Q(eta) contains, of the 28 DIFFERENTIABLE
  // traits, ONLY a_l1/a_l2 (k_I and raw eta are strategy constants, not in the trait
  // vector) -- so this graft is COMPLETE for the differentiable set; the other 26 do
  // not enter the light at fixed geometry, hence resident == frozen for them exactly.
  // One forward pass per direction (no tape) gives the total for all metrics at once.
  if (want_resident) {
    const auto names = field_names();
    for (std::size_t k = 0; k < nT; ++k) {
      const std::string& tn = names[idx[k]];
      if (tn != "a_l1" && tn != "a_l2") continue;
      auto pf  = lift<fwd_t>(pd);
      auto fpf = field_ptrs<fwd_t>(pf);
      xad::derivative(*fpf[idx[k]]) = 1.0;                 // seed this trait direction
      fwd_t h0f = fwd_t(h0v); xad::derivative(h0f) = dh0[k];  // IFT seedling-size channel
      std::vector<fwd_t> Jf = assemble_metrics<fwd_t>(pf, F, metrics, birth_rate, kI, h0f);
      for (std::size_t m = 0; m < M; ++m) jac(m, k) = xad::derivative(Jf[m]);
    }
    // For feedback="resident_noanchor", report the genuine recon VALUES (theta-
    // dependent, ~1e-4 off the SCM) so a finite-difference can validate the forward
    // gradient; feedback="resident" keeps the exact (frozen) values.
    if (!F.anchor) {
      std::vector<double> Jd = assemble_metrics<double>(pd, F, metrics, birth_rate, kI, h0v);
      for (std::size_t m = 0; m < M; ++m) values[m] = Jd[m];
    }
  }

  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac,
                            Rcpp::Named("values")   = values);
}

// [[Rcpp::export]]
Rcpp::List ff16_stand_gradient_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits,
    std::vector<std::string> metrics, double birth_rate,
    std::string feedback, Rcpp::List sh_h_list, Rcpp::List sh_c_list,
    double patch_area, double al1_base, double al2_base) {
  auto s = make_strategy(pp);
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  return ff16_stand_gradient_core(F, s, traits, metrics, birth_rate, feedback,
                                  sh_h_list, sh_c_list, patch_area, al1_base, al2_base);
}

// FULLY native-SCM variant: builds the WHOLE harvest (env + schedule + birth steps +
// weights + per-stage survival) from the live Patch -- no R-side ff16_harvest at all,
// so the O(stand) pr_survival loop (R/emergent_gradient.R:70) runs once in C++ on the
// native patch. `species` is 0-based; birth_rate < 0 recovers the constant rate from
// the run (offspring_production / net_reproduction_ratio). sh_h_list/sh_c_list are the
// (rare) resident_noanchor stand-stage harvest; empty for the frozen path. `scm_` is
// the RcppR6 SCM<FF16,FF16_Env>; Rcpp::as only unwraps its .ptr (no serialization).
// [[Rcpp::export]]
Rcpp::List ff16_stand_gradient_native(
    SEXP scm_, Rcpp::NumericVector pp, int species,
    std::vector<std::string> traits, std::vector<std::string> metrics,
    double birth_rate, std::string feedback, Rcpp::List sh_h_list,
    Rcpp::List sh_c_list, double patch_area, double al1_base, double al2_base) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::FF16_Strategy, plant::FF16_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  auto s = make_strategy(pp);
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  Frozen F = build_frozen_scm(s, patch, (std::size_t)species, br);
  return ff16_stand_gradient_core(F, s, traits, metrics, br, feedback,
                                  sh_h_list, sh_c_list, patch_area, al1_base, al2_base);
}

// R0 gate for the COUPLED resident replay (#472 scope B): a double-precision whole-
// stand re-evolution that reconstructs the active canopy per RK stage and reduces the
// emergent metrics. Returns the reconstructed metric VALUES and `env_err` = the worst
// |reconstructed knot light - SCM knot light| over all stages (the coupled-drift
// gauge -- should be ~3e-14 if the re-evolution reproduces the SCM stand). nn_h_list/
// nn_c_list are the boundary new_node height/effect per cached stage (R harvest).
// [[Rcpp::export]]
Rcpp::List ff16_coupled_metrics_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> metrics, double birth_rate,
    Rcpp::List nn_h_list, Rcpp::List nn_c_list, double patch_area,
    bool active_birthenv = true, double birth_rate0 = -1.0) {
  auto s  = make_strategy(pp);
  auto pd = s.prod_pars();
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  attach_coupled(F, s.pars.k_I, patch_area, nn_h_list, nn_c_list, active_birthenv);
  // net_reproduction_ratio unit weights are tw/birth_rate0: pass birth_rate0 = the FIXED
  // harvest br so an FD that perturbs `birth_rate` (the canopy/density driver) holds the
  // per-seed weight fixed -- matching the AD. Defaults to birth_rate (census FD callers).
  F.birth_rate0 = (birth_rate0 > 0.0) ? birth_rate0 : birth_rate;
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment" &&
        nm != "net_reproduction_ratio")
      Rcpp::stop("coupled metrics: expected LAI / biomass / size_moment / "
                 "net_reproduction_ratio, got " + nm);
  const double h0v = s.initial_height();
  double env_err = 0.0, env_err_z = -1.0;
  std::vector<double> J = assemble_metrics_coupled<double>(pd, F, metrics, birth_rate,
                                                           h0v, &env_err, &env_err_z);
  Rcpp::NumericVector values(J.size());
  for (std::size_t m = 0; m < J.size(); ++m) values[m] = J[m];
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("values") = values,
                            Rcpp::Named("env_err") = env_err,
                            Rcpp::Named("env_err_z") = env_err_z);
}

// R1 of the COUPLED resident replay (#472 scope B): the resident TOTAL trait gradient
// of the emergent census metrics, one reverse sweep per metric over the coupled
// whole-stand re-evolution. Unlike the frozen-geometry graft (a_l1/a_l2 only), EVERY
// trait feeds back here -- a trait that changes growth/mortality moves the cohorts'
// heights/densities, which re-shade the active canopy every cohort reads. One forward
// replay onto a single adjoint tape, then a cheap reverse sweep per metric gives the
// full metrics x traits Jacobian. Validate vs an FD over ff16_coupled_metrics_impl
// (the SAME coupled reconstruction, frozen geometry) -- AD and FD differentiate one
// function. nn_h_list/nn_c_list: boundary new_node height/effect per cached stage.
// Shared core of the coupled-gradient entries: F already has the coupled canopy
// attached. The R-list (`_impl`) and native-SCM (`_native`) entries differ only in
// how F is built (lossy eh_list + R harvest vs faithful native Patch read).
Rcpp::List ff16_coupled_gradient_core(
    Frozen& F, const plant::FF16_Strategy& s, std::vector<std::string> traits,
    std::vector<std::string> metrics, double birth_rate) {
  auto pd = s.prod_pars();
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("coupled gradient: expected LAI / biomass / size_moment, got " + nm);
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);
  const std::size_t M = metrics.size(), nT = idx.size();

  Rcpp::NumericMatrix jac(M, nT);
  Rcpp::NumericVector values(M);
  {
    // One coupled tape, scoped so it is destroyed before returning. Reverse mode:
    // many traits in, few metrics out -> one recording, one sweep per metric.
    ad::tape_type tape;
    auto pa = lift<ad_t>(pd);
    auto fp = field_ptrs<ad_t>(pa);
    for (auto i : idx) tape.registerInput(*fp[i]);
    tape.newRecording();
    ad_t h0 = inject_h0(h0v, fp, idx, dh0);            // IFT seedling-size channel
    std::vector<ad_t> J = assemble_metrics_coupled<ad_t>(pa, F, metrics, birth_rate, h0);
    for (std::size_t m = 0; m < M; ++m) tape.registerOutput(J[m]);
    for (std::size_t m = 0; m < M; ++m) {
      tape.clearDerivatives();
      xad::derivative(J[m]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < nT; ++k) jac(m, k) = xad::derivative(*fp[idx[k]]);
      values[m] = as_double(J[m]);
    }
  }
  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac,
                            Rcpp::Named("values")   = values);
}

// FULLY native coupled gradient: the whole harvest + boundary new_node history come
// from the live Patch (no R ff16_harvest, no Rcpp::as<> env). `species` 0-based;
// birth_rate<0 recovers natively. This is the resident TOTAL trait gradient (every
// trait re-shades the canopy) -- the selection-gradient path.
// [[Rcpp::export]]
Rcpp::List ff16_coupled_gradient_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, std::vector<std::string> traits,
    std::vector<std::string> metrics, double birth_rate, double patch_area,
    bool active_birthenv = true) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::FF16_Strategy, plant::FF16_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  if (patch.stand_newnode_height_stage_history.size() < 1)
    Rcpp::stop("feedback = 'resident' needs the per-RK-stage boundary-node harvest; "
               "re-run the resident SCM with control(save_RK45_cache = TRUE)");
  auto s = make_strategy(pp);
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  Frozen F = build_frozen_scm(s, patch, (std::size_t)species, br);
  attach_coupled_native(F, s.pars.k_I, patch_area,
                        patch.stand_newnode_height_stage_history,
                        patch.stand_newnode_competition_stage_history, active_birthenv);
  return ff16_coupled_gradient_core(F, s, traits, metrics, br);
}

// Resident-coupled birth_rate gradient core (#472 scope B, the birth_rate AXIS, not a
// trait). d(census metric)/d(birth_rate) with the FULL resident feedback: birth_rate
// scales every cohort's establishment density (logd0 = log(birth_rate*pr_estab/g0)), so
// it re-shades the active canopy the whole stand reads -- the same coupled re-evolution
// the trait gradient uses, but with birth_rate as the SOLE registered tape input (the
// traits stay constant). The use case is evolving toward demographic equilibrium: the
// derivative gives the Newton step that drives an emergent census target. (The FROZEN /
// invasion reading is trivially analytic -- every density is linear in birth_rate, so
// d(metric)/d(birth_rate) = metric/birth_rate -- which is why only the resident-feedback
// axis needs a tape.) One reverse sweep per metric.
Rcpp::List ff16_birth_rate_gradient_core(
    Frozen& F, const plant::FF16_Strategy& s, std::vector<std::string> metrics,
    double birth_rate) {
  auto pd = s.prod_pars();
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment" &&
        nm != "net_reproduction_ratio")
      Rcpp::stop("birth_rate gradient: expected LAI / biomass / size_moment / "
                 "net_reproduction_ratio, got " + nm);
  const double h0v = s.initial_height();
  const std::size_t M = metrics.size();
  Rcpp::NumericVector dbr(M), values(M);
  {
    // One coupled tape; birth_rate is the only input. Traits are lifted constant.
    ad::tape_type tape;
    ad_t br = birth_rate;
    tape.registerInput(br);
    tape.newRecording();
    auto pa = lift<ad_t>(pd);
    std::vector<ad_t> J = assemble_metrics_coupled<ad_t>(pa, F, metrics, br, ad_t(h0v));
    for (std::size_t m = 0; m < M; ++m) tape.registerOutput(J[m]);
    for (std::size_t m = 0; m < M; ++m) {
      tape.clearDerivatives();
      xad::derivative(J[m]) = 1.0;
      tape.computeAdjoints();
      dbr[m] = xad::derivative(br);
      values[m] = as_double(J[m]);
    }
  }
  dbr.attr("names") = Rcpp::wrap(metrics);
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("d_birth_rate") = dbr,
                            Rcpp::Named("values")       = values);
}

// FULLY native resident-coupled birth_rate gradient: harvest + boundary new_node history
// from the live Patch (no R ff16_harvest, no Rcpp::as<> env). `species` 0-based;
// birth_rate<0 recovers natively. Returns d(census metric)/d(birth_rate) (the resident
// TOTAL, every cohort's density re-shading the canopy) + the reconstructed values.
// [[Rcpp::export]]
Rcpp::List ff16_birth_rate_gradient_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, std::vector<std::string> metrics,
    double birth_rate, double patch_area, bool active_birthenv = true) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::FF16_Strategy, plant::FF16_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  if (patch.stand_newnode_height_stage_history.size() < 1)
    Rcpp::stop("the resident birth_rate gradient needs the per-RK-stage boundary-node "
               "harvest; re-run the resident SCM with control(save_RK45_cache = TRUE)");
  auto s = make_strategy(pp);
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  Frozen F = build_frozen_scm(s, patch, (std::size_t)species, br);
  attach_coupled_native(F, s.pars.k_I, patch_area,
                        patch.stand_newnode_height_stage_history,
                        patch.stand_newnode_competition_stage_history, active_birthenv);
  return ff16_birth_rate_gradient_core(F, s, metrics, br);
}

// Cross-species resident birth_rate gradient core: d(TOTAL-stand census metric)/
// d(birth_rate of species `target`). The joint canopy is re-evolved from EVERY species'
// cohorts; only the target species' birth_rate is registered as the tape input, so the
// reverse sweep captures the genuinely new CROSS term -- the target's recruitment density
// re-shades the canopy every species reads (the frozen reading is the diagonal identity
// metric_s/birth_rate_s with zero cross term). One reverse sweep per metric; returns the
// gradient + reconstructed values + env_err (the caller gates the stiff-schedule case).
Rcpp::List ff16_birth_rate_gradient_ms_core(FrozenMS& F, std::vector<std::string> metrics,
    int target) {
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("birth_rate MS gradient: expected LAI / biomass / size_moment, got " + nm);
  const std::size_t tgt = (std::size_t)(target - 1);
  if (tgt >= F.nS) Rcpp::stop("target species out of range");
  const std::size_t M = metrics.size(), nS = F.nS;

  // R0 double pass for env_err (the joint re-evolution must reproduce the resident stand).
  double env_err = 0.0, env_err_z = -1.0;
  (void)assemble_metrics_coupled_ms<double>(F.pd, F, metrics, F.h0, 0.0,
                                            &env_err, &env_err_z);

  Rcpp::NumericVector dbr(M), values(M);
  {
    ad::tape_type tape;
    ad_t br_t = F.birth_rate[tgt];
    tape.registerInput(br_t);
    tape.newRecording();
    std::vector<plant::FF16ProdPars<ad_t>> pas(nS);
    for (std::size_t s = 0; s < nS; ++s) pas[s] = lift<ad_t>(F.pd[s]);
    std::vector<ad_t> h0s(nS), bra(nS);
    for (std::size_t s = 0; s < nS; ++s) {
      h0s[s] = ad_t(F.h0[s]); bra[s] = ad_t(F.birth_rate[s]);
    }
    bra[tgt] = br_t;                                   // only the target's birth_rate active
    std::vector<ad_t> J = assemble_metrics_coupled_ms<ad_t>(pas, F, metrics, h0s, 0.0,
                            nullptr, nullptr, &bra);
    for (std::size_t m = 0; m < M; ++m) tape.registerOutput(J[m]);
    for (std::size_t m = 0; m < M; ++m) {
      tape.clearDerivatives();
      xad::derivative(J[m]) = 1.0;
      tape.computeAdjoints();
      dbr[m] = xad::derivative(br_t);
      values[m] = as_double(J[m]);
    }
  }
  dbr.attr("names") = Rcpp::wrap(metrics);
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("d_birth_rate") = dbr,
                            Rcpp::Named("values")       = values,
                            Rcpp::Named("env_err")      = env_err);
}

// FULLY native cross-species resident birth_rate gradient (mirror of
// ff16_coupled_gradient_ms_native): joint env + schedule + birth steps + all-species
// boundary harvest from the live Patch. birth_rate the per-species driver vector; target
// 1-based. Returns d(total-stand metric)/d(birth_rate_target) + values + env_err.
// [[Rcpp::export]]
Rcpp::List ff16_birth_rate_gradient_ms_native(
    SEXP scm_, Rcpp::List pp_list, std::vector<std::string> metrics,
    std::vector<double> birth_rate, double patch_area, int target,
    bool active_birthenv = true) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::FF16_Strategy, plant::FF16_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  if (patch.stand_newnode_height_stage_history_all.size() < 1)
    Rcpp::stop("the multi-species resident birth_rate gradient needs the all-species "
               "per-RK-stage harvest; re-run the resident SCM with "
               "control(save_RK45_cache = TRUE)");
  FrozenMS F = build_frozen_ms_scm(pp_list, patch, birth_rate, patch_area, active_birthenv);
  const std::size_t tgt = (std::size_t)(target - 1);
  if (tgt >= F.nS) Rcpp::stop("target species out of range");
  Rcpp::NumericVector pp0 = pp_list[(R_xlen_t)tgt];
  plant::FF16_Strategy s0 = make_strategy(pp0);   // owns the GK integrator
  F.integ = &s0.function_integrator;
  return ff16_birth_rate_gradient_ms_core(F, metrics, target);
}

// R0 gate for the MULTI-SPECIES coupled replay (#472 scope B, R2): a double-precision
// joint whole-stand re-evolution reconstructing the joint canopy per RK stage. Returns
// the reconstructed TOTAL-stand metric VALUES and `env_err` = worst |reconstructed knot
// light - SCM knot light| over all stages (the coupled-drift gauge). pp_list/birth_list
// are per-species; nn_h_list/nn_c_list the all-species boundary harvest [step][stage][species].
// [[Rcpp::export]]
Rcpp::List ff16_coupled_metrics_ms_impl(
    Rcpp::List pp_list, Rcpp::List eh_list, std::vector<double> sh,
    Rcpp::List birth_list, std::vector<std::string> metrics,
    std::vector<double> birth_rate, Rcpp::List nn_h_list, Rcpp::List nn_c_list,
    double patch_area, bool active_birthenv = true) {
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("coupled MS metrics: expected LAI / biomass / size_moment, got " + nm);
  FrozenMS F = build_frozen_ms(pp_list, eh_list, sh, birth_list, birth_rate,
                               nn_h_list, nn_c_list, patch_area, active_birthenv);
  // Own a strategy so its (strategy-independent) GK integrator outlives the assemble.
  Rcpp::NumericVector pp0 = pp_list[0];
  plant::FF16_Strategy s0 = make_strategy(pp0);
  F.integ = &s0.function_integrator;
  std::vector<double> h0s = F.h0;
  double env_err = 0.0, env_err_z = -1.0;
  std::vector<double> J = assemble_metrics_coupled_ms<double>(F.pd, F, metrics, h0s,
                           0.0, &env_err, &env_err_z);
  Rcpp::NumericVector values(J.size());
  for (std::size_t m = 0; m < J.size(); ++m) values[m] = J[m];
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("values") = values,
                            Rcpp::Named("env_err") = env_err,
                            Rcpp::Named("env_err_z") = env_err_z);
}

// R1 of the MULTI-SPECIES coupled replay (#472 scope B, R2): the CROSS-SPECIES resident
// total trait gradient of the TOTAL-stand census metrics, w.r.t. the traits of ONE
// species (`target` is 1-based). All species' cohorts are lifted to the adjoint tape and
// re-evolved jointly; only the target species' trait fields (+ its seedling-size IFT
// channel) are registered as inputs, so one reverse sweep per metric yields
// d(total metric)/d(theta of the target species) -- including the cross term whereby the
// target's traits re-shade the joint canopy that every species reads. nn_h_list/nn_c_list:
// the all-species boundary harvest [step][stage][species].
// Shared core of the MS coupled-gradient entries: F already carries the joint harvest +
// integrator. Runs the R0 double pass for the env-drift gate (env_err = worst joint-knot
// light drift -- a stiff schedule diverges and the gradient is meaningless), then the AD
// tape for the cross-species gradient. The R-list (`_impl`) and native-SCM (`_native`)
// entries differ only in how F is built (lossy eh_list + R harvest vs faithful Patch read).
Rcpp::List ff16_coupled_gradient_ms_core(FrozenMS& F, std::vector<std::string> traits,
    std::vector<std::string> metrics, int target) {
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("coupled MS gradient: expected LAI / biomass / size_moment, got " + nm);
  const std::size_t tgt = (std::size_t)(target - 1);
  if (tgt >= F.nS) Rcpp::stop("target species out of range");
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = F.h0[tgt];
  std::vector<double> dh0 = compute_dh0(F.pd[tgt], h0v, idx);
  const std::size_t M = metrics.size(), nT = idx.size(), nS = F.nS;

  // R0 double pass: the joint re-evolution must reproduce the resident stand (env_err
  // small) before the gradient is trusted; the caller gates on it.
  double env_err = 0.0, env_err_z = -1.0;
  (void)assemble_metrics_coupled_ms<double>(F.pd, F, metrics, F.h0, 0.0,
                                            &env_err, &env_err_z);

  Rcpp::NumericMatrix jac(M, nT);
  Rcpp::NumericVector values(M);
  {
    // One joint adjoint tape, scoped. Lift every species' prod_pars to ad_t (the canopy
    // couples them) but register ONLY the target species' selected traits as inputs.
    ad::tape_type tape;
    std::vector<plant::FF16ProdPars<ad_t>> pas(nS);
    for (std::size_t s = 0; s < nS; ++s) pas[s] = lift<ad_t>(F.pd[s]);
    auto fp = field_ptrs<ad_t>(pas[tgt]);
    for (auto i : idx) tape.registerInput(*fp[i]);
    tape.newRecording();
    // Per-species seedling height: active (IFT injection) for the target, constant else.
    std::vector<ad_t> h0s(nS);
    for (std::size_t s = 0; s < nS; ++s) h0s[s] = ad_t(F.h0[s]);
    h0s[tgt] = inject_h0(h0v, fp, idx, dh0);
    std::vector<ad_t> J = assemble_metrics_coupled_ms<ad_t>(pas, F, metrics, h0s, 0.0);
    for (std::size_t m = 0; m < M; ++m) tape.registerOutput(J[m]);
    for (std::size_t m = 0; m < M; ++m) {
      tape.clearDerivatives();
      xad::derivative(J[m]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < nT; ++k) jac(m, k) = xad::derivative(*fp[idx[k]]);
      values[m] = as_double(J[m]);
    }
  }
  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac,
                            Rcpp::Named("values")   = values,
                            Rcpp::Named("env_err")  = env_err);
}

// [[Rcpp::export]]
Rcpp::List ff16_coupled_gradient_ms_impl(
    Rcpp::List pp_list, Rcpp::List eh_list, std::vector<double> sh,
    Rcpp::List birth_list, std::vector<std::string> traits,
    std::vector<std::string> metrics, std::vector<double> birth_rate,
    Rcpp::List nn_h_list, Rcpp::List nn_c_list, double patch_area, int target,
    bool active_birthenv = true) {
  FrozenMS F = build_frozen_ms(pp_list, eh_list, sh, birth_list, birth_rate,
                               nn_h_list, nn_c_list, patch_area, active_birthenv);
  const std::size_t tgt = (std::size_t)(target - 1);
  if (tgt >= F.nS) Rcpp::stop("target species out of range");
  Rcpp::NumericVector pp0 = pp_list[(R_xlen_t)tgt];
  plant::FF16_Strategy s0 = make_strategy(pp0);   // owns the GK integrator
  F.integ = &s0.function_integrator;
  return ff16_coupled_gradient_ms_core(F, traits, metrics, target);
}

// FULLY native multi-species coupled cross-species gradient: the joint env + schedule +
// birth steps + all-species boundary harvest come from the live Patch (no R
// ff16_harvest_ms, no Rcpp::as<> env). pp_list is cheap ($parameters$strategies);
// birth_rate the per-species constant driver; target 1-based. Returns the cross-species
// Jacobian + values + env_err (the caller gates the stiff-schedule case on env_err).
// [[Rcpp::export]]
Rcpp::List ff16_coupled_gradient_ms_native(
    SEXP scm_, Rcpp::List pp_list, std::vector<std::string> traits,
    std::vector<std::string> metrics, std::vector<double> birth_rate,
    double patch_area, int target, bool active_birthenv = true) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::FF16_Strategy, plant::FF16_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  if (patch.stand_newnode_height_stage_history_all.size() < 1)
    Rcpp::stop("feedback = 'resident' on a multi-species stand needs the all-species "
               "per-RK-stage harvest; re-run the resident SCM with "
               "control(save_RK45_cache = TRUE)");
  FrozenMS F = build_frozen_ms_scm(pp_list, patch, birth_rate, patch_area, active_birthenv);
  const std::size_t tgt = (std::size_t)(target - 1);
  if (tgt >= F.nS) Rcpp::stop("target species out of range");
  Rcpp::NumericVector pp0 = pp_list[(R_xlen_t)tgt];
  plant::FF16_Strategy s0 = make_strategy(pp0);   // owns the GK integrator
  F.integ = &s0.function_integrator;
  return ff16_coupled_gradient_ms_core(F, traits, metrics, target);
}

// Escape hatch (#472 scope B, build-order step 1): the per-cohort state x trait
// Jacobian. For ANY emergent metric that is NOT a simple reduction -- quantiles,
// ratios, bespoke statistics a downstream package invents -- expose d(state_i,c)/
// d(theta_k) for every cohort i, every final-state component c, and let the metric
// gradient compose downstream by chain rule (the same boundary as "likelihoods
// live downstream"). Each cohort's final state is INDEPENDENT (it depends only on
// its own replay), so we tape ONE cohort at a time and take one reverse sweep per
// state component -- nC small sweeps, not one giant tape. Returns the cohort final
// states + the [cohort, component, trait] Jacobian.
// [[Rcpp::export]]
Rcpp::List ff16_state_jacobian_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits) {
  auto s  = make_strategy(pp);
  auto pd = s.prod_pars();
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);
  const std::size_t nC = F.birth.size(), nT = idx.size();

  const std::vector<std::string> comp =
    {"height","mortality","fecundity","area_heartwood","mass_heartwood","offspring"};
  const std::size_t nS = comp.size();

  Rcpp::NumericMatrix states(nC, nS);
  // 3D Jacobian, R column-major dim [nC, nS, nT].
  Rcpp::NumericVector jac(nC * nS * nT);
  auto JAC = [&](std::size_t i, std::size_t c, std::size_t k) -> double& {
    return jac[i + nC * (c + nS * k)];
  };

  for (std::size_t i = 0; i < nC; ++i) {
    ad::tape_type tape;
    auto pa = lift<ad_t>(pd);
    auto fp = field_ptrs<ad_t>(pa);
    for (auto j : idx) tape.registerInput(*fp[j]);
    tape.newRecording();
    ad_t h0 = h0v;
    for (std::size_t k = 0; k < nT; ++k)
      h0 = h0 + ad_t(dh0[k]) * (*fp[idx[k]] - ad_t(xad::value(*fp[idx[k]])));

    plant::FF16LifeState<ad_t> y = replay_cohort_final<ad_t>(pa, F, i, h0);
    ad_t out[6] = {y.demog.height, y.demog.mortality, y.demog.fecundity,
                   y.demog.area_heartwood, y.demog.mass_heartwood, y.offspring};
    for (std::size_t c = 0; c < nS; ++c) {
      states(i, c) = as_double(out[c]);
      tape.registerOutput(out[c]);
    }
    for (std::size_t c = 0; c < nS; ++c) {
      tape.clearDerivatives();
      xad::derivative(out[c]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < nT; ++k) JAC(i, c, k) = xad::derivative(*fp[idx[k]]);
    }
  }

  states.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp));
  jac.attr("dim") = Rcpp::IntegerVector::create((int)nC, (int)nS, (int)nT);
  jac.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp),
                                            Rcpp::wrap(traits));
  return Rcpp::List::create(Rcpp::Named("states") = states,
                            Rcpp::Named("jacobian") = jac);
}

// Trait gradient of grow_individual_to_size (FF16): given the harvested adaptive step
// schedule `sh` (step times t[0..M] from grow_individual_bracket), the initial ODE
// state `y0v`, the FIXED environment `env`, and a vector of `targets` for the size
// component `sidx` (0=height,...), return per (size x trait) the derivative of the
// stopping time t* and per (size x component x trait) the TOTAL derivative of the ODE
// state at t*. active_h0 injects the seedling-size IFT d(h0)/d(theta) (matching the
// live grow path, which starts each plant at its trait-dependent initial_height()).
// [[Rcpp::export]]
Rcpp::List ff16_grow_to_size_gradient_impl(
    Rcpp::NumericVector pp, plant::FF16_Environment env, Rcpp::NumericVector y0v,
    std::vector<double> sh, std::vector<double> targets, int sidx,
    std::vector<std::string> traits, bool active_h0) {
  auto s = make_strategy(pp);
  auto pd = s.prod_pars();
  const double eta = s.pars.eta;
  const plant::quadrature::QK* integ = &s.function_integrator;
  const plant::FF16_Environment* e = &env;
  std::vector<std::size_t> idx = resolve_traits(traits);
  const std::size_t nT = idx.size();
  if (sh.size() < 2) Rcpp::stop("schedule must have at least two step times");
  std::vector<double> step_h(sh.size() - 1);
  for (std::size_t n = 0; n + 1 < sh.size(); ++n) step_h[n] = sh[n + 1] - sh[n];

  const double h0v = y0v[0];
  std::vector<double> dh0 = active_h0 ? compute_dh0(pd, h0v, idx) : std::vector<double>(nT, 0.0);
  plant::FF16State<double> y0{y0v[0], y0v[1], y0v[2], y0v[3], y0v[4]};

  const std::vector<std::string> comp =
    {"height","mortality","fecundity","area_heartwood","mass_heartwood"};
  const std::size_t nS = comp.size(), nG = targets.size();

  Rcpp::NumericVector tstar(nG);
  Rcpp::NumericMatrix state(nG, nS);
  Rcpp::NumericMatrix dtime(nG, nT);                  // d(t*)/d(theta)
  Rcpp::NumericVector dstate(nG * nS * nT);           // [nG, nS, nT] total d y_c/d theta
  auto DS = [&](std::size_t g, std::size_t c, std::size_t k) -> double& {
    return dstate[g + nG * (c + nS * k)]; };

  for (std::size_t g = 0; g < nG; ++g) {
    int nfull; double dt_final; plant::FF16State<double> fin;
    grow_discover(pd, integ, eta, e, y0, step_h, targets[g], sidx, nfull, dt_final, fin);
    if (nfull < 0) Rcpp::stop("target size not reached within the schedule");
    tstar[g] = sh[nfull] + dt_final;
    double fs[5] = {fin.height, fin.mortality, fin.fecundity, fin.area_heartwood, fin.mass_heartwood};
    for (std::size_t c = 0; c < nS; ++c) state(g, c) = fs[c];
    // Rates at t* (double) for the IFT correction.
    plant::FF16State<double> rate = grow_deriv_fixed<double>(pd, integ, eta, e, fin);
    double yd[5] = {rate.height, rate.mortality, rate.fecundity, rate.area_heartwood, rate.mass_heartwood};
    const double sdot = yd[sidx];

    // Frozen schedule to t*: nfull full steps + the single partial step dt_final.
    std::vector<double> sched(step_h.begin(), step_h.begin() + nfull);
    sched.push_back(dt_final);

    // AD: partial d y_c/d theta holding t* fixed (one reverse sweep per component).
    ad::tape_type tape;
    auto pa = lift<ad_t>(pd);
    auto fp = field_ptrs<ad_t>(pa);
    for (auto j : idx) tape.registerInput(*fp[j]);
    tape.newRecording();
    ad_t h0 = h0v;
    for (std::size_t k = 0; k < nT; ++k)
      h0 = h0 + ad_t(dh0[k]) * (*fp[idx[k]] - ad_t(xad::value(*fp[idx[k]])));
    plant::FF16State<ad_t> yad{h0, ad_t(y0v[1]), ad_t(y0v[2]), ad_t(y0v[3]), ad_t(y0v[4])};
    plant::FF16State<ad_t> out = grow_replay_fixed<ad_t>(pa, integ, eta, e, yad, sched);
    ad_t oc[5] = {out.height, out.mortality, out.fecundity, out.area_heartwood, out.mass_heartwood};
    for (std::size_t c = 0; c < nS; ++c) tape.registerOutput(oc[c]);
    std::vector<std::vector<double>> P(nS, std::vector<double>(nT, 0.0));
    for (std::size_t c = 0; c < nS; ++c) {
      tape.clearDerivatives();
      xad::derivative(oc[c]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < nT; ++k) P[c][k] = xad::derivative(*fp[idx[k]]);
    }
    // IFT: d(t*)/d(theta_k) = -(d size/d theta_k|t*) / size_dt(t*); total = P + ydot*dt*.
    for (std::size_t k = 0; k < nT; ++k) {
      double dtk = (sdot != 0.0) ? -P[sidx][k] / sdot : 0.0;
      dtime(g, k) = dtk;
      for (std::size_t c = 0; c < nS; ++c) DS(g, c, k) = P[c][k] + yd[c] * dtk;
    }
  }
  state.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp));
  dtime.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(traits));
  dstate.attr("dim") = Rcpp::IntegerVector::create((int)nG, (int)nS, (int)nT);
  dstate.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp), Rcpp::wrap(traits));
  return Rcpp::List::create(Rcpp::Named("time") = tstar, Rcpp::Named("state") = state,
                            Rcpp::Named("d_time") = dtime, Rcpp::Named("d_state") = dstate);
}
