// Census-metric reconstruction for TF24f (#472 scope B, build-order step 1 -- the
// TF24f frozen-census R0 GATE). No AD: a pure double-precision replay that proves
// the collar-state trajectory + the census number-density (log_density) reproduce
// the SCM's stored stand before any reverse-mode tape is built (the R1 follow-up).
//
// Why TF24f, not TF24, gets the census gradient (the scope decision): TF24's census
// number density needs the
// SECOND-order leaf-optimiser sensitivity d(growth-rate-gradient)/d(theta), which the
// linearised leaf-opt harvest of tf24_emergent.cpp does not differentiate faithfully.
// TF24f removes the per-step golden-section optimiser entirely: the optimal root-collar
// potential is a 6th ODE state (opt_root_psi_state) that relaxes toward its optimum by
// gradient ascent (dpsi/dt = k_acclim * dprofit_dpsi), and the leaf is evaluated AT the
// tracked collar by an analytic, IFT-able path (Leaf::dprofit_droot_collar_psi). So the
// trajectory is a clean 6-state ODE with differentiable rates and the census g' needs no
// curvature harvest -- carrying the collar as a replayed state propagates the height
// sensitivity directly (the R1 tape).
//
// This R0 gate establishes the prerequisite the scope note flags first (§7): that the
// collar-state replay is FAITHFUL. It re-evolves, per cohort, the 5 demographic states +
// the tracked collar + log_density over the SCM's frozen Cash-Karp schedule, reading the
// frozen resident environment and driving the REAL TF24f leaf evaluation at the tracked
// collar (crown-centre: one analytic leaf solve per stage, no optimiser). The census
// number density uses the SCM's own growth-rate-gradient scheme exactly: a backward
// finite difference of height_dt with step node_gradient_eps = 1e-6, holding the collar
// fixed and perturbing only height (Node::growth_rate_gradient with the default
// node_gradient_direction = -1, exact_ad off). It returns the reconstructed per-cohort
// heights / collar states / log-densities and the census metric values (LAI / biomass /
// size_moment); the R gate checks these against the SCM's stored stand and
// compute_competition(0).
#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <XAD/XAD.hpp>             // reverse-mode adjoint tape (the R1 census gradient)
#include <plant.h>                 // RcppR6 as<>/wrap for TF24_Environment
#include <plant/models/tf24f_strategy.h>
#include <plant/models/tf24_production_kernel.h>
#include <plant/models/ff16_production_kernel.h>   // ff16_cashkarp_replay, FF16State
#include <plant/gradient/scm_harvest.h>            // shared birth_steps / recover_birth_rate / census_trapezium
#include <plant/gradient/coupled_canopy.h>         // shared Yokozawa canopy_comp_at
#include <odelia/interpolator.hpp>                 // basic_interpolator<S> (coupled canopy)

namespace {

// Build a TF24f strategy from the SCM's parameter vector. Mirrors the TF24
// make_strategy in tf24_emergent.cpp but constructs the fast-acclimation variant and
// sets its acclimation knobs. crown-centre + GSS_tol_abs match the resident run.
plant::TF24f_Strategy make_tf24f(const Rcpp::NumericVector& pp, double k_acclim,
                                 bool use_ad_gradient,
                                 const std::string& shading = "crown-centre",
                                 double gss_tol = 1e-9) {
  plant::TF24f_Strategy s;
  s.control.shading_model = shading; s.control.GSS_tol_abs = gss_tol;
  s.k_acclim = k_acclim; s.use_ad_gradient = use_ad_gradient;
  auto& q = s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.a_d0=pp["a_d0"];
  q.recruitment_decay=pp["recruitment_decay"];
  q.vcmax_25=pp["vcmax_25"];q.K_s=pp["K_s"];q.b=pp["b"];q.c=pp["c"];q.beta2=pp["beta2"];
  q.jmax_25=pp["jmax_25"];q.a=pp["a"];q.curv_fact_elec_trans=pp["curv_fact_elec_trans"];
  q.curv_fact_colim=pp["curv_fact_colim"];
  s.prepare_strategy();
  return s;
}

// The 7-state cohort: TF24's 5 demographic states + the tracked root-collar potential
// (TF24f's 6th ODE state) + the census number density log_density.
struct St7 { plant::FF16State<double> demog; double collar; double log_density; };

// Evaluate the leaf at a GIVEN tracked collar potential, returning net production and
// (out) the profit gradient dprofit/dpsi that drives the collar's gradient-ascent rate.
// Drives the real TF24f leaf path: setting tracked_root_psi_ then calling the inherited
// net_mass_production_dt makes solve_leaf() evaluate at the tracked collar (clamped to
// the feasible interval) rather than re-optimise -- the analytic, IFT-able TF24f path.
// crown-centre runs exactly one leaf solve, so dprofit_dpsi_ is set cleanly.
double net_at_tracked(plant::TF24f_Strategy& s, const plant::TF24_Environment& e,
                      double h, double collar, double& dprofit_dpsi) {
  s.tracked_root_psi_ = collar;
  const double net = s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0 / h);
  dprofit_dpsi = s.dprofit_dpsi_;
  return net;
}

} // namespace



// Compiled core of tf24f_census_recon(). Re-evolves every cohort's {5 demog, collar,
// log_density} over the harvested frozen schedule (per-RK-stage TF24_Environment in
// eh_list, step times sh, per-cohort birth steps), driving the real TF24f leaf at the
// tracked collar in double precision. Returns the reconstructed per-cohort final heights
// / collar states / log-densities and the census metric values (the height-trapezium of
// density*kI*area_leaf etc., matching Patch::compute_competition / FF16's census reduce).
// Shared core: takes the per-RK-stage env as a NATIVE vector (EH) so the R-list
// (`_impl`) and live-SCM (`_native`) entries differ only in how EH is sourced
// (lossy Rcpp::as<> vs faithful native copy). EH[n][s] is the TF24_Environment at
// step n, Cash-Karp stage s; sh the step times.
Rcpp::List tf24f_census_recon_core(
    Rcpp::NumericVector pp,
    const std::vector<std::vector<plant::TF24_Environment>>& EH,
    std::vector<double> sh, std::vector<int> birth, double birth_rate,
    double k_acclim, bool use_ad_gradient, std::vector<std::string> metrics,
    bool exact_ad_gprime) {
  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double kI = s.pars.k_I, recruitment_decay = s.pars.recruitment_decay;
  const double GEPS = 1e-6;                     // Control::node_gradient_eps (backward FD)

  // Validate metrics up front.
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm +
                 " (LAI / biomass / size_moment)");

  const std::size_t N = EH.size(), nC = birth.size();
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];

  // Replay one cohort to its final state over the frozen schedule.
  auto replay_cohort = [&](std::size_t i) -> St7 {
    using std::exp; using std::log;
    const std::size_t b = (std::size_t)birth[i];
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];

    // --- Seed (Node::compute_initial_conditions) -----------------------------------
    // 1. collar_0 = the optimum at the birth env (set_initial_states runs the base
    //    golden-section optimiser ONCE via initializing_).
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    const double collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    // 2. Establishment + seedling growth rate, evaluated at the tracked collar_0 (the
    //    operating point the SCM uses: compute_rates set tracked_root_psi_ = collar_0
    //    before establishment_probability).
    const double decay = std::exp(-recruitment_decay * sh[b]);
    double dpsi0 = 0.0;
    const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, eb, h0v, collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double mort0 = (pr_estab > 0.0) ? -log(pr_estab)
                                          : std::numeric_limits<double>::infinity();
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    const double logd0 = (g0 > 0.0 && pr_estab > 0.0)
                           ? log(birth_rate * pr_estab / g0)
                           : -std::numeric_limits<double>::infinity();

    // --- Cash-Karp replay over the frozen schedule ---------------------------------
    auto deriv = [&](const St7& y, std::size_t n, int stage) -> St7 {
      const plant::TF24_Environment& e =
        (stage == 0) ? ((n > 0) ? EH[n - 1][5] : EH[0][0]) : EH[n][stage - 1];
      const double h = y.demog.height, collar = y.collar;
      double dpsi = 0.0;
      const double net = net_at_tracked(s, e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      const double collar_dt = k_acclim * dpsi;
      // g' matches the SCM's Node::growth_rate_gradient, collar held fixed. Either the
      // exact-AD gradient (node_gradient_exact_ad = TRUE, the analytic kmax+light+E_up
      // channels) or the backward FD (dir = -1, step GEPS) -- whichever the resident
      // run used -- so the reconstructed density reproduces the SCM's stored density.
      double gprime;
      if (exact_ad_gprime) {
        s.tracked_root_psi_ = collar;
        gprime = s.growth_rate_gradient_height_ad(h, e);
        // restore the operating point at (h, collar) for the rate fill below
        net_at_tracked(s, e, h, collar, dpsi);
      } else {
        double dpsi_b = 0.0;
        const double net_b = net_at_tracked(s, e, h - GEPS, collar, dpsi_b);
        const double al_b = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h - GEPS);
        const double g_back = plant::tf24_height_dt_from_net<double>(pd, h - GEPS, al_b, net_b);
        gprime = (r.height_dt - g_back) / GEPS;
      }
      const double log_density_dt = -gprime - r.mortality_dt;
      return St7{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, log_density_dt};
    };
    auto axpy = [](const St7& a, double c, const St7& k) -> St7 {
      return St7{plant::FF16State<double>{
        a.demog.height + c * k.demog.height, a.demog.mortality + c * k.demog.mortality,
        a.demog.fecundity + c * k.demog.fecundity,
        a.demog.area_heartwood + c * k.demog.area_heartwood,
        a.demog.mass_heartwood + c * k.demog.mass_heartwood},
        a.collar + c * k.collar, a.log_density + c * k.log_density};
    };
    St7 y{plant::FF16State<double>{h0v, mort0, 0.0, 0.0, 0.0}, collar0, logd0};
    return plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
  };

  std::vector<St7> finals(nC);
  for (std::size_t i = 0; i < nC; ++i) finals[i] = replay_cohort(i);

  // The pending-seed (new_node) census tail term: a seed born into the FINAL-time env
  // (eh.back()[5]), giving density_new for the trapezium tail from h_last down to h0.
  double dens_new = 0.0;
  {
    const plant::TF24_Environment& ef = EH[N - 1][5];
    // collar seed at the final env
    s.initializing_ = true;
    s.net_mass_production_dt(ef, h0v, s.area_leaf(h0v), 1.0 / h0v);
    const double collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    const double decay = std::exp(-recruitment_decay * sh[N]);
    double dpsi0 = 0.0;
    const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, ef, h0v, collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    dens_new = (g0 > 0.0 && pr_estab > 0.0) ? birth_rate * pr_estab / g0 : 0.0;
  }

  // --- Census reduction (height-trapezium, descending height + pending-seed tail) ---
  // Mirrors ff16_emergent.cpp's census_reduce / Species::compute_competition exactly.
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return finals[a].demog.height > finals[b].demog.height; });

  auto census_reduce = [&](auto psi) -> double {
    return plant::gradient::census_trapezium<double>(nC, ord, h0v, dens_new,
      [&](std::size_t i) -> double { return finals[i].demog.height; },
      [&](std::size_t i) -> double { return std::exp(finals[i].log_density); },
      [&](std::size_t i) -> double { return finals[i].demog.mass_heartwood; }, psi);
  };

  Rcpp::NumericVector values(metrics.size());
  for (std::size_t m = 0; m < metrics.size(); ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      values[m] = census_reduce([&](double h, double dens, double) -> double {
        return dens * kI * plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h); });
    } else if (nm == "biomass") {
      // live mass + heartwood (mirrors FF16 biomass); TF24 mass_live via the cascade.
      values[m] = census_reduce([&](double h, double dens, double mhw) -> double {
        const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
        const double mass_leaf    = al * pd.lma;
        const double area_sapwood = al * pd.theta;
        const double mass_sapwood = area_sapwood * h * pd.eta_c * pd.rho;
        const double area_bark    = pd.a_b1 * al * pd.theta;
        const double mass_bark    = area_bark * h * pd.eta_c * pd.rho;
        const double mass_root    = pd.a_r1 * al;
        return dens * (mass_leaf + mass_sapwood + mass_bark + mass_root + mhw); });
    } else { // size_moment
      values[m] = census_reduce([&](double h, double dens, double) -> double {
        return dens * h; });
    }
  }
  values.attr("names") = Rcpp::wrap(metrics);

  Rcpp::NumericVector heights(nC), collar(nC), log_density(nC);
  for (std::size_t i = 0; i < nC; ++i) {
    heights[i]     = finals[i].demog.height;
    collar[i]      = finals[i].collar;
    log_density[i] = finals[i].log_density;
  }
  return Rcpp::List::create(
    Rcpp::Named("values") = values,
    Rcpp::Named("heights") = heights,
    Rcpp::Named("collar") = collar,
    Rcpp::Named("log_density") = log_density,
    Rcpp::Named("density_new") = dens_new);
}

namespace {
// Lossy: rebuild the per-RK-stage env from an R list (Rcpp::as<TF24_Environment>
// reconstructs each light spline from its serialized knots via r_init_interpolators
// -- NOT faithful for the crown-sampled light, the TF24f long-horizon census floor).
std::vector<std::vector<plant::TF24_Environment>> eh_from_list(Rcpp::List eh_list) {
  const std::size_t N = eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));
  }
  return EH;
}
// Per-cohort birth steps from a live TF24f patch: the step time nearest each node
// introduction time (0-based), matching tf24f_harvest's which.min. Computing this in
// C++ lets the census native entries skip the R ff16/tf24f_harvest entirely (whose
// O(stand) ppsurv loop is dead work for the census metrics).
std::vector<int> tf24f_birth_steps(
    const plant::Patch<plant::TF24f_Strategy, plant::TF24_Environment>& patch,
    std::size_t species) {
  return plant::gradient::birth_steps(patch, species);
}
} // namespace

// Native-SCM recon: faithful per-RK-stage env read directly from the live Patch.
// scm_ is the RcppR6 SCM<TF24f,TF24_Env>; the wrapper + patch ref stay alive for the
// (synchronous) core call, so the borrowed env history never dangles.
// [[Rcpp::export]]
Rcpp::List tf24f_census_recon_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, double birth_rate,
    double k_acclim, bool use_ad_gradient, std::vector<std::string> metrics,
    bool exact_ad_gprime) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::TF24f_Strategy, plant::TF24_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  const std::vector<int> birth = tf24f_birth_steps(patch, (std::size_t)species);
  return tf24f_census_recon_core(pp, patch.environment_history, patch.step_history,
                                 birth, br, k_acclim, use_ad_gradient,
                                 metrics, exact_ad_gprime);
}


// ===========================================================================
// TF24f frozen-census reverse-mode AD trait gradient (#472 scope B, build-order
// step 2 -- the R1 tape, the "refine" step). One reverse sweep over the 7-state
// replay {5 demog, tracked collar, log_density} per metric, replacing the
// per-trait FD over the recon (tf24f_census_gradient_fd) it must reproduce.
//
// The structure extends the TF24 OFFSPRING tape (src/tf24_emergent.cpp): a
// double DISCOVERY pass runs the real TF24f leaf at the tracked collar along the
// frozen-schedule trajectory and HARVESTS the trait-independent operating point
// per RK stage; a second pass replays a leaf-opt-free, fully tapeable expression.
//
// The one ingredient beyond the offspring tape:
// TF24f's tracked collar is a strongly theta-dependent STATE that at k_acclim=1
// LAGS the optimum, so the envelope theorem does NOT zero its contribution. We
// therefore (a) carry the collar as a TAPED state whose gradient-ascent rate
// k_acclim*dprofit_dpsi is linearised with a CURVATURE harvest (the second
// derivatives d2profit/dpsi2, d2profit/dpsi dh, d2profit/dpsi dtheta_k), and
// (b) linearise profit itself in (h, collar, theta) -- the collar term
// dprofit_dpsi0*(collar-collar0) is what the offspring tape omits. The census
// number density's g' = d(height_dt)/d(height) reproduces the SCM's backward-FD
// scheme by harvesting a SECOND operating point at h-GEPS and differencing the
// two tapeable height_dt expressions (as ff16_emergent.cpp's census tape does).
//
// All harvested derivatives are taken by finite difference in the double pass
// (no leaf templating): profit/dprofit_dpsi are differenced w.r.t. collar, height
// and -- via independently-built perturbed strategies -- each requested trait.
// The trait set is a prototype subset (the seed: validate a few traits, not all
// 27); every requested trait must name a strategy `pars` entry (matching the FD
// gate, which perturbs the same parameter vector).
// ===========================================================================

namespace {

using ad   = xad::adj<double>;
using ad_t = ad::active_type;

double as_dbl(const ad_t& v)  { return xad::value(v); }

// Evaluate the real TF24f leaf at a GIVEN (height, tracked collar) in env e,
// returning the leaf profit and (out) the analytic dprofit/dpsi the acclimation
// rate uses. Drives the same path as net_at_tracked / the recon, but exposes
// leaf.profit_ (the quantity the tape linearises and feeds to the kernel).
void eval_leaf(plant::TF24f_Strategy& s, const plant::TF24_Environment& e,
               double h, double collar, double& profit, double& dpsi) {
  s.tracked_root_psi_ = collar;
  s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0 / h);
  profit = s.leaf.profit_;
  dpsi   = s.dprofit_dpsi_;
}

// Trait-independent leaf-operating-point harvest at one (h, collar) point: the
// profit value, its first derivatives w.r.t. height (collar held) and collar, the
// per-trait profit sensitivity inj_k (the leaf channel theta -> profit, the TF24
// `inj` term done by FD over a perturbed strategy), and the curvature the collar
// gradient-ascent rate needs: d(dprofit_dpsi)/d{collar, height, theta_k}.
struct LH {
  double h, collar;                 // the harvested operating point
  double profit, dprofit_dh, dprofit_dpsi;
  double d2p_dpsi2, d2p_dpsidh;     // collar-rate curvature (collar / height)
  std::vector<double> inj;          // d(profit)/d(theta_k)            [per trait]
  std::vector<double> d2p_dpsidth;  // d(dprofit_dpsi)/d(theta_k)      [per trait]
};

// Harvest at (h, collar) in env e. s is the base strategy; sp[k]/sm[k] are the
// +/- trait-perturbed strategies (built once) and dk[k] the matching abs steps.
LH harvest_point(plant::TF24f_Strategy& s,
                 std::vector<plant::TF24f_Strategy>& sp,
                 std::vector<plant::TF24f_Strategy>& sm,
                 const std::vector<double>& dk,
                 const plant::TF24_Environment& e, double h, double collar) {
  const std::size_t T = sp.size();
  LH H; H.h = h; H.collar = collar; H.inj.resize(T); H.d2p_dpsidth.resize(T);
  double dpsi0;
  eval_leaf(s, e, h, collar, H.profit, dpsi0);
  H.dprofit_dpsi = dpsi0;

  // height channel (collar held): central FD of profit + dprofit_dpsi.
  const double dh = 1e-5 * h;
  double pph, dph, pmh, dmh;
  eval_leaf(s, e, h + dh, collar, pph, dph);
  eval_leaf(s, e, h - dh, collar, pmh, dmh);
  H.dprofit_dh = (pph - pmh) / (2.0 * dh);
  H.d2p_dpsidh = (dph - dmh) / (2.0 * dh);

  // collar channel (height held): d2profit/dpsi2 from FD of dprofit_dpsi.
  const double dc = 1e-5 * std::abs(collar) + 1e-8;
  double ppc, dpc, pmc, dmc;
  eval_leaf(s, e, h, collar + dc, ppc, dpc);
  eval_leaf(s, e, h, collar - dc, pmc, dmc);
  H.d2p_dpsi2 = (dpc - dmc) / (2.0 * dc);

  // trait channels: inj_k = d(profit)/d(theta_k); d2p_dpsidth_k = d(dpsi)/d(theta_k),
  // both at the SAME (h, collar) via the prebuilt perturbed strategies.
  for (std::size_t k = 0; k < T; ++k) {
    double pp_, dpp, pm_, dpm;
    eval_leaf(sp[k], e, h, collar, pp_, dpp);
    eval_leaf(sm[k], e, h, collar, pm_, dpm);
    H.inj[k]         = (pp_ - pm_) / (2.0 * dk[k]);
    H.d2p_dpsidth[k] = (dpp - dpm) / (2.0 * dk[k]);
  }
  return H;
}

// Per-stage harvest: the forward operating point (h, collar) and the backward
// point (h - GEPS, collar) the census g' (backward FD of height_dt) reads.
struct StageH { LH fwd, back; };

// Linearised leaf profit at (active height h, active collar) about a harvested
// point Hp, with the requested traits active: profit0 + dprofit_dh*(h - h0)
// + dprofit_dpsi*(collar - collar0) + sum_k inj_k*(theta_k - theta_k0). The
// collar term is the offspring tape's omission (collar is a taped state here).
// h is the active evaluation height; Hp.h the height the point was harvested at
// (for the backward point both already carry the -GEPS offset, so the difference
// h - Hp.h is the active displacement from the harvested trajectory either way).
ad_t profit_lin(const LH& Hp, ad_t h, ad_t collar,
                const std::vector<ad_t>& tr, const std::vector<double>& tr0) {
  ad_t p = ad_t(Hp.profit) + ad_t(Hp.dprofit_dh) * (h - ad_t(Hp.h))
         + ad_t(Hp.dprofit_dpsi) * (collar - ad_t(Hp.collar));
  for (std::size_t k = 0; k < tr.size(); ++k)
    if (Hp.inj[k] != 0.0) p += ad_t(Hp.inj[k]) * (tr[k] - ad_t(tr0[k]));
  return p;
}

// Linearised dprofit_dpsi (the gradient-ascent rate driver) about Hp:
// dpsi0 + d2p_dpsi2*(collar - collar0) + d2p_dpsidh*(h - h0)
//       + sum_k d2p_dpsidth_k*(theta_k - theta_k0).
ad_t dpsi_lin(const LH& Hp, ad_t h, ad_t collar,
              const std::vector<ad_t>& tr, const std::vector<double>& tr0) {
  ad_t d = ad_t(Hp.dprofit_dpsi) + ad_t(Hp.d2p_dpsi2) * (collar - ad_t(Hp.collar))
         + ad_t(Hp.d2p_dpsidh) * (h - ad_t(Hp.h));
  for (std::size_t k = 0; k < tr.size(); ++k)
    if (Hp.d2p_dpsidth[k] != 0.0) d += ad_t(Hp.d2p_dpsidth[k]) * (tr[k] - ad_t(tr0[k]));
  return d;
}

// Build TF24ProdPars<ad_t> from base doubles pd, overwriting any requested trait
// that names a kernel (cascade / area / allometry) field with its active scalar.
// Leaf-only traits (vcmax_25, K_s, b, c, ...) have no kernel field -- they enter
// solely through the harvested profit channel (inj), so they are skipped here.
void set_kernel_trait(plant::TF24ProdPars<ad_t>& p, const std::string& t, const ad_t& v) {
  if      (t=="lma")  p.lma=v;  else if (t=="rho")  p.rho=v;  else if (t=="theta") p.theta=v;
  else if (t=="a_b1") p.a_b1=v; else if (t=="a_r1") p.a_r1=v; else if (t=="r_l")   p.r_l=v;
  else if (t=="r_s")  p.r_s=v;  else if (t=="r_b")  p.r_b=v;  else if (t=="r_r")   p.r_r=v;
  else if (t=="k_l")  p.k_l=v;  else if (t=="k_b")  p.k_b=v;  else if (t=="k_s")   p.k_s=v;
  else if (t=="k_r")  p.k_r=v;  else if (t=="a_bio")p.a_bio=v;else if (t=="a_y")   p.a_y=v;
  else if (t=="a_l1") p.a_l1=v; else if (t=="a_l2") p.a_l2=v; else if (t=="a_f1")  p.a_f1=v;
  else if (t=="a_f2") p.a_f2=v; else if (t=="hmat") p.hmat=v; else if (t=="omega") p.omega=v;
  else if (t=="a_f3") p.a_f3=v; else if (t=="d_I")  p.d_I=v;  else if (t=="a_dG1") p.a_dG1=v;
  else if (t=="a_dG2")p.a_dG2=v;
  // (else: leaf-only trait -> profit channel only)
}

plant::TF24ProdPars<ad_t> pf_active(const plant::TF24ProdPars<double>& pd,
    const std::vector<std::string>& traits, const std::vector<ad_t>& tr) {
  plant::TF24ProdPars<ad_t> p;
  p.lma=pd.lma;p.rho=pd.rho;p.theta=pd.theta;p.a_b1=pd.a_b1;p.a_r1=pd.a_r1;p.eta_c=pd.eta_c;
  p.r_l=pd.r_l;p.r_s=pd.r_s;p.r_b=pd.r_b;p.r_r=pd.r_r;p.k_l=pd.k_l;p.k_b=pd.k_b;p.k_s=pd.k_s;
  p.k_r=pd.k_r;p.a_bio=pd.a_bio;p.a_y=pd.a_y;p.a_l1=pd.a_l1;p.a_l2=pd.a_l2;p.a_f1=pd.a_f1;
  p.a_f2=pd.a_f2;p.hmat=pd.hmat;p.omega=pd.omega;p.a_f3=pd.a_f3;p.d_I=pd.d_I;
  p.a_dG1=pd.a_dG1;p.a_dG2=pd.a_dG2;
  for (std::size_t k = 0; k < traits.size(); ++k) set_kernel_trait(p, traits[k], tr[k]);
  return p;
}

// The 7-state cohort on the tape: 5 demog + collar + log_density.
template <typename S> struct St7ad { plant::FF16State<S> demog; S collar; S log_density; };

} // namespace


// Compiled core of tf24f_census_gradient_ad(). Returns {jacobian = metrics x
// traits, values}. Harvests every cohort's per-stage leaf operating point (fwd +
// backward-for-g') in a double pass, then replays the 7-state system onto ONE
// tape per cohort and takes one reverse sweep per metric. The census metrics
// couple cohorts through the height-trapezium, so the reduction is taped over all
// cohorts' final states together (one tape spanning the stand, sweeping per
// metric). FROZEN resident light (the rare-mutant / invasion census gradient).
// Shared core (native EH; see tf24f_census_recon_core). The R-list (`_impl`) and
// live-SCM (`_native`) entries below differ only in how EH is sourced.
Rcpp::List tf24f_census_gradient_ad_core(
    Rcpp::NumericVector pp,
    const std::vector<std::vector<plant::TF24_Environment>>& EH,
    std::vector<double> sh, std::vector<int> birth, double birth_rate,
    double k_acclim, bool use_ad_gradient, std::vector<std::string> traits,
    std::vector<std::string> metrics, double trait_rel_step) {
  using std::exp; using std::log;
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm);

  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double kI = s.pars.k_I, recruitment_decay = s.pars.recruitment_decay;
  const double GEPS = 1e-6;
  const std::size_t T = traits.size(), M = metrics.size();

  // Base trait values + per-trait FD steps; build the +/- perturbed strategies once
  // (reused across every stage harvest) and the IFT injection d(h0)/d(theta_k) and
  // d(collar0)/d(theta_k) (the latter taken at the seed operating point below).
  Rcpp::CharacterVector ppn = pp.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in strategy pars): " + traits[k]);
    tr0[k] = pp[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(pp); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(pp); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim, use_ad_gradient));
    sm.push_back(make_tf24f(q2, k_acclim, use_ad_gradient));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  const std::size_t N = EH.size(), nC = birth.size();
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];

  auto env_at = [&](std::size_t n, int stage) -> const plant::TF24_Environment& {
    return (stage == 0) ? ((n > 0) ? EH[n - 1][5] : EH[0][0]) : EH[n][stage - 1];
  };

  // --- Per-cohort harvest of the double trajectory (mirrors tf24f_census_recon_core's
  // replay, recording the per-stage operating point + the seed). ----------------
  struct CohortH {
    std::size_t b; double collar0; LH seed;            // birth seed (h0, collar0)
    std::vector<StageH> stages;                        // per RK stage
  };
  auto harvest_cohort = [&](std::size_t i) -> CohortH {
    CohortH C; const std::size_t b = (std::size_t)birth[i]; C.b = b;
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];
    // collar0 = optimum at birth (the base optimiser run once).
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    C.collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    C.seed = harvest_point(s, sp, sm, dk, eb, h0v, C.collar0);

    // Re-evolve the 7-state trajectory in double, harvesting per stage. We need only
    // the (h, collar) sequence the recon visits; reuse its deriv/axpy structure.
    auto deriv = [&](const St7& y, std::size_t n, int stage) -> St7 {
      const plant::TF24_Environment& e = env_at(n, stage);
      const double h = y.demog.height, collar = y.collar;
      C.stages.push_back(StageH{
        harvest_point(s, sp, sm, dk, e, h, collar),
        harvest_point(s, sp, sm, dk, e, h - GEPS, collar)});
      double dpsi = 0.0;
      const double net = net_at_tracked(s, e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      const double collar_dt = k_acclim * dpsi;
      double dpsi_b = 0.0;
      const double net_b = net_at_tracked(s, e, h - GEPS, collar, dpsi_b);
      const double al_b = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h - GEPS);
      const double g_back = plant::tf24_height_dt_from_net<double>(pd, h - GEPS, al_b, net_b);
      const double gprime = (r.height_dt - g_back) / GEPS;
      const double log_density_dt = -gprime - r.mortality_dt;
      return St7{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, log_density_dt};
    };
    auto axpy = [](const St7& a, double c, const St7& k) -> St7 {
      return St7{plant::FF16State<double>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar,
        a.log_density+c*k.log_density};
    };
    // Seed the same way the recon does (so the harvested trajectory matches).
    const double decay = std::exp(-recruitment_decay * sh[b]);
    double dpsi0 = 0.0;
    const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, eb, h0v, C.collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double mort0 = (pr_estab > 0.0) ? -log(pr_estab)
                                          : std::numeric_limits<double>::infinity();
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    const double logd0 = (g0 > 0.0 && pr_estab > 0.0) ? log(birth_rate * pr_estab / g0)
                                          : -std::numeric_limits<double>::infinity();
    St7 y{plant::FF16State<double>{h0v, mort0, 0.0, 0.0, 0.0}, C.collar0, logd0};
    plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    return C;
  };

  std::vector<CohortH> CH(nC);
  for (std::size_t i = 0; i < nC; ++i) CH[i] = harvest_cohort(i);

  // Pending-seed (new_node) tail: seed harvested in the FINAL-time env.
  const plant::TF24_Environment& ef = EH[N - 1][5];
  double collar0_new;
  { s.initializing_ = true; s.net_mass_production_dt(ef, h0v, s.area_leaf(h0v), 1.0 / h0v);
    collar0_new = -s.leaf.root_collar_psi_; s.initializing_ = false; }
  LH seed_new = harvest_point(s, sp, sm, dk, ef, h0v, collar0_new);
  const double decay_new = std::exp(-recruitment_decay * sh[N]);

  // --- Tape: replay the 7-state system for cohort i with active traits, returning
  // the final (height, log_density, mass_heartwood) the census reduction needs. ----
  auto replay_ad = [&](std::size_t i, const std::vector<ad_t>& tr,
                       const plant::TF24ProdPars<ad_t>& pf,
                       ad_t& H_out, ad_t& logd_out, ad_t& mhw_out) {
    const CohortH& C = CH[i];
    const std::size_t b = C.b;
    // Active h0 (IFT) and collar0 (= optimum: d(collar0)/dtheta = -d2p_dpsidth/d2p_dpsi2).
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
    ad_t collar0 = ad_t(C.collar0);
    const double inv_c2 = (C.seed.d2p_dpsi2 != 0.0) ? 1.0 / C.seed.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dcollar0 = -C.seed.d2p_dpsidth[k] * inv_c2;
      if (dcollar0 != 0.0) collar0 += ad_t(dcollar0) * (tr[k] - ad_t(tr0[k]));
    }
    // Establishment + seedling growth -> initial mortality, log_density.
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0);
    ad_t profit0 = profit_lin(C.seed, h0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0, area0, profit0);
    ad_t uu = ad_t(a_d0) * area0 / net0;
    const double decay = std::exp(-recruitment_decay * sh[b]);
    ad_t pr_estab = ad_t(decay) / (uu * uu + ad_t(1.0));
    ad_t mort0 = -log(pr_estab);
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf, h0, area0, net0);
    ad_t logd0 = log(ad_t(birth_rate) * pr_estab / g0);

    std::size_t idx = 0;
    auto deriv = [&](const St7ad<ad_t>& y, std::size_t n, int stage) -> St7ad<ad_t> {
      const StageH& S = C.stages[idx++];
      ad_t h = y.demog.height, collar = y.collar;
      ad_t profit = profit_lin(S.fwd, h, collar, tr, tr0);
      ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf, h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf, h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim) * dpsi_lin(S.fwd, h, collar, tr, tr0);
      // g' = backward FD of height_dt (collar held), the SCM's scheme: the backward
      // profit is linearised about S.back (harvested at h-GEPS), height shift -GEPS.
      ad_t profit_b = profit_lin(S.back, h - ad_t(GEPS), collar, tr, tr0);
      ad_t al_b = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h - ad_t(GEPS));
      ad_t net_b = plant::tf24_net_mass_production<ad_t>(pf, h - ad_t(GEPS), al_b, profit_b);
      ad_t g_back = plant::tf24_height_dt_from_net<ad_t>(pf, h - ad_t(GEPS), al_b, net_b);
      ad_t gprime = (r.height_dt - g_back) / ad_t(GEPS);
      ad_t log_density_dt = -gprime - r.mortality_dt;
      return St7ad<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, log_density_dt};
    };
    auto axpy = [](const St7ad<ad_t>& a, double c, const St7ad<ad_t>& k) -> St7ad<ad_t> {
      return St7ad<ad_t>{plant::FF16State<ad_t>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar,
        a.log_density+c*k.log_density};
    };
    St7ad<ad_t> y{plant::FF16State<ad_t>{h0, mort0, ad_t(0), ad_t(0), ad_t(0)}, collar0, logd0};
    St7ad<ad_t> f = plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    H_out = f.demog.height; logd_out = f.log_density; mhw_out = f.demog.mass_heartwood;
  };

  // Active pending-seed density at the final env (the trapezium tail term).
  auto dens_new_ad = [&](const std::vector<ad_t>& tr,
                         const plant::TF24ProdPars<ad_t>& pf) -> ad_t {
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
    ad_t collar0 = ad_t(collar0_new);
    const double inv_c2 = (seed_new.d2p_dpsi2 != 0.0) ? 1.0 / seed_new.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -seed_new.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0 += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    }
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0);
    ad_t profit0 = profit_lin(seed_new, h0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0, area0, profit0);
    ad_t uu = ad_t(a_d0) * area0 / net0;
    ad_t pr_estab = ad_t(decay_new) / (uu * uu + ad_t(1.0));
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf, h0, area0, net0);
    return ad_t(birth_rate) * pr_estab / g0;
  };

  // ONE tape over all cohorts (census couples them via the trapezium); one reverse
  // sweep per metric. Record the forward replay, then reduce each metric.
  ad::tape_type tape;
  std::vector<ad_t> tr(T);
  for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
  for (auto& x : tr) tape.registerInput(x);
  tape.newRecording();
  plant::TF24ProdPars<ad_t> pf = pf_active(pd, traits, tr);

  std::vector<ad_t> Hf(nC), Lf(nC), Mf(nC);
  for (std::size_t i = 0; i < nC; ++i) replay_ad(i, tr, pf, Hf[i], Lf[i], Mf[i]);
  ad_t h0a = ad_t(h0v);
  for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0a += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
  ad_t dens_new = dens_new_ad(tr, pf);

  // Descending-height order (frozen, on the double values) for the trapezium walk.
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return as_dbl(Hf[a]) > as_dbl(Hf[b]); });

  auto census_reduce = [&](auto psi) -> ad_t {
    std::vector<ad_t> phi(nC);
    for (std::size_t i = 0; i < nC; ++i) phi[i] = psi(Hf[i], exp(Lf[i]), Mf[i]);
    ad_t J = ad_t(0.0);
    for (std::size_t j = 0; j + 1 < nC; ++j) {
      const std::size_t a = ord[j], b = ord[j + 1];
      J += ad_t(0.5) * (Hf[a] - Hf[b]) * (phi[a] + phi[b]);
    }
    if (nC > 0) {
      const std::size_t last = ord[nC - 1];
      ad_t phi_new = psi(h0a, dens_new, ad_t(0.0));
      J += ad_t(0.5) * (Hf[last] - h0a) * (phi[last] + phi_new);
    }
    return J;
  };

  std::vector<ad_t> J(M);
  for (std::size_t m = 0; m < M; ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t {
        return dens * ad_t(kI) * plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h); });
    } else if (nm == "biomass") {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t mhw) -> ad_t {
        ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
        ad_t mass_leaf    = al * pf.lma;
        ad_t area_sapwood = al * pf.theta;
        ad_t mass_sapwood = area_sapwood * h * pf.eta_c * pf.rho;
        ad_t area_bark    = pf.a_b1 * al * pf.theta;
        ad_t mass_bark    = area_bark * h * pf.eta_c * pf.rho;
        ad_t mass_root    = pf.a_r1 * al;
        return dens * (mass_leaf + mass_sapwood + mass_bark + mass_root + mhw); });
    } else { // size_moment
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t { return dens * h; });
    }
  }

  // Reduce each metric: register, seed, sweep (clearing derivatives between metrics).
  Rcpp::NumericMatrix jac(M, T);
  Rcpp::NumericVector values(M);
  for (std::size_t m = 0; m < M; ++m) { values[m] = as_dbl(J[m]); tape.registerOutput(J[m]); }
  for (std::size_t m = 0; m < M; ++m) {
    tape.clearDerivatives();
    xad::derivative(J[m]) = 1.0;
    tape.computeAdjoints();
    for (std::size_t k = 0; k < T; ++k) jac(m, k) = xad::derivative(tr[k]);
  }
  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac, Rcpp::Named("values") = values);
}

// Native-SCM census gradient: faithful per-RK-stage env + birth steps from the live
// Patch -- the whole harvest is native, so the R wrapper never calls tf24f_harvest
// (whose ppsurv/tw loop is dead work for census). birth_rate<0 recovers natively.
// [[Rcpp::export]]
Rcpp::List tf24f_census_gradient_ad_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, double birth_rate,
    double k_acclim, bool use_ad_gradient, std::vector<std::string> traits,
    std::vector<std::string> metrics, double trait_rel_step) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::TF24f_Strategy, plant::TF24_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  const std::vector<int> birth = tf24f_birth_steps(patch, (std::size_t)species);
  return tf24f_census_gradient_ad_core(pp, patch.environment_history,
    patch.step_history, birth, br, k_acclim, use_ad_gradient, traits,
    metrics, trait_rel_step);
}


// ===========================================================================
// TF24f individual grow-to-size reverse-mode AD trait gradient (#472 scope B,
// build-order step 4 -- the "individuals" surface). A single TF24f plant grown
// in a FIXED environment to target size(s), differentiated w.r.t. traits: the
// AD refine of tf24f_grow_individual_to_size_gradient_fd. The TF24f analogue of
// ff16_grow_to_size_gradient_impl, but the trajectory carries the tracked collar
// (opt_root_psi_state) as a 6th state -- the SAME curvature-linearised collar
// machinery the census tape uses -- and there is NO canopy / density / g' (the
// env is given), so it is the lightest TF24f AD surface.
//
// No resident feedback, so the only machinery beyond the replay is the stopping-
// time implicit-function step. Pass 1 (R) runs the live grow (grow_individual_bracket)
// to harvest its adaptive Cash-Karp schedule + initial state. Pass 2 (here):
// discover t* (the size component crosses target) over the FROZEN schedule reading
// the FIXED env, harvest the per-RK-stage leaf operating point at the tracked
// collar along the schedule-to-t*, replay the 6-state system as a leaf-opt-free
// tapeable expression, and take one reverse sweep per state component. The collar
// starts at the individual's birth value (0, clamped into the feasible interval on
// the first eval) -- theta-independent, so unlike the census seed it needs no
// collar-IFT injection; only the seedling height h0 carries d(h0)/d(theta).
// The stopping time responds via the IFT on size(t*, theta) = target:
//   d(t*)/d(theta) = -(d size/d theta|t*) / size_dt(t*),
// and the TOTAL derivative of each component is the partial plus rate*d(t*)/d(theta).

namespace {

// 6-state grow cohort on the tape: 5 demog + tracked collar.
template <typename S> struct St6 { plant::FF16State<S> demog; S collar; };

template <typename S>
St6<S> st6_axpy(const St6<S>& a, double c, const St6<S>& k) {
  return St6<S>{plant::FF16State<S>{
    a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
    a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
    a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar};
}
double st6_at(const St6<double>& s, int c) {
  switch (c) { case 0: return s.demog.height; case 1: return s.demog.mortality;
    case 2: return s.demog.fecundity; case 3: return s.demog.area_heartwood;
    case 4: return s.demog.mass_heartwood; default: return s.collar; }
}

} // namespace


// Compiled core of tf24f_grow_individual_to_size_gradient(). Returns, per target
// size: the stopping time t*, the ODE state at t*, d(t*)/d(theta) and the TOTAL
// d(state at t*)/d(theta) (a sizes x component x trait array). y0v is the live
// individual's initial ODE state {height, mortality, fecundity, area_heartwood,
// mass_heartwood, opt_root_psi_state}; sh the harvested adaptive step-time schedule
// (grow_individual_bracket$time); sidx the size component (0 = height). shading /
// gss_tol / k_acclim / use_ad_gradient mirror the live individual's strategy.
// [[Rcpp::export]]
Rcpp::List tf24f_grow_to_size_gradient_impl(
    Rcpp::NumericVector pp, plant::TF24_Environment env, Rcpp::NumericVector y0v,
    std::vector<double> sh, std::vector<double> targets, int sidx,
    std::vector<std::string> traits, double k_acclim, bool use_ad_gradient,
    std::string shading, double gss_tol, double trait_rel_step) {
  using std::exp; using std::log;
  if (sh.size() < 2) Rcpp::stop("schedule must have at least two step times");
  if (y0v.size() < 6) Rcpp::stop("y0 must have the 6 TF24f ODE states");

  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient, shading, gss_tol);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const std::size_t T = traits.size();

  // Base trait values + per-trait FD steps; build the +/- perturbed strategies once
  // (reused across targets/stages) and the seedling-size IFT d(h0)/d(theta).
  Rcpp::CharacterVector ppn = pp.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  const double h0v = y0v[0];
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in strategy pars): " + traits[k]);
    tr0[k] = pp[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(pp); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(pp); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim, use_ad_gradient, shading, gss_tol));
    sm.push_back(make_tf24f(q2, k_acclim, use_ad_gradient, shading, gss_tol));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  std::vector<double> step_h(sh.size() - 1);
  for (std::size_t n = 0; n + 1 < sh.size(); ++n) step_h[n] = sh[n + 1] - sh[n];

  // Double 6-state derivative reading the fixed env (real leaf at tracked collar).
  auto deriv_d = [&](const St6<double>& y, std::size_t, int) -> St6<double> {
    double dpsi = 0.0;
    const double net = net_at_tracked(s, env, y.demog.height, y.collar, dpsi);
    const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, y.demog.height);
    plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, y.demog.height, al, net, true);
    return St6<double>{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
      r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim * dpsi};
  };
  auto axpy_d = [](const St6<double>& a, double c, const St6<double>& k){ return st6_axpy<double>(a,c,k); };
  auto replay_one_d = [&](St6<double> y, const std::vector<double>& sched) -> St6<double> {
    return plant::ff16_cashkarp_replay(y, sched, (std::size_t)0, deriv_d, axpy_d);
  };

  St6<double> y0{plant::FF16State<double>{y0v[0], y0v[1], y0v[2], y0v[3], y0v[4]}, y0v[5]};

  const std::vector<std::string> comp =
    {"height","mortality","fecundity","area_heartwood","mass_heartwood","opt_root_psi_state"};
  const std::size_t nS = comp.size(), nG = targets.size();
  Rcpp::NumericVector tstar(nG);
  Rcpp::NumericMatrix state(nG, nS);
  Rcpp::NumericMatrix dtime(nG, T);
  Rcpp::NumericVector dstate(nG * nS * T);
  auto DS = [&](std::size_t g, std::size_t c, std::size_t k) -> double& {
    return dstate[g + nG * (c + nS * k)]; };

  for (std::size_t g = 0; g < nG; ++g) {
    // --- Discover t*: walk the schedule one RKCK step at a time until the size
    // component crosses target, then bisect the partial final step. -------------
    int nfull = -1; double dt_final = 0.0; St6<double> fin{};
    {
      St6<double> y = y0;
      for (std::size_t n = 0; n < step_h.size(); ++n) {
        St6<double> ynext = replay_one_d(y, {step_h[n]});
        if (st6_at(ynext, sidx) >= targets[g]) {
          nfull = (int)n;
          double lo = 0.0, hi = step_h[n];
          for (int it = 0; it < 80; ++it) {
            double mid = 0.5 * (lo + hi);
            double sm_ = st6_at(replay_one_d(y, {mid}), sidx);
            if (sm_ < targets[g]) lo = mid; else hi = mid;
          }
          dt_final = 0.5 * (lo + hi);
          fin = replay_one_d(y, {dt_final});
          break;
        }
        y = ynext;
      }
    }
    if (nfull < 0) Rcpp::stop("target size not reached within the schedule");
    tstar[g] = sh[nfull] + dt_final;
    for (std::size_t c = 0; c < nS; ++c) state(g, c) = st6_at(fin, (int)c);

    // Rates at t* (double) for the IFT correction.
    St6<double> rate = deriv_d(fin, 0, 0);
    double yd[6] = {rate.demog.height, rate.demog.mortality, rate.demog.fecundity,
                    rate.demog.area_heartwood, rate.demog.mass_heartwood, rate.collar};
    const double sdot = yd[sidx];

    // Frozen schedule to t*: nfull full steps + the single partial step dt_final.
    std::vector<double> sched(step_h.begin(), step_h.begin() + nfull);
    sched.push_back(dt_final);

    // --- Harvest the per-RK-stage leaf operating point along sched (double). ----
    std::vector<LH> stages;
    {
      auto deriv_h = [&](const St6<double>& y, std::size_t, int) -> St6<double> {
        stages.push_back(harvest_point(s, sp, sm, dk, env, y.demog.height, y.collar));
        return deriv_d(y, 0, 0);
      };
      St6<double> y = y0;
      plant::ff16_cashkarp_replay(y, sched, (std::size_t)0, deriv_h, axpy_d);
    }

    // --- Tape replay over sched: partial d(state)/d(theta) holding t* fixed. -----
    ad::tape_type tape;
    std::vector<ad_t> tr(T);
    for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
    for (auto& x : tr) tape.registerInput(x);
    tape.newRecording();
    plant::TF24ProdPars<ad_t> pf = pf_active(pd, traits, tr);
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));

    std::size_t idx = 0;
    auto deriv_ad = [&](const St6<ad_t>& y, std::size_t, int) -> St6<ad_t> {
      const LH& H = stages[idx++];
      ad_t h = y.demog.height, collar = y.collar;
      ad_t profit = profit_lin(H, h, collar, tr, tr0);
      ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf, h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf, h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim) * dpsi_lin(H, h, collar, tr, tr0);
      return St6<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt};
    };
    auto axpy_ad = [](const St6<ad_t>& a, double c, const St6<ad_t>& k){ return st6_axpy<ad_t>(a,c,k); };
    St6<ad_t> yad{plant::FF16State<ad_t>{h0, ad_t(y0v[1]), ad_t(y0v[2]), ad_t(y0v[3]),
      ad_t(y0v[4])}, ad_t(y0v[5])};
    St6<ad_t> out = plant::ff16_cashkarp_replay(yad, sched, (std::size_t)0, deriv_ad, axpy_ad);
    ad_t oc[6] = {out.demog.height, out.demog.mortality, out.demog.fecundity,
                  out.demog.area_heartwood, out.demog.mass_heartwood, out.collar};
    for (std::size_t c = 0; c < nS; ++c) tape.registerOutput(oc[c]);
    std::vector<std::vector<double>> P(nS, std::vector<double>(T, 0.0));
    for (std::size_t c = 0; c < nS; ++c) {
      tape.clearDerivatives();
      xad::derivative(oc[c]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < T; ++k) P[c][k] = xad::derivative(tr[k]);
    }
    // IFT: d(t*)/d(theta_k) = -(d size/d theta_k|t*) / size_dt(t*); total = P + ydot*dt*.
    for (std::size_t k = 0; k < T; ++k) {
      double dtk = (sdot != 0.0) ? -P[sidx][k] / sdot : 0.0;
      dtime(g, k) = dtk;
      for (std::size_t c = 0; c < nS; ++c) DS(g, c, k) = P[c][k] + yd[c] * dtk;
    }
  }

  state.attr("dimnames")  = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp));
  dtime.attr("dimnames")  = Rcpp::List::create(R_NilValue, Rcpp::wrap(traits));
  dstate.attr("dim")      = Rcpp::IntegerVector::create((int)nG, (int)nS, (int)T);
  dstate.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp), Rcpp::wrap(traits));
  return Rcpp::List::create(Rcpp::Named("time") = tstar, Rcpp::Named("state") = state,
                            Rcpp::Named("d_time") = dtime, Rcpp::Named("d_state") = dstate);
}


// ===========================================================================
// TF24f RESIDENT (coupled) census reverse-mode AD trait gradient (#472 scope B,
// build-order step 5 -- the resident-feedback surface). The TOTAL stand-level
// d(census metric)/d(theta): all cohorts are re-evolved TOGETHER over the frozen
// schedule and the canopy light each RK stage is reconstructed from the ACTIVE
// stand (heights AND densities respond to theta), so EVERY trait re-shades the
// canopy that every cohort reads -- the genuine resident feedback that routinely
// flips the sign of d(census)/d(theta) relative to the frozen (rare-mutant)
// gradient. Validates against tf24f_resident_census_gradient_fd (full SCM re-runs).
//
// Mirrors the FF16 coupled engine (assemble_metrics_coupled in ff16_emergent.cpp):
// the canopy reconstruction is structurally identical (the Yokozawa height-trapezium
// of density*kI*area_leaf*Q, Q = (1-(z/h)^eta)^2 -- TF24's [eqn 10] == FF16's form),
// with TWO differences:
//   * the cohort carries the 7-state {5 demog, tracked collar, log_density}; and
//   * the per-cohort RATE is the harvested+linearised TF24f leaf at the tracked
//     collar (FF16's closed-form deep_net has no TF24f analogue), the SAME curvature
//     harvest the frozen census tape uses, PLUS a resident LIGHT channel.
//
// TF24f crown-centre reads light at the SINGLE height z = h*eta_c (not a crown
// integral), so the resident feedback needs only one canopy read L_i = canopy(z_i).
// The light channel is harvested as dprofit_dL (and the collar-rate curvature
// d2profit/dpsi dL) by a finite difference over a FLAT environment
// (set_fixed_environment) at the operating point -- the clean way to perturb the
// crown light the leaf sees while holding every other input fixed. On the tape the
// resident light correction is ANCHORED to its own value:
//   profit += dprofit_dL * (canopy_active(zv) - value(canopy_active(zv)))
// at zv = the operating-point crown height (a tape CONSTANT). The value part is zero
// (clean baseline at every theta) and the derivative is purely the resident
// reshaping; the focal plant's OWN height->light slope stays inside the frozen-census
// dprofit_dh (harvested in the frozen env), so there is no double count. This is the
// FF16 "anchored resident light" trick (resident_light_anchored) adapted.
// ===========================================================================

namespace {

double dval(double v)        { return v; }
double dval(const ad_t& v)   { return xad::value(v); }

// Competition at a frozen knot z over the active stand: see
// plant::gradient::canopy_comp_at (shared Yokozawa trapezium == FF16's; TF24's
// [eqn 10] is identical). Used by the single- and multi-species coupled canopies.

// Build the active canopy light interpolator at the frozen knot x-positions kx from
// the current active stand (alive cohorts) + the frozen boundary node (nn_h, nn_c).
// light(z) = exp(-competition(z)/area). Active in heights / densities / allometry.
template <typename S>
odelia::interpolator::basic_interpolator<S>
build_canopy(const std::vector<St7ad<S>>& stand, const std::vector<char>& alive,
             const std::vector<double>& kx, double nn_h, double nn_c,
             const plant::TF24ProdPars<S>& pf, double kI, double eta, double area,
             std::vector<double>* ly_out = nullptr) {
  using std::exp;
  std::vector<S> hv, gv;
  hv.reserve(stand.size() + 1); gv.reserve(stand.size() + 1);
  for (std::size_t i = 0; i < stand.size(); ++i) if (alive[i]) {
    S hi = stand[i].demog.height;
    S dens = exp(stand[i].log_density);
    S al = plant::tf24_area_leaf<S>(pf.a_l1, pf.a_l2, hi);
    hv.push_back(hi); gv.push_back(dens * S(kI) * al);
  }
  hv.push_back(S(nn_h)); gv.push_back(S(nn_c));
  std::vector<std::size_t> ord(hv.size());
  for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return dval(hv[a]) > dval(hv[b]); });
  std::vector<S> hs(hv.size()), gs(hv.size());
  for (std::size_t k = 0; k < ord.size(); ++k) { hs[k] = hv[ord[k]]; gs[k] = gv[ord[k]]; }
  std::vector<S> ly(kx.size());
  for (std::size_t k = 0; k < kx.size(); ++k)
    ly[k] = exp(-plant::gradient::canopy_comp_at<S>(kx[k], hs, gs, eta) * S(1.0 / area));
  if (ly_out) { ly_out->resize(kx.size());
    for (std::size_t k = 0; k < kx.size(); ++k) (*ly_out)[k] = dval(ly[k]); }
  odelia::interpolator::basic_interpolator<S> interp; interp.init(kx, ly);
  return interp;
}

// Evaluate the real TF24f leaf at a GIVEN crown light L (canopy openness 0-1),
// holding everything else from the per-stage env e_tmpl: a FLAT light spline at L
// makes crown-centre's single get_environment_at_height(h*eta_c) read return L while
// the soil / vpd / ca / PPFD state is preserved. Returns the leaf profit, the
// acclimation rate driver dprofit_dpsi, and (out) the kernel net at the operating
// point. The clean primitive for the light-channel FD harvest and the R0 re-evolution.
void leaf_at_light(plant::TF24f_Strategy& s, const plant::TF24_Environment& e_tmpl,
                   double L, double h, double collar,
                   double& profit, double& dpsi) {
  plant::TF24_Environment ef = e_tmpl;
  ef.set_fixed_environment(L, 1.0e4);
  s.tracked_root_psi_ = collar;
  s.net_mass_production_dt(ef, h, s.area_leaf(h), 1.0 / h);
  profit = s.leaf.profit_;
  dpsi   = s.dprofit_dpsi_;
}

// The frozen-census operating-point harvest (LH) extended with the resident light
// channel: dprofit_dL and the collar-rate light curvature d2profit/dpsi dL, both by
// a central FD over the crown light at the operating point (flat-env evaluator).
struct LHc { LH base; double dprofit_dL; double d2p_dpsidL; };

LHc harvest_point_c(plant::TF24f_Strategy& s,
                    std::vector<plant::TF24f_Strategy>& sp,
                    std::vector<plant::TF24f_Strategy>& sm,
                    const std::vector<double>& dk,
                    const plant::TF24_Environment& e, double h, double collar,
                    double eta_c) {
  LHc H; H.base = harvest_point(s, sp, sm, dk, e, h, collar);
  const double L0 = e.get_environment_at_height(h * eta_c);
  const double dL = 1e-5 * std::max(std::abs(L0), 1e-4);
  double pP, dpsiP, pM, dpsiM;
  leaf_at_light(s, e, L0 + dL, h, collar, pP, dpsiP);
  leaf_at_light(s, e, L0 - dL, h, collar, pM, dpsiM);
  H.dprofit_dL = (pP - pM) / (2.0 * dL);
  H.d2p_dpsidL = (dpsiP - dpsiM) / (2.0 * dL);
  return H;
}

// Per-stage harvest: forward operating point (h, collar) + the backward point
// (h - GEPS, collar) the census g' (backward FD of height_dt) reads.
struct StageHc { LHc fwd, back; };

// Anchored resident light correction at the active crown light Lact: value 0 (clean
// baseline at every theta), derivative = the resident reshaping d(canopy)/d(theta).
ad_t anchor(ad_t Lact) { return Lact - ad_t(xad::value(Lact)); }

// Single Cash-Karp RKCK step over the whole-stand state (peeled from ff16's driver
// so cohorts can be introduced between steps).
template <typename State, typename DerivFn, typename AxpyFn>
State rkck_one_step_tf(State y, double h, std::size_t rn, DerivFn&& deriv, AxpyFn&& axpy) {
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

template <typename S>
std::vector<St7ad<S>> stand_axpy(const std::vector<St7ad<S>>& a, double c,
                                 const std::vector<St7ad<S>>& k) {
  std::vector<St7ad<S>> r(a.size());
  for (std::size_t i = 0; i < a.size(); ++i)
    r[i] = St7ad<S>{plant::FF16State<S>{
      a[i].demog.height+c*k[i].demog.height, a[i].demog.mortality+c*k[i].demog.mortality,
      a[i].demog.fecundity+c*k[i].demog.fecundity, a[i].demog.area_heartwood+c*k[i].demog.area_heartwood,
      a[i].demog.mass_heartwood+c*k[i].demog.mass_heartwood}, a[i].collar+c*k[i].collar,
      a[i].log_density+c*k[i].log_density};
  return r;
}

} // namespace


// R0 gate for the TF24f COUPLED replay (#472 scope B step 5): a double-precision
// whole-stand re-evolution reconstructing the ACTIVE canopy each RK stage and driving
// the real TF24f leaf at the tracked collar at the reconstructed crown light. Returns
// the reconstructed TOTAL-stand census metric VALUES and `env_err` = worst
// |reconstructed knot light - SCM knot light| over all stages (the coupled-drift
// gauge), confirming the active re-evolution reproduces the resident stand before the
// AD tape. Seeds each cohort from the FROZEN birth env (the active-birthenv channel is
// deferred, as in FF16); the resident feedback enters through the growth-phase canopy.
// [[Rcpp::export]]
Rcpp::List tf24f_coupled_metrics_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, double birth_rate, double k_acclim, bool use_ad_gradient,
    std::vector<std::string> metrics, Rcpp::List nn_h_list, Rcpp::List nn_c_list,
    double patch_area) {
  using std::exp; using std::log;
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm);
  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double kI = s.pars.k_I, eta = s.pars.eta, eta_c = s.eta_c;
  const double recruitment_decay = s.pars.recruitment_decay, GEPS = 1e-6;
  const std::size_t N = eh_list.size(), nC = birth.size();

  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));
  }
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];
  std::vector<std::vector<double>> NNH(N), NNC(N);
  for (std::size_t n = 0; n < N; ++n) {
    NNH[n] = Rcpp::as<std::vector<double>>(nn_h_list[n]);
    NNC[n] = Rcpp::as<std::vector<double>>(nn_c_list[n]);
  }
  auto env_idx = [&](std::size_t rn, int rs, std::size_t& en, int& es) {
    if (rs == 0) { if (rn > 0) { en = rn - 1; es = 5; } else { en = 0; es = 0; } }
    else { en = rn; es = rs - 1; }
  };

  std::vector<St7ad<double>> stand(nC);
  std::vector<char> alive(nC, 0);
  double env_err = 0.0;

  auto deriv = [&](const std::vector<St7ad<double>>& y, std::size_t rn, int rs)
      -> std::vector<St7ad<double>> {
    std::size_t en; int es; env_idx(rn, rs, en, es);
    const std::vector<double>& kx = EH[en][es].light_availability.spline.get_x();
    std::vector<double> ly0;
    auto interp = build_canopy<double>(y, alive, kx, NNH[en][es], NNC[en][es],
                                       pd, kI, eta, patch_area, &ly0);
    const std::vector<double>& y0 = EH[en][es].light_availability.spline.get_y();
    for (std::size_t k = 0; k < ly0.size() && k < y0.size(); ++k)
      env_err = std::max(env_err, std::abs(ly0[k] - y0[k]));
    const double cap = kx.back();
    std::vector<St7ad<double>> dy(nC);
    for (std::size_t i = 0; i < nC; ++i) {
      if (!alive[i]) { dy[i] = St7ad<double>{plant::FF16State<double>{0,0,0,0,0}, 0, 0}; continue; }
      const double h = y[i].demog.height, collar = y[i].collar;
      const double z = h * eta_c;
      const double L = (z > cap) ? 1.0 : interp(z);
      double profit, dpsi;
      leaf_at_light(s, EH[en][es], L, h, collar, profit, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      const double net = plant::tf24_net_mass_production<double>(pd, h, al, profit);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      // g' = backward FD of height_dt (collar held), reading the canopy at h-GEPS.
      const double zb = (h - GEPS) * eta_c;
      const double Lb = (zb > cap) ? 1.0 : interp(zb);
      double pb, dpb; leaf_at_light(s, EH[en][es], Lb, h - GEPS, collar, pb, dpb);
      const double alb = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h - GEPS);
      const double netb = plant::tf24_net_mass_production<double>(pd, h - GEPS, alb, pb);
      const double gback = plant::tf24_height_dt_from_net<double>(pd, h - GEPS, alb, netb);
      const double gprime = (r.height_dt - gback) / GEPS;
      dy[i] = St7ad<double>{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim * dpsi, -gprime - r.mortality_dt};
    }
    return dy;
  };
  auto axpy = [](const std::vector<St7ad<double>>& a, double c,
                 const std::vector<St7ad<double>>& k){ return stand_axpy<double>(a, c, k); };

  // Establish cohort i in its frozen birth env at step index `bstep` (clamped to the
  // last step for a final-boundary cohort, birth == N).
  auto establish = [&](std::size_t i, std::size_t bstep) {
    const std::size_t en = (bstep > 0) ? std::min(bstep, N) - 1 : 0;
    const plant::TF24_Environment& eb = (bstep > 0) ? EH[en][5] : EH[0][0];
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    const double collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    const double decay = std::exp(-recruitment_decay * sh[bstep]);
    double dpsi0; const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, eb, h0v, collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double mort0 = (pr_estab > 0.0) ? -log(pr_estab)
                                          : std::numeric_limits<double>::infinity();
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    const double logd0 = (g0 > 0.0 && pr_estab > 0.0) ? log(birth_rate * pr_estab / g0)
                                          : -std::numeric_limits<double>::infinity();
    stand[i] = St7ad<double>{plant::FF16State<double>{h0v, mort0, 0, 0, 0}, collar0, logd0};
    alive[i] = 1;
  };
  // March: introduce cohorts born at rn (frozen birth env seed), then step the stand.
  for (std::size_t rn = 0; rn < N; ++rn) {
    for (std::size_t i = 0; i < nC; ++i) if ((std::size_t)birth[i] == rn) establish(i, rn);
    stand = rkck_one_step_tf(stand, step_h[rn], rn, deriv, axpy);
  }
  // Final-boundary cohorts (birth >= N): established at h0, never stepped (see the AD
  // core note) -- else the zero-initialised height feeds tf24_area_leaf as 0*log(0)=NaN.
  for (std::size_t i = 0; i < nC; ++i) if (!alive[i]) establish(i, (std::size_t)birth[i]);

  // Pending-seed density at the final env (frozen tail), then census reduction.
  double dens_new = 0.0;
  {
    const plant::TF24_Environment& ef = EH[N - 1][5];
    s.initializing_ = true; s.net_mass_production_dt(ef, h0v, s.area_leaf(h0v), 1.0 / h0v);
    const double collar0 = -s.leaf.root_collar_psi_; s.initializing_ = false;
    const double decay = std::exp(-recruitment_decay * sh[N]);
    double dpsi0; const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, ef, h0v, collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    dens_new = (g0 > 0.0 && pr_estab > 0.0) ? birth_rate * pr_estab / g0 : 0.0;
  }
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return stand[a].demog.height > stand[b].demog.height; });
  auto census_reduce = [&](auto psi) -> double {
    std::vector<double> phi(nC);
    for (std::size_t i = 0; i < nC; ++i)
      phi[i] = psi(stand[i].demog.height, std::exp(stand[i].log_density), stand[i].demog.mass_heartwood);
    double J = 0.0;
    for (std::size_t j = 0; j + 1 < nC; ++j) {
      const std::size_t a = ord[j], b = ord[j + 1];
      J += 0.5 * (stand[a].demog.height - stand[b].demog.height) * (phi[a] + phi[b]);
    }
    if (nC > 0) {
      const std::size_t last = ord[nC - 1];
      const double phi_new = psi(h0v, dens_new, 0.0);
      J += 0.5 * (stand[last].demog.height - h0v) * (phi[last] + phi_new);
    }
    return J;
  };
  Rcpp::NumericVector values(metrics.size());
  for (std::size_t m = 0; m < metrics.size(); ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      values[m] = census_reduce([&](double h, double dens, double){
        return dens * kI * plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h); });
    } else if (nm == "biomass") {
      values[m] = census_reduce([&](double h, double dens, double mhw){
        const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
        const double mass_leaf = al * pd.lma, area_sapwood = al * pd.theta;
        const double mass_sapwood = area_sapwood * h * pd.eta_c * pd.rho;
        const double area_bark = pd.a_b1 * al * pd.theta;
        const double mass_bark = area_bark * h * pd.eta_c * pd.rho;
        const double mass_root = pd.a_r1 * al;
        return dens * (mass_leaf + mass_sapwood + mass_bark + mass_root + mhw); });
    } else {
      values[m] = census_reduce([&](double h, double dens, double){ return dens * h; });
    }
  }
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("values") = values,
                            Rcpp::Named("env_err") = env_err);
}


// R1 of the TF24f COUPLED replay (#472 scope B step 5): the resident TOTAL trait
// gradient of the emergent census metrics, one reverse sweep per metric over the
// coupled whole-stand re-evolution. A double DISCOVERY pass harvests every cohort's
// per-stage leaf operating point (the frozen-census curvature harvest + the resident
// LIGHT channel) along the frozen-env trajectory; the AD pass re-evolves all cohorts
// together on ONE tape, reconstructing the active canopy each RK stage so every trait
// re-shades the light every cohort reads. Validate vs tf24f_resident_census_gradient_fd.
// [[Rcpp::export]]
Rcpp::List tf24f_coupled_gradient_core(
    Rcpp::NumericVector pp,
    const std::vector<std::vector<plant::TF24_Environment>>& EH,
    std::vector<double> sh, std::vector<int> birth, double birth_rate,
    double k_acclim, bool use_ad_gradient, std::vector<std::string> traits,
    std::vector<std::string> metrics,
    const std::vector<std::vector<double>>& NNH,
    const std::vector<std::vector<double>>& NNC,
    double patch_area, double trait_rel_step) {
  using std::exp; using std::log;
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm);
  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double kI = s.pars.k_I, eta = s.pars.eta, eta_c = s.eta_c;
  const double recruitment_decay = s.pars.recruitment_decay, GEPS = 1e-6;
  const std::size_t T = traits.size(), M = metrics.size();

  // Base trait values + FD steps; prebuilt +/- perturbed strategies; IFT d(h0)/dtheta.
  Rcpp::CharacterVector ppn = pp.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in strategy pars): " + traits[k]);
    tr0[k] = pp[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(pp); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(pp); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim, use_ad_gradient));
    sm.push_back(make_tf24f(q2, k_acclim, use_ad_gradient));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  const std::size_t N = EH.size(), nC = birth.size();
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];
  auto env_idx = [&](std::size_t rn, int rs, std::size_t& en, int& es) {
    if (rs == 0) { if (rn > 0) { en = rn - 1; es = 5; } else { en = 0; es = 0; } }
    else { en = rn; es = rs - 1; }
  };
  auto env_at = [&](std::size_t rn, int rs) -> const plant::TF24_Environment& {
    std::size_t en; int es; env_idx(rn, rs, en, es); return EH[en][es];
  };

  // --- DISCOVERY: per-cohort frozen-env replay, harvesting the operating point
  // (curvature + light channel) per stage + the birth seed (mirrors the frozen
  // census tape's harvest_cohort, with harvest_point_c and the StageHc store). ----
  struct CohortH { std::size_t b; double collar0; LHc seed; std::vector<StageHc> stages; };
  auto harvest_cohort = [&](std::size_t i) -> CohortH {
    CohortH C; const std::size_t b = (std::size_t)birth[i]; C.b = b;
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    C.collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    C.seed = harvest_point_c(s, sp, sm, dk, eb, h0v, C.collar0, eta_c);
    auto deriv = [&](const St7& y, std::size_t n, int stage) -> St7 {
      const plant::TF24_Environment& e = env_at(n, stage);
      const double h = y.demog.height, collar = y.collar;
      C.stages.push_back(StageHc{
        harvest_point_c(s, sp, sm, dk, e, h, collar, eta_c),
        harvest_point_c(s, sp, sm, dk, e, h - GEPS, collar, eta_c)});
      double dpsi = 0.0;
      const double net = net_at_tracked(s, e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      double dpsi_b = 0.0;
      const double net_b = net_at_tracked(s, e, h - GEPS, collar, dpsi_b);
      const double al_b = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h - GEPS);
      const double g_back = plant::tf24_height_dt_from_net<double>(pd, h - GEPS, al_b, net_b);
      const double gprime = (r.height_dt - g_back) / GEPS;
      return St7{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim * dpsi, -gprime - r.mortality_dt};
    };
    auto axpy = [](const St7& a, double c, const St7& k) -> St7 {
      return St7{plant::FF16State<double>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar,
        a.log_density+c*k.log_density};
    };
    const double decay = std::exp(-recruitment_decay * sh[b]);
    double dpsi0; const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, eb, h0v, C.collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double mort0 = (pr_estab > 0.0) ? -log(pr_estab)
                                          : std::numeric_limits<double>::infinity();
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    const double logd0 = (g0 > 0.0 && pr_estab > 0.0) ? log(birth_rate * pr_estab / g0)
                                          : -std::numeric_limits<double>::infinity();
    St7 y{plant::FF16State<double>{h0v, mort0, 0, 0, 0}, C.collar0, logd0};
    plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    return C;
  };
  std::vector<CohortH> CH(nC);
  for (std::size_t i = 0; i < nC; ++i) CH[i] = harvest_cohort(i);

  // Pending-seed (new_node) tail: seed harvested in the FINAL-time env (no light).
  const plant::TF24_Environment& ef = EH[N - 1][5];
  double collar0_new;
  { s.initializing_ = true; s.net_mass_production_dt(ef, h0v, s.area_leaf(h0v), 1.0 / h0v);
    collar0_new = -s.leaf.root_collar_psi_; s.initializing_ = false; }
  LHc seed_new = harvest_point_c(s, sp, sm, dk, ef, h0v, collar0_new, eta_c);
  const double decay_new = std::exp(-recruitment_decay * sh[N]);

  // --- AD tape: whole-stand coupled re-evolution with the harvested linearised rate
  // + the anchored resident light channel; one reverse sweep per metric. -----------
  ad::tape_type tape;
  std::vector<ad_t> tr(T);
  for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
  for (auto& x : tr) tape.registerInput(x);
  tape.newRecording();
  plant::TF24ProdPars<ad_t> pf = pf_active(pd, traits, tr);

  std::vector<St7ad<ad_t>> stand(nC);
  std::vector<char> alive(nC, 0);
  std::vector<std::size_t> sidx(nC, 0);          // per-cohort stage-visit counter

  auto seed_h0 = [&](void) -> ad_t {
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
    return h0;
  };

  auto deriv = [&](const std::vector<St7ad<ad_t>>& y, std::size_t rn, int rs)
      -> std::vector<St7ad<ad_t>> {
    std::size_t en; int es; env_idx(rn, rs, en, es);
    const std::vector<double>& kx = EH[en][es].light_availability.spline.get_x();
    auto interp = build_canopy<ad_t>(y, alive, kx, NNH[en][es], NNC[en][es],
                                     pf, kI, eta, patch_area);
    const double cap = kx.back();
    std::vector<St7ad<ad_t>> dy(nC);
    for (std::size_t i = 0; i < nC; ++i) {
      if (!alive[i]) { dy[i] = St7ad<ad_t>{plant::FF16State<ad_t>{ad_t(0),ad_t(0),ad_t(0),ad_t(0),ad_t(0)},
        ad_t(0), ad_t(0)}; continue; }
      const StageHc& S = CH[i].stages[sidx[i]++];
      ad_t h = y[i].demog.height, collar = y[i].collar;
      // forward resident light correction (anchored: value 0, derivative = reshaping).
      ad_t z = h * ad_t(eta_c); double zv = xad::value(z);
      ad_t dLr = (zv > cap) ? ad_t(0.0) : anchor(interp(zv));
      ad_t profit = profit_lin(S.fwd.base, h, collar, tr, tr0) + ad_t(S.fwd.dprofit_dL) * dLr;
      ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf, h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf, h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim) * (dpsi_lin(S.fwd.base, h, collar, tr, tr0)
                                         + ad_t(S.fwd.d2p_dpsidL) * dLr);
      // g' = backward FD of height_dt (collar held), canopy read at h-GEPS.
      ad_t zb = (h - ad_t(GEPS)) * ad_t(eta_c); double zbv = xad::value(zb);
      ad_t dLrb = (zbv > cap) ? ad_t(0.0) : anchor(interp(zbv));
      ad_t profit_b = profit_lin(S.back.base, h - ad_t(GEPS), collar, tr, tr0)
                      + ad_t(S.back.dprofit_dL) * dLrb;
      ad_t al_b = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h - ad_t(GEPS));
      ad_t net_b = plant::tf24_net_mass_production<ad_t>(pf, h - ad_t(GEPS), al_b, profit_b);
      ad_t g_back = plant::tf24_height_dt_from_net<ad_t>(pf, h - ad_t(GEPS), al_b, net_b);
      ad_t gprime = (r.height_dt - g_back) / ad_t(GEPS);
      dy[i] = St7ad<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, -gprime - r.mortality_dt};
    }
    return dy;
  };
  auto axpy = [](const std::vector<St7ad<ad_t>>& a, double c,
                 const std::vector<St7ad<ad_t>>& k){ return stand_axpy<ad_t>(a, c, k); };

  // Establish cohort i at its birth env (frozen, baked into seed.base). `bstep` is the
  // step time index for the establishment decay; for a final-boundary cohort (CH[i].b
  // == N) it is N, whose sh[] entry is the final time.
  auto establish = [&](std::size_t i, std::size_t bstep) {
    const CohortH& C = CH[i];
    ad_t h0 = seed_h0();
    ad_t collar0 = ad_t(C.collar0);
    const double inv_c2 = (C.seed.base.d2p_dpsi2 != 0.0) ? 1.0 / C.seed.base.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -C.seed.base.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0 += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    }
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0);
    ad_t profit0 = profit_lin(C.seed.base, h0, collar0, tr, tr0);   // frozen birth env (no light)
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0, area0, profit0);
    ad_t uu = ad_t(a_d0) * area0 / net0;
    const double decay = std::exp(-recruitment_decay * sh[bstep]);
    ad_t pr_estab = ad_t(decay) / (uu * uu + ad_t(1.0));
    ad_t mort0 = -log(pr_estab);
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf, h0, area0, net0);
    ad_t logd0 = log(ad_t(birth_rate) * pr_estab / g0);
    stand[i] = St7ad<ad_t>{plant::FF16State<ad_t>{h0, mort0, ad_t(0), ad_t(0), ad_t(0)}, collar0, logd0};
    alive[i] = 1;
  };
  for (std::size_t rn = 0; rn < N; ++rn) {
    for (std::size_t i = 0; i < nC; ++i) if (CH[i].b == rn) establish(i, rn);
    stand = rkck_one_step_tf(stand, step_h[rn], rn, deriv, axpy);
  }
  // Final-boundary cohorts (CH[i].b >= N): introduced at the last node time but never
  // stepped (the discovery replay runs zero steps for them), so they sit at h0. Without
  // this they keep the zero-initialised height and tf24_area_leaf's a_l2 channel hits
  // 0*log(0)=NaN, which then poisons the whole shared census tape. Mirrors FF16's
  // assemble_metrics_coupled final-boundary establishment.
  for (std::size_t i = 0; i < nC; ++i) if (!alive[i]) establish(i, CH[i].b);

  // Pending-seed density at the final env (frozen tail; collar0 IFT injected).
  ad_t h0a = seed_h0();
  ad_t collar0n = ad_t(collar0_new);
  { const double inv_c2 = (seed_new.base.d2p_dpsi2 != 0.0) ? 1.0 / seed_new.base.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -seed_new.base.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0n += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    } }
  ad_t area0n = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0a);
  ad_t profit0n = profit_lin(seed_new.base, h0a, collar0n, tr, tr0);
  ad_t net0n = plant::tf24_net_mass_production<ad_t>(pf, h0a, area0n, profit0n);
  ad_t uun = ad_t(a_d0) * area0n / net0n;
  ad_t pr_estabn = ad_t(decay_new) / (uun * uun + ad_t(1.0));
  ad_t g0n = plant::tf24_height_dt_from_net<ad_t>(pf, h0a, area0n, net0n);
  ad_t dens_new = ad_t(birth_rate) * pr_estabn / g0n;

  // Census reduction (descending-height trapezium + pending-seed tail), one tape.
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return as_dbl(stand[a].demog.height) > as_dbl(stand[b].demog.height); });
  auto census_reduce = [&](auto psi) -> ad_t {
    std::vector<ad_t> phi(nC);
    for (std::size_t i = 0; i < nC; ++i)
      phi[i] = psi(stand[i].demog.height, exp(stand[i].log_density), stand[i].demog.mass_heartwood);
    ad_t J = ad_t(0.0);
    for (std::size_t j = 0; j + 1 < nC; ++j) {
      const std::size_t a = ord[j], b = ord[j + 1];
      J += ad_t(0.5) * (stand[a].demog.height - stand[b].demog.height) * (phi[a] + phi[b]);
    }
    if (nC > 0) {
      const std::size_t last = ord[nC - 1];
      ad_t phi_new = psi(h0a, dens_new, ad_t(0.0));
      J += ad_t(0.5) * (stand[last].demog.height - h0a) * (phi[last] + phi_new);
    }
    return J;
  };
  std::vector<ad_t> J(M);
  for (std::size_t m = 0; m < M; ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t {
        return dens * ad_t(kI) * plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h); });
    } else if (nm == "biomass") {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t mhw) -> ad_t {
        ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
        ad_t mass_leaf = al * pf.lma, area_sapwood = al * pf.theta;
        ad_t mass_sapwood = area_sapwood * h * pf.eta_c * pf.rho;
        ad_t area_bark = pf.a_b1 * al * pf.theta;
        ad_t mass_bark = area_bark * h * pf.eta_c * pf.rho;
        ad_t mass_root = pf.a_r1 * al;
        return dens * (mass_leaf + mass_sapwood + mass_bark + mass_root + mhw); });
    } else {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t { return dens * h; });
    }
  }

  Rcpp::NumericMatrix jac(M, T);
  Rcpp::NumericVector values(M);
  for (std::size_t m = 0; m < M; ++m) { values[m] = as_dbl(J[m]); tape.registerOutput(J[m]); }
  for (std::size_t m = 0; m < M; ++m) {
    tape.clearDerivatives();
    xad::derivative(J[m]) = 1.0;
    tape.computeAdjoints();
    for (std::size_t k = 0; k < T; ++k) jac(m, k) = xad::derivative(tr[k]);
  }
  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac, Rcpp::Named("values") = values);
}

namespace {
// Boundary new_node history per RK stage from an R list (lossy env list companion).
std::vector<std::vector<double>> nn_from_list(Rcpp::List l) {
  std::vector<std::vector<double>> v(l.size());
  for (R_xlen_t n = 0; n < l.size(); ++n) v[n] = Rcpp::as<std::vector<double>>(l[n]);
  return v;
}
} // namespace

// Fully native TF24f coupled (resident) census gradient: env + birth steps + boundary
// new_node history from the live Patch (no tf24f_harvest, no Rcpp::as<> env). species
// 0-based; birth_rate<0 recovers natively.
// [[Rcpp::export]]
Rcpp::List tf24f_coupled_gradient_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, double birth_rate, double k_acclim,
    bool use_ad_gradient, std::vector<std::string> traits,
    std::vector<std::string> metrics, double patch_area, double trait_rel_step) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::TF24f_Strategy, plant::TF24_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  if (patch.stand_newnode_height_stage_history.size() < 1)
    Rcpp::stop("the TF24f coupled (resident) census gradient needs the per-RK-stage "
               "boundary-node harvest; re-run with control(save_RK45_cache = TRUE)");
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  const std::vector<int> birth = tf24f_birth_steps(patch, (std::size_t)species);
  return tf24f_coupled_gradient_core(pp, patch.environment_history, patch.step_history,
    birth, br, k_acclim, use_ad_gradient, traits, metrics,
    patch.stand_newnode_height_stage_history,
    patch.stand_newnode_competition_stage_history, patch_area, trait_rel_step);
}


// ===========================================================================
// TF24f offspring_production reverse-mode AD trait gradient (#472 scope B, the
// offspring surface). The seed-rain integral offspring_production = sum_i tw_i *
// offspring_i, where offspring_i is the survival-weighted lifetime fecundity of
// cohort i (the FROZEN rare-mutant / invasion gradient, as for TF24/FF16). The
// scope (§5) calls this "inherits TF24's; minor (tracked-state seed)": it is the
// TF24 offspring tape (tf24_emergent.cpp) with the ONE TF24f difference the census
// tape already solves -- the tracked collar is a theta-dependent STATE that lags the
// optimum, so it is carried on the tape with the curvature-linearised gradient-ascent
// rate (NOT zeroed by the envelope theorem as TF24's optimised collar is). No canopy,
// no density / g': a per-cohort {5 demog, tracked collar, offspring} replay over the
// frozen schedule, offspring_dt = fecundity_dt * exp(-mortality) * (ppsurv/ppsab); one
// tape over all cohorts (offspring_production is linear in the per-cohort offspring).
// ===========================================================================

namespace {

template <typename S> struct St7off { plant::FF16State<S> demog; S collar; S offspring; };
template <typename S>
St7off<S> st7off_axpy(const St7off<S>& a, double c, const St7off<S>& k) {
  return St7off<S>{plant::FF16State<S>{
    a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
    a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
    a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar,
    a.offspring+c*k.offspring};
}

} // namespace


// Compiled core of tf24f_offspring_production_gradient(). Returns d(offspring_production)
// /d(trait) (a 1 x traits row) + the reconstructed value. Pass 1 harvests each cohort's
// per-stage leaf operating point (the frozen-census curvature harvest -- profit + its
// h/collar/theta sensitivities) along the frozen-env trajectory; pass 2 replays the
// {5 demog, collar, offspring} system onto ONE tape and takes one reverse sweep over the
// tw-weighted offspring sum. The tracked collar is carried as a taped state (the lag the
// envelope theorem does not zero), seeded at the birth optimum (IFT-injected).
// [[Rcpp::export]]
Rcpp::List tf24f_offspring_gradient_core(
    Rcpp::NumericVector pp,
    const std::vector<std::vector<plant::TF24_Environment>>& EH,
    std::vector<double> sh, std::vector<int> birth, Rcpp::NumericMatrix ppsurv,
    std::vector<double> ppsab, std::vector<double> tw, double k_acclim,
    bool use_ad_gradient, std::vector<std::string> traits, double trait_rel_step) {
  using std::exp; using std::log;
  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double recruitment_decay = s.pars.recruitment_decay;
  const std::size_t T = traits.size();

  Rcpp::CharacterVector ppn = pp.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in strategy pars): " + traits[k]);
    tr0[k] = pp[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(pp); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(pp); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim, use_ad_gradient));
    sm.push_back(make_tf24f(q2, k_acclim, use_ad_gradient));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  const std::size_t N = EH.size(), nC = birth.size();
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];
  auto env_at = [&](std::size_t n, int stage) -> const plant::TF24_Environment& {
    return (stage == 0) ? ((n > 0) ? EH[n - 1][5] : EH[0][0]) : EH[n][stage - 1];
  };

  // --- DISCOVERY: per-cohort frozen-env trajectory, harvesting the operating point
  // (no g' / backward point -- offspring has no log_density). ----------------------
  struct CohortHo { std::size_t b; double collar0; LH seed; std::vector<LH> stages; };
  auto harvest_cohort = [&](std::size_t i) -> CohortHo {
    CohortHo C; const std::size_t b = (std::size_t)birth[i]; C.b = b;
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    C.collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    C.seed = harvest_point(s, sp, sm, dk, eb, h0v, C.collar0);
    auto deriv = [&](const St6<double>& y, std::size_t n, int stage) -> St6<double> {
      const plant::TF24_Environment& e = env_at(n, stage);
      const double h = y.demog.height, collar = y.collar;
      C.stages.push_back(harvest_point(s, sp, sm, dk, e, h, collar));
      double dpsi = 0.0;
      const double net = net_at_tracked(s, e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      return St6<double>{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim * dpsi};
    };
    auto axpy = [](const St6<double>& a, double c, const St6<double>& k){ return st6_axpy<double>(a, c, k); };
    St6<double> y{plant::FF16State<double>{h0v, 0, 0, 0, 0}, C.collar0};
    plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    return C;
  };
  std::vector<CohortHo> CH(nC);
  for (std::size_t i = 0; i < nC; ++i) CH[i] = harvest_cohort(i);

  // --- Tape: replay {5 demog, collar, offspring} for each cohort; sum tw_i*offspring_i.
  auto replay_off = [&](std::size_t i, const std::vector<ad_t>& tr,
                        const plant::TF24ProdPars<ad_t>& pf) -> ad_t {
    const CohortHo& C = CH[i];
    const double ppsab_i = ppsab[i];
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
    ad_t collar0 = ad_t(C.collar0);
    const double inv_c2 = (C.seed.d2p_dpsi2 != 0.0) ? 1.0 / C.seed.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -C.seed.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0 += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    }
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0);
    ad_t profit0 = profit_lin(C.seed, h0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0, area0, profit0);
    ad_t uu = ad_t(a_d0) * area0 / net0;
    const double decay = std::exp(-recruitment_decay * sh[C.b]);
    ad_t mort0 = -log(ad_t(decay) / (uu * uu + ad_t(1.0)));

    std::size_t idx = 0;
    auto deriv = [&](const St7off<ad_t>& y, std::size_t n, int stage) -> St7off<ad_t> {
      const LH& Hh = C.stages[idx++];
      ad_t h = y.demog.height, collar = y.collar;
      ad_t profit = profit_lin(Hh, h, collar, tr, tr0);
      ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf, h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf, h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim) * dpsi_lin(Hh, h, collar, tr, tr0);
      ad_t off_dt = r.fecundity_dt * exp(-y.demog.mortality) * ad_t(ppsurv(n, stage) / ppsab_i);
      return St7off<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, off_dt};
    };
    auto axpy = [](const St7off<ad_t>& a, double c, const St7off<ad_t>& k){ return st7off_axpy<ad_t>(a, c, k); };
    St7off<ad_t> y{plant::FF16State<ad_t>{h0, mort0, ad_t(0), ad_t(0), ad_t(0)}, collar0, ad_t(0)};
    return plant::ff16_cashkarp_replay(y, step_h, C.b, deriv, axpy).offspring;
  };

  ad::tape_type tape;
  std::vector<ad_t> tr(T);
  for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
  for (auto& x : tr) tape.registerInput(x);
  tape.newRecording();
  plant::TF24ProdPars<ad_t> pf = pf_active(pd, traits, tr);
  ad_t acc = ad_t(0.0);
  for (std::size_t i = 0; i < nC; ++i) acc += ad_t(tw[i]) * replay_off(i, tr, pf);
  tape.registerOutput(acc);
  xad::derivative(acc) = 1.0;
  tape.computeAdjoints();

  Rcpp::NumericVector out(T);
  for (std::size_t k = 0; k < T; ++k) out[k] = xad::derivative(tr[k]);
  out.attr("names") = Rcpp::wrap(traits);
  out.attr("offspring_production") = as_dbl(acc);
  return Rcpp::List::create(Rcpp::Named("gradient") = out,
                            Rcpp::Named("value") = as_dbl(acc));
}

// [[Rcpp::export]]
Rcpp::List tf24f_offspring_gradient_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, double k_acclim, bool use_ad_gradient,
    std::vector<std::string> traits, double trait_rel_step) {
  return tf24f_offspring_gradient_core(pp, eh_from_list(eh_list), sh, birth, ppsurv,
    ppsab, tw, k_acclim, use_ad_gradient, traits, trait_rel_step);
}

// Fully native TF24f offspring (invasion/selection) gradient: env + the offspring
// harvest (birth steps, trapezoid weights, per-RK-stage + at-birth survival) built
// from the live Patch -- no tf24f_harvest, mirroring FF16's offspring path.
// [[Rcpp::export]]
Rcpp::List tf24f_offspring_gradient_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, double birth_rate, double k_acclim,
    bool use_ad_gradient, std::vector<std::string> traits, double trait_rel_step) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::TF24f_Strategy, plant::TF24_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  auto w = plant::gradient::offspring_weights(patch, (std::size_t)species,
                                              pp["S_D"], br);
  Rcpp::NumericMatrix ppsurv((int)w.N, 6);
  for (std::size_t k = 0; k < w.N; ++k)
    for (int st = 0; st < 6; ++st) ppsurv((int)k, st) = w.ppsurv[k * 6 + (std::size_t)st];
  return tf24f_offspring_gradient_core(pp, patch.environment_history, patch.step_history,
    w.birth, ppsurv, w.ppsab, w.tw, k_acclim, use_ad_gradient, traits, trait_rel_step);
}


// ===========================================================================
// TF24f MULTI-SPECIES coupled census reverse-mode AD trait gradient (#472 scope B,
// the cross-species resident Jacobian). All species' cohorts are re-evolved TOGETHER
// over the frozen schedule; the canopy light each RK stage is the JOINT reconstruction
//   competition(z) = (1/area) * sum_species [ trapezium over species s's active cohorts
//                    + its boundary node, with species s's eta ]
// (matching Patch::compute_competition = sum_s Species_s::compute_competition/area).
// Each cohort reads the joint canopy at ITS crown centre (h*eta_c of its species).
// Differentiating w.r.t. ONE species' traits gives the cross-species total: the target
// species' traits move its cohorts, re-shading the joint canopy that EVERY species reads,
// so the other species' contributions respond too (the cross term the frozen gradient
// zeroes). The TF24f mirror of ff16_emergent.cpp's assemble_metrics_coupled_ms; the leaf
// rate is the harvested+linearised tracked-collar leaf + the anchored crown-centre light
// channel (the same single-species machinery), with the theta-injection harvested ONLY
// for the target species (non-target cohorts respond to the canopy but carry no
// target-trait derivative -- their inj / d2p_dpsidth are zero).
// ===========================================================================

namespace {

// Joint canopy interpolator at frozen knots kx: per-species trapezium summed (each with
// its own eta / area_leaf / kI / boundary node). nnh[s]/nnc[s] = species s's boundary at
// this stage. Mirror of ff16's multi-species build_interp.
template <typename S>
odelia::interpolator::basic_interpolator<S>
build_canopy_ms(const std::vector<St7ad<S>>& stand, const std::vector<char>& alive,
                const std::vector<std::size_t>& gspec, const std::vector<double>& kx,
                const std::vector<double>& nnh, const std::vector<double>& nnc,
                const std::vector<plant::TF24ProdPars<S>>& pf,
                const std::vector<double>& kI, const std::vector<double>& eta,
                double area, std::size_t nS, std::vector<double>* ly_out = nullptr) {
  using std::exp;
  std::vector<std::vector<S>> hv(nS), gv(nS);
  for (std::size_t g = 0; g < stand.size(); ++g) if (alive[g]) {
    const std::size_t s = gspec[g];
    S hi = stand[g].demog.height;
    S dens = exp(stand[g].log_density);
    S al = plant::tf24_area_leaf<S>(pf[s].a_l1, pf[s].a_l2, hi);
    hv[s].push_back(hi); gv[s].push_back(dens * S(kI[s]) * al);
  }
  std::vector<S> comp(kx.size(), S(0.0));
  for (std::size_t s = 0; s < nS; ++s) {
    hv[s].push_back(S(nnh[s])); gv[s].push_back(S(nnc[s]));
    std::vector<std::size_t> ord(hv[s].size());
    for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(),
              [&](std::size_t a, std::size_t b){ return dval(hv[s][a]) > dval(hv[s][b]); });
    std::vector<S> hs(hv[s].size()), gs(hv[s].size());
    for (std::size_t k = 0; k < ord.size(); ++k) { hs[k] = hv[s][ord[k]]; gs[k] = gv[s][ord[k]]; }
    for (std::size_t k = 0; k < kx.size(); ++k)
      comp[k] = comp[k] + plant::gradient::canopy_comp_at<S>(kx[k], hs, gs, eta[s]);
  }
  std::vector<S> ly(kx.size());
  for (std::size_t k = 0; k < kx.size(); ++k) ly[k] = exp(-comp[k] * S(1.0 / area));
  if (ly_out) { ly_out->resize(kx.size());
    for (std::size_t k = 0; k < kx.size(); ++k) (*ly_out)[k] = dval(ly[k]); }
  odelia::interpolator::basic_interpolator<S> interp; interp.init(kx, ly);
  return interp;
}

// Harvest one cohort's operating point for the multi-species coupled tape: the leaf
// sensitivities (h / collar / curvature + light channel) via the cohort's OWN species
// strategy `s_sp`; the theta-injection (inj / d2p_dpsidth) only when this cohort belongs
// to the differentiated (target) species -- otherwise zeros sized to T_target.
LHc harvest_point_ms(plant::TF24f_Strategy& s_sp,
                     std::vector<plant::TF24f_Strategy>& sp_t,
                     std::vector<plant::TF24f_Strategy>& sm_t,
                     const std::vector<double>& dk, std::size_t T,
                     const plant::TF24_Environment& e, double h, double collar,
                     double eta_c, bool is_target) {
  LHc H;
  if (is_target) {
    H = harvest_point_c(s_sp, sp_t, sm_t, dk, e, h, collar, eta_c);
  } else {
    std::vector<plant::TF24f_Strategy> none; std::vector<double> noned;
    H.base = harvest_point(s_sp, none, none, noned, e, h, collar);
    H.base.inj.assign(T, 0.0); H.base.d2p_dpsidth.assign(T, 0.0);
    const double L0 = e.get_environment_at_height(h * eta_c);
    const double dL = 1e-5 * std::max(std::abs(L0), 1e-4);
    double pP, dpsiP, pM, dpsiM;
    leaf_at_light(s_sp, e, L0 + dL, h, collar, pP, dpsiP);
    leaf_at_light(s_sp, e, L0 - dL, h, collar, pM, dpsiM);
    H.dprofit_dL = (pP - pM) / (2.0 * dL);
    H.d2p_dpsidL = (dpsiP - dpsiM) / (2.0 * dL);
  }
  return H;
}

} // namespace


// Compiled core of the TF24f MULTI-SPECIES coupled census gradient. Returns {jacobian =
// metrics x traits (the target species' traits), values = TOTAL-stand metric values}.
// pp_list / birth_list / k_acclim / use_ad_gradient are per species; nn_h_list/nn_c_list
// the all-species per-RK-stage boundary harvest [step][stage][species]; target is 1-based.
// If `gate_only` (R0), runs a double joint re-evolution and returns values + env_err (the
// coupled-drift gauge) with no tape; else the AD cross-species sweep.
// Shared core of the TF24f cross-species coupled gradient (the analog of FF16's
// ff16_coupled_gradient_ms_core). Takes the joint env / step schedule / per-species
// cohort birth steps / all-species per-RK-stage boundary harvest as ready-built C++
// structures, so the R-list `_impl` (the FD-reference surface + lossy Rcpp::as<> env)
// and the live-Patch `_native` entry feed it identically. gate_only=TRUE runs the
// double R0 pass and returns {values, env_err}; FALSE runs the cross-species AD sweep
// and returns {jacobian, values}.
Rcpp::List tf24f_coupled_gradient_ms_core(
    Rcpp::List pp_list,
    const std::vector<std::vector<plant::TF24_Environment>>& EH,
    const std::vector<double>& sh,
    const std::vector<std::vector<int>>& birth,
    std::vector<double> birth_rate, std::vector<double> k_acclim,
    std::vector<int> use_ad_gradient, std::vector<std::string> traits,
    std::vector<std::string> metrics,
    const std::vector<std::vector<std::vector<double>>>& NNH,
    const std::vector<std::vector<std::vector<double>>>& NNC,
    double patch_area, int target, double trait_rel_step, bool gate_only) {
  using std::exp; using std::log;
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm);
  const std::size_t nS = pp_list.size();
  const std::size_t tgt = (std::size_t)(target - 1);
  if (tgt >= nS) Rcpp::stop("target species out of range");
  const std::size_t T = traits.size(), M = metrics.size();
  const double GEPS = 1e-6;

  // Per-species strategies + derived scalars.
  std::vector<plant::TF24f_Strategy> S_(nS);
  std::vector<plant::TF24ProdPars<double>> pd(nS);
  std::vector<double> h0(nS), a_d0(nS), kI(nS), eta(nS), eta_c(nS), rdecay(nS);
  for (std::size_t s = 0; s < nS; ++s) {
    Rcpp::NumericVector pp = pp_list[s];
    S_[s] = make_tf24f(pp, k_acclim[s], use_ad_gradient[s] != 0);
    pd[s] = S_[s].prod_pars();
    h0[s] = S_[s].initial_height(); a_d0[s] = S_[s].pars.a_d0;
    kI[s] = S_[s].pars.k_I; eta[s] = S_[s].pars.eta; eta_c[s] = S_[s].eta_c;
    rdecay[s] = S_[s].pars.recruitment_decay;
  }

  // Target species' base trait values + FD steps + perturbed strategies + IFT d(h0)/dth.
  Rcpp::NumericVector ppt = pp_list[(R_xlen_t)tgt];
  Rcpp::CharacterVector ppn = ppt.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in target pars): " + traits[k]);
    tr0[k] = ppt[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(ppt); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(ppt); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim[tgt], use_ad_gradient[tgt] != 0));
    sm.push_back(make_tf24f(q2, k_acclim[tgt], use_ad_gradient[tgt] != 0));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  // Joint env (EH) + per-stage boundary harvest (NNH/NNC, [step][stage][species]) arrive
  // ready-built: the `_impl` wrapper Rcpp::as<>-rebuilds them from R lists; the `_native`
  // entry borrows them from the live Patch. Only the per-step durations are derived here.
  const std::size_t N = EH.size();
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];
  auto env_idx = [&](std::size_t rn, int rs, std::size_t& en, int& es) {
    if (rs == 0) { if (rn > 0) { en = rn - 1; es = 5; } else { en = 0; es = 0; } }
    else { en = rn; es = rs - 1; }
  };

  // Flatten cohorts across species: gspec / glocal / gbirth (birth steps arrive built).
  std::vector<std::size_t> gspec, glocal;
  std::vector<int> gbirth;
  for (std::size_t s = 0; s < nS; ++s)
    for (std::size_t i = 0; i < birth[s].size(); ++i) {
      gspec.push_back(s); glocal.push_back(i); gbirth.push_back(birth[s][i]);
    }
  const std::size_t nG = gspec.size();

  // ---------- R0 gate: double joint re-evolution + env_err ----------
  if (gate_only) {
    std::vector<St7ad<double>> stand(nG);
    std::vector<char> alive(nG, 0);
    double env_err = 0.0;
    std::vector<plant::TF24ProdPars<double>> pfd = pd;   // double prod pars per species
    auto deriv = [&](const std::vector<St7ad<double>>& y, std::size_t rn, int rs)
        -> std::vector<St7ad<double>> {
      std::size_t en; int es; env_idx(rn, rs, en, es);
      const std::vector<double>& kx = EH[en][es].light_availability.spline.get_x();
      std::vector<double> ly0;
      auto interp = build_canopy_ms<double>(y, alive, gspec, kx, NNH[en][es], NNC[en][es],
                                            pfd, kI, eta, patch_area, nS, &ly0);
      const std::vector<double>& y0 = EH[en][es].light_availability.spline.get_y();
      for (std::size_t k = 0; k < ly0.size() && k < y0.size(); ++k)
        env_err = std::max(env_err, std::abs(ly0[k] - y0[k]));
      const double cap = kx.back();
      std::vector<St7ad<double>> dy(nG);
      for (std::size_t g = 0; g < nG; ++g) {
        if (!alive[g]) { dy[g] = St7ad<double>{plant::FF16State<double>{0,0,0,0,0},0,0}; continue; }
        const std::size_t s = gspec[g];
        const double h = y[g].demog.height, collar = y[g].collar;
        const double z = h * eta_c[s];
        const double L = (z > cap) ? 1.0 : interp(z);
        double profit, dpsi;
        leaf_at_light(S_[s], EH[en][es], L, h, collar, profit, dpsi);
        const double al = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h);
        const double net = plant::tf24_net_mass_production<double>(pd[s], h, al, profit);
        plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd[s], h, al, net, true);
        const double zb = (h - GEPS) * eta_c[s];
        const double Lb = (zb > cap) ? 1.0 : interp(zb);
        double pb, dpb; leaf_at_light(S_[s], EH[en][es], Lb, h - GEPS, collar, pb, dpb);
        const double alb = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h - GEPS);
        const double netb = plant::tf24_net_mass_production<double>(pd[s], h - GEPS, alb, pb);
        const double gback = plant::tf24_height_dt_from_net<double>(pd[s], h - GEPS, alb, netb);
        const double gprime = (r.height_dt - gback) / GEPS;
        dy[g] = St7ad<double>{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
          r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim[s] * dpsi, -gprime - r.mortality_dt};
      }
      return dy;
    };
    auto axpy = [](const std::vector<St7ad<double>>& a, double c,
                   const std::vector<St7ad<double>>& k){ return stand_axpy<double>(a, c, k); };
    auto establish = [&](std::size_t g, std::size_t bstep) {
      const std::size_t s = gspec[g];
      const std::size_t en = (bstep > 0) ? std::min(bstep, N) - 1 : 0;
      const plant::TF24_Environment& eb = (bstep > 0) ? EH[en][5] : EH[0][0];
      S_[s].initializing_ = true;
      S_[s].net_mass_production_dt(eb, h0[s], S_[s].area_leaf(h0[s]), 1.0 / h0[s]);
      const double collar0 = -S_[s].leaf.root_collar_psi_;
      S_[s].initializing_ = false;
      const double decay = std::exp(-rdecay[s] * sh[bstep]);
      double dpsi0; const double area0 = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h0[s]);
      const double net0 = net_at_tracked(S_[s], eb, h0[s], collar0, dpsi0);
      double pr_estab = 0.0;
      if (net0 > 0.0) { const double uu = a_d0[s] * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
      const double mort0 = (pr_estab > 0.0) ? -log(pr_estab) : std::numeric_limits<double>::infinity();
      const double g0 = plant::tf24_height_dt_from_net<double>(pd[s], h0[s], area0, net0);
      const double logd0 = (g0 > 0.0 && pr_estab > 0.0) ? log(birth_rate[s] * pr_estab / g0)
                                            : -std::numeric_limits<double>::infinity();
      stand[g] = St7ad<double>{plant::FF16State<double>{h0[s], mort0, 0, 0, 0}, collar0, logd0};
      alive[g] = 1;
    };
    for (std::size_t rn = 0; rn < N; ++rn) {
      for (std::size_t g = 0; g < nG; ++g) if ((std::size_t)gbirth[g] == rn) establish(g, rn);
      stand = rkck_one_step_tf(stand, step_h[rn], rn, deriv, axpy);
    }
    // Final-boundary cohorts (birth >= N): established at h0, never stepped (see the
    // single-species note) -- prevents the a_l2 0*log(0)=NaN from a zero census height.
    for (std::size_t g = 0; g < nG; ++g) if (!alive[g]) establish(g, (std::size_t)gbirth[g]);
    // Total-stand reduction: sum over species of that species' trapezium + tail.
    Rcpp::NumericVector values(M);
    for (std::size_t s = 0; s < nS; ++s) {
      std::vector<std::size_t> gs;
      for (std::size_t g = 0; g < nG; ++g) if (gspec[g] == s) gs.push_back(g);
      const std::size_t nc = gs.size();
      std::vector<std::size_t> ord(nc);
      for (std::size_t i = 0; i < nc; ++i) ord[i] = i;
      std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){
        return stand[gs[a]].demog.height > stand[gs[b]].demog.height; });
      const plant::TF24_Environment& ef = EH[N - 1][5];
      S_[s].initializing_ = true; S_[s].net_mass_production_dt(ef, h0[s], S_[s].area_leaf(h0[s]), 1.0 / h0[s]);
      const double collar0 = -S_[s].leaf.root_collar_psi_; S_[s].initializing_ = false;
      const double decay = std::exp(-rdecay[s] * sh[N]);
      double dpsi0; const double area0 = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h0[s]);
      const double net0 = net_at_tracked(S_[s], ef, h0[s], collar0, dpsi0);
      double pr_estab = 0.0;
      if (net0 > 0.0) { const double uu = a_d0[s] * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
      const double g0 = plant::tf24_height_dt_from_net<double>(pd[s], h0[s], area0, net0);
      const double dens_new = (g0 > 0.0 && pr_estab > 0.0) ? birth_rate[s] * pr_estab / g0 : 0.0;
      auto reduce = [&](auto psi) -> double {
        std::vector<double> phi(nc);
        for (std::size_t i = 0; i < nc; ++i)
          phi[i] = psi(stand[gs[i]].demog.height, std::exp(stand[gs[i]].log_density), stand[gs[i]].demog.mass_heartwood);
        double Js = 0.0;
        for (std::size_t j = 0; j + 1 < nc; ++j) {
          const std::size_t a = ord[j], b = ord[j + 1];
          Js += 0.5 * (stand[gs[a]].demog.height - stand[gs[b]].demog.height) * (phi[a] + phi[b]);
        }
        if (nc > 0) { const std::size_t last = ord[nc - 1];
          Js += 0.5 * (stand[gs[last]].demog.height - h0[s]) * (phi[last] + psi(h0[s], dens_new, 0.0)); }
        return Js;
      };
      for (std::size_t m = 0; m < M; ++m) {
        const std::string& nm = metrics[m];
        if (nm == "LAI") values[m] += reduce([&](double h, double dens, double){
          return dens * kI[s] * plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h); });
        else if (nm == "size_moment") values[m] += reduce([&](double h, double dens, double){ return dens * h; });
        else values[m] += reduce([&](double h, double dens, double mhw){
          const double al = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h);
          const double ml = al * pd[s].lma, as_ = al * pd[s].theta;
          const double ms = as_ * h * pd[s].eta_c * pd[s].rho;
          const double ab = pd[s].a_b1 * al * pd[s].theta, mb = ab * h * pd[s].eta_c * pd[s].rho;
          const double mr = pd[s].a_r1 * al;
          return dens * (ml + ms + mb + mr + mhw); });
      }
    }
    values.attr("names") = Rcpp::wrap(metrics);
    return Rcpp::List::create(Rcpp::Named("values") = values, Rcpp::Named("env_err") = env_err);
  }

  // ---------- R1: cross-species AD sweep ----------
  // DISCOVERY: per-cohort frozen-env harvest (target cohorts carry the theta channel).
  struct CohortHms { std::size_t s, b; double collar0; LHc seed; std::vector<StageHc> stages; };
  std::vector<CohortHms> CH(nG);
  for (std::size_t g = 0; g < nG; ++g) {
    const std::size_t s = gspec[g], b = (std::size_t)gbirth[g];
    const bool is_t = (s == tgt);
    CohortHms C; C.s = s; C.b = b;
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];
    S_[s].initializing_ = true;
    S_[s].net_mass_production_dt(eb, h0[s], S_[s].area_leaf(h0[s]), 1.0 / h0[s]);
    C.collar0 = -S_[s].leaf.root_collar_psi_;
    S_[s].initializing_ = false;
    // seed: a StageHc-less harvest (no g'); reuse harvest_point_ms (fwd point only).
    C.seed = harvest_point_ms(S_[s], sp, sm, dk, T, eb, h0[s], C.collar0, eta_c[s], is_t);
    auto deriv = [&](const St7& y, std::size_t n, int stage) -> St7 {
      std::size_t en; int es; env_idx(n, stage, en, es);
      const plant::TF24_Environment& e = EH[en][es];
      const double h = y.demog.height, collar = y.collar;
      C.stages.push_back(StageHc{
        harvest_point_ms(S_[s], sp, sm, dk, T, e, h, collar, eta_c[s], is_t),
        harvest_point_ms(S_[s], sp, sm, dk, T, e, h - GEPS, collar, eta_c[s], is_t)});
      double dpsi = 0.0;
      const double net = net_at_tracked(S_[s], e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd[s], h, al, net, true);
      double dpsi_b = 0.0;
      const double net_b = net_at_tracked(S_[s], e, h - GEPS, collar, dpsi_b);
      const double al_b = plant::tf24_area_leaf<double>(pd[s].a_l1, pd[s].a_l2, h - GEPS);
      const double g_back = plant::tf24_height_dt_from_net<double>(pd[s], h - GEPS, al_b, net_b);
      const double gprime = (r.height_dt - g_back) / GEPS;
      return St7{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim[s] * dpsi, -gprime - r.mortality_dt};
    };
    auto axpy = [](const St7& a, double c, const St7& k) -> St7 {
      return St7{plant::FF16State<double>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar, a.log_density+c*k.log_density};
    };
    St7 y{plant::FF16State<double>{h0[s], 0, 0, 0, 0}, C.collar0, 0};
    plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    CH[g] = C;
  }
  // Per-species pending-seed tail harvest at the final joint env.
  std::vector<double> collar0n(nS); std::vector<LHc> seed_new(nS); std::vector<double> decay_new(nS);
  const plant::TF24_Environment& ef = EH[N - 1][5];
  for (std::size_t s = 0; s < nS; ++s) {
    S_[s].initializing_ = true; S_[s].net_mass_production_dt(ef, h0[s], S_[s].area_leaf(h0[s]), 1.0 / h0[s]);
    collar0n[s] = -S_[s].leaf.root_collar_psi_; S_[s].initializing_ = false;
    seed_new[s] = harvest_point_ms(S_[s], sp, sm, dk, T, ef, h0[s], collar0n[s], eta_c[s], s == tgt);
    decay_new[s] = std::exp(-rdecay[s] * sh[N]);
  }

  // Tape: lift each species' prod_pars; register only the TARGET's traits.
  ad::tape_type tape;
  std::vector<ad_t> tr(T);
  for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
  for (auto& x : tr) tape.registerInput(x);
  tape.newRecording();
  std::vector<plant::TF24ProdPars<ad_t>> pf(nS);
  for (std::size_t s = 0; s < nS; ++s)
    pf[s] = (s == tgt) ? pf_active(pd[s], traits, tr) : pf_active(pd[s], {}, {});
  // target seedling height (IFT); others constant.
  std::vector<ad_t> h0a(nS);
  for (std::size_t s = 0; s < nS; ++s) h0a[s] = ad_t(h0[s]);
  for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0a[tgt] += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));

  std::vector<St7ad<ad_t>> stand(nG);
  std::vector<char> alive(nG, 0);
  std::vector<std::size_t> vidx(nG, 0);
  auto deriv = [&](const std::vector<St7ad<ad_t>>& y, std::size_t rn, int rs)
      -> std::vector<St7ad<ad_t>> {
    std::size_t en; int es; env_idx(rn, rs, en, es);
    const std::vector<double>& kx = EH[en][es].light_availability.spline.get_x();
    auto interp = build_canopy_ms<ad_t>(y, alive, gspec, kx, NNH[en][es], NNC[en][es],
                                        pf, kI, eta, patch_area, nS);
    const double cap = kx.back();
    std::vector<St7ad<ad_t>> dy(nG);
    for (std::size_t g = 0; g < nG; ++g) {
      if (!alive[g]) { dy[g] = St7ad<ad_t>{plant::FF16State<ad_t>{ad_t(0),ad_t(0),ad_t(0),ad_t(0),ad_t(0)},
        ad_t(0), ad_t(0)}; continue; }
      const std::size_t s = gspec[g];
      const StageHc& St = CH[g].stages[vidx[g]++];
      ad_t h = y[g].demog.height, collar = y[g].collar;
      ad_t z = h * ad_t(eta_c[s]); double zv = xad::value(z);
      ad_t dLr = (zv > cap) ? ad_t(0.0) : anchor(interp(zv));
      ad_t profit = profit_lin(St.fwd.base, h, collar, tr, tr0) + ad_t(St.fwd.dprofit_dL) * dLr;
      ad_t al = plant::tf24_area_leaf<ad_t>(pf[s].a_l1, pf[s].a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf[s], h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf[s], h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim[s]) * (dpsi_lin(St.fwd.base, h, collar, tr, tr0)
                                            + ad_t(St.fwd.d2p_dpsidL) * dLr);
      ad_t zb = (h - ad_t(GEPS)) * ad_t(eta_c[s]); double zbv = xad::value(zb);
      ad_t dLrb = (zbv > cap) ? ad_t(0.0) : anchor(interp(zbv));
      ad_t profit_b = profit_lin(St.back.base, h - ad_t(GEPS), collar, tr, tr0)
                      + ad_t(St.back.dprofit_dL) * dLrb;
      ad_t al_b = plant::tf24_area_leaf<ad_t>(pf[s].a_l1, pf[s].a_l2, h - ad_t(GEPS));
      ad_t net_b = plant::tf24_net_mass_production<ad_t>(pf[s], h - ad_t(GEPS), al_b, profit_b);
      ad_t g_back = plant::tf24_height_dt_from_net<ad_t>(pf[s], h - ad_t(GEPS), al_b, net_b);
      ad_t gprime = (r.height_dt - g_back) / ad_t(GEPS);
      dy[g] = St7ad<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, -gprime - r.mortality_dt};
    }
    return dy;
  };
  auto axpy = [](const std::vector<St7ad<ad_t>>& a, double c,
                 const std::vector<St7ad<ad_t>>& k){ return stand_axpy<ad_t>(a, c, k); };

  auto establish = [&](std::size_t g, std::size_t bstep) {
    const std::size_t s = gspec[g]; const CohortHms& C = CH[g];
    ad_t hh0 = h0a[s];
    ad_t collar0 = ad_t(C.collar0);
    const double inv_c2 = (C.seed.base.d2p_dpsi2 != 0.0) ? 1.0 / C.seed.base.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -C.seed.base.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0 += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    }
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf[s].a_l1, pf[s].a_l2, hh0);
    ad_t profit0 = profit_lin(C.seed.base, hh0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf[s], hh0, area0, profit0);
    ad_t uu = ad_t(a_d0[s]) * area0 / net0;
    const double decay = std::exp(-rdecay[s] * sh[bstep]);
    ad_t pr_estab = ad_t(decay) / (uu * uu + ad_t(1.0));
    ad_t mort0 = -log(pr_estab);
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf[s], hh0, area0, net0);
    ad_t logd0 = log(ad_t(birth_rate[s]) * pr_estab / g0);
    stand[g] = St7ad<ad_t>{plant::FF16State<ad_t>{hh0, mort0, ad_t(0), ad_t(0), ad_t(0)}, collar0, logd0};
    alive[g] = 1;
  };
  for (std::size_t rn = 0; rn < N; ++rn) {
    for (std::size_t g = 0; g < nG; ++g) if (CH[g].b == rn) establish(g, rn);
    stand = rkck_one_step_tf(stand, step_h[rn], rn, deriv, axpy);
  }
  // Final-boundary cohorts (birth >= N): established at h0, never stepped (see the AD
  // core note) -- else the zero census height poisons the shared tape via 0*log(0).
  for (std::size_t g = 0; g < nG; ++g) if (!alive[g]) establish(g, CH[g].b);

  // Total-stand reduction: sum over species of its trapezium + pending-seed tail.
  std::vector<ad_t> J(M, ad_t(0.0));
  for (std::size_t s = 0; s < nS; ++s) {
    std::vector<std::size_t> gs;
    for (std::size_t g = 0; g < nG; ++g) if (gspec[g] == s) gs.push_back(g);
    const std::size_t nc = gs.size();
    std::vector<std::size_t> ord(nc);
    for (std::size_t i = 0; i < nc; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){
      return as_dbl(stand[gs[a]].demog.height) > as_dbl(stand[gs[b]].demog.height); });
    // pending seed (frozen final joint env).
    ad_t hh0 = h0a[s];
    ad_t collar0 = ad_t(collar0n[s]);
    const double inv_c2 = (seed_new[s].base.d2p_dpsi2 != 0.0) ? 1.0 / seed_new[s].base.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -seed_new[s].base.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0 += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    }
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf[s].a_l1, pf[s].a_l2, hh0);
    ad_t profit0 = profit_lin(seed_new[s].base, hh0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf[s], hh0, area0, profit0);
    ad_t uu = ad_t(a_d0[s]) * area0 / net0;
    ad_t pr_estab = ad_t(decay_new[s]) / (uu * uu + ad_t(1.0));
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf[s], hh0, area0, net0);
    ad_t dens_new = ad_t(birth_rate[s]) * pr_estab / g0;
    auto reduce = [&](auto psi) -> ad_t {
      std::vector<ad_t> phi(nc);
      for (std::size_t i = 0; i < nc; ++i)
        phi[i] = psi(stand[gs[i]].demog.height, exp(stand[gs[i]].log_density), stand[gs[i]].demog.mass_heartwood);
      ad_t Js = ad_t(0.0);
      for (std::size_t j = 0; j + 1 < nc; ++j) {
        const std::size_t a = ord[j], b = ord[j + 1];
        Js += ad_t(0.5) * (stand[gs[a]].demog.height - stand[gs[b]].demog.height) * (phi[a] + phi[b]);
      }
      if (nc > 0) { const std::size_t last = ord[nc - 1];
        Js += ad_t(0.5) * (stand[gs[last]].demog.height - hh0) * (phi[last] + psi(hh0, dens_new, ad_t(0.0))); }
      return Js;
    };
    for (std::size_t m = 0; m < M; ++m) {
      const std::string& nm = metrics[m];
      if (nm == "LAI") J[m] += reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t {
        return dens * ad_t(kI[s]) * plant::tf24_area_leaf<ad_t>(pf[s].a_l1, pf[s].a_l2, h); });
      else if (nm == "size_moment") J[m] += reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t { return dens * h; });
      else J[m] += reduce([&](ad_t h, ad_t dens, ad_t mhw) -> ad_t {
        ad_t al = plant::tf24_area_leaf<ad_t>(pf[s].a_l1, pf[s].a_l2, h);
        ad_t ml = al * pf[s].lma, as_ = al * pf[s].theta;
        ad_t ms = as_ * h * pf[s].eta_c * pf[s].rho;
        ad_t ab = pf[s].a_b1 * al * pf[s].theta, mb = ab * h * pf[s].eta_c * pf[s].rho;
        ad_t mr = pf[s].a_r1 * al;
        return dens * (ml + ms + mb + mr + mhw); });
    }
  }

  Rcpp::NumericMatrix jac(M, T);
  Rcpp::NumericVector values(M);
  for (std::size_t m = 0; m < M; ++m) { values[m] = as_dbl(J[m]); tape.registerOutput(J[m]); }
  for (std::size_t m = 0; m < M; ++m) {
    tape.clearDerivatives();
    xad::derivative(J[m]) = 1.0;
    tape.computeAdjoints();
    for (std::size_t k = 0; k < T; ++k) jac(m, k) = xad::derivative(tr[k]);
  }
  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac, Rcpp::Named("values") = values);
}

// R-list entry (the FD-reference surface the test perturbs by re-passing pp_list, and
// the lossy Rcpp::as<>-env path): rebuild the joint env + boundary harvest + per-species
// birth steps from R lists, then call the shared core. Bit-identical to the `_native`
// entry below (same arithmetic, faithful data either way).
// [[Rcpp::export]]
Rcpp::List tf24f_coupled_gradient_ms_impl(
    Rcpp::List pp_list, Rcpp::List eh_list, std::vector<double> sh, Rcpp::List birth_list,
    std::vector<double> birth_rate, std::vector<double> k_acclim,
    std::vector<int> use_ad_gradient, std::vector<std::string> traits,
    std::vector<std::string> metrics, Rcpp::List nn_h_list, Rcpp::List nn_c_list,
    double patch_area, int target, double trait_rel_step, bool gate_only) {
  const std::size_t N = eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));
  }
  std::vector<std::vector<std::vector<double>>> NNH(N), NNC(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List hns = nn_h_list[n], cns = nn_c_list[n];
    const std::size_t ns = hns.size();
    NNH[n].resize(ns); NNC[n].resize(ns);
    for (std::size_t k = 0; k < ns; ++k) {
      NNH[n][k] = Rcpp::as<std::vector<double>>(hns[k]);
      NNC[n][k] = Rcpp::as<std::vector<double>>(cns[k]);
    }
  }
  const std::size_t nS = pp_list.size();
  std::vector<std::vector<int>> birth(nS);
  for (std::size_t s = 0; s < nS; ++s) birth[s] = Rcpp::as<std::vector<int>>(birth_list[s]);
  return tf24f_coupled_gradient_ms_core(pp_list, EH, sh, birth, birth_rate, k_acclim,
    use_ad_gradient, traits, metrics, NNH, NNC, patch_area, target, trait_rel_step,
    gate_only);
}

// FULLY native multi-species coupled cross-species gradient (the TF24f mirror of
// ff16_coupled_gradient_ms_native): the joint env + step schedule + per-species birth
// steps + all-species per-RK-stage boundary harvest come from the live Patch (no R
// tf24f_harvest_ms, no Rcpp::as<> env round-trip). pp_list / birth_rate / k_acclim /
// use_ad_gradient are cheap per-species scalars read from $parameters in R. gate_only
// switches the core between the cheap double R0 gate (returns env_err) and the AD sweep;
// the R wrapper runs the gate first then the sweep, exactly as the `_impl` path did.
// [[Rcpp::export]]
Rcpp::List tf24f_coupled_gradient_ms_native(
    SEXP scm_, Rcpp::List pp_list, std::vector<double> birth_rate,
    std::vector<double> k_acclim, std::vector<int> use_ad_gradient,
    std::vector<std::string> traits, std::vector<std::string> metrics,
    double patch_area, int target, double trait_rel_step, bool gate_only) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::TF24f_Strategy, plant::TF24_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  if (patch.stand_newnode_height_stage_history_all.size() < 1)
    Rcpp::stop("feedback = 'resident' on a multi-species TF24f stand needs the "
               "all-species per-RK-stage harvest; re-run the resident SCM with "
               "control(save_RK45_cache = TRUE)");
  const std::size_t nS = pp_list.size();
  std::vector<std::vector<int>> birth(nS);
  for (std::size_t s = 0; s < nS; ++s) birth[s] = plant::gradient::birth_steps(patch, s);
  return tf24f_coupled_gradient_ms_core(pp_list, patch.environment_history,
    patch.step_history, birth, birth_rate, k_acclim, use_ad_gradient, traits, metrics,
    patch.stand_newnode_height_stage_history_all,
    patch.stand_newnode_competition_stage_history_all, patch_area, target,
    trait_rel_step, gate_only);
}
