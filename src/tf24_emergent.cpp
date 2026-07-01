// Reverse-mode AD emergent-gradient routine for TF24, compiled into plant.so
// (#472 scope B, Phase F1-full). The TF24 counterpart of ff16_emergent.cpp.
//
// FF16 can tape the whole cohort trajectory directly because its net production is a
// closed form of light. TF24 cannot: net comes from the hydraulic LEAF OPTIMISATION
// (a maximisation over collar potential nesting a psi_stem->ci root-find), which has
// no adjoint tape. The resolution is a per-cohort HARVEST: a first (double) pass runs
// the real leaf opt at every RK stage and records the trait-INDEPENDENT operating
// point -- the optimised profit, the leaf sensitivities Leaf::dprofit_d*, k_max,
// E_up_, and the d(profit)/d(height) Jacobian. That turns the propagation into a
// leaf-opt-FREE, fully tapeable expression
//   profit(h, theta) = profit_0 + dprofit_dh*(h - h0_stage)
//                      + sum_k dprofit_dtheta_k * (theta_k - theta_k0),
// the cascade/area traits entering the committed kernel directly. ONE reverse sweep
// per cohort then yields d(w_i offspring_i)/d(all traits); summed over cohorts ->
// the emergent gradient. The XAD adjoint tape (xad::adj) is odelia's single compiled
// copy, resolved at load (see src/Makevars; odelia imported first).
//
// Differentiated: the full mass cascade + leaf path, the seedling size height_0 (by
// the implicit function theorem), and the recruitment filter (establishment, through
// the seedling net production). Resident light is frozen (the rare-mutant / invasion
// gradient). Takes the harvested data as plain arrays, carrying no templating into
// the SCM; the R-facing tf24_offspring_production_gradient() gathers these.
#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <XAD/XAD.hpp>
#include <plant.h>                 // RcppR6 as<>/wrap for TF24_Environment
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_production_kernel.h>
#include <plant/models/ff16_production_kernel.h>   // ff16_cashkarp_replay, FF16(Life)State
#include <plant/gradient/scm_harvest.h>            // shared native offspring harvest

using ad   = xad::adj<double>;
using ad_t = ad::active_type;

namespace {

const std::vector<std::string> TRAITS = {
  "vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec","curv_colim",
  "lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b","k_s","k_r","a_bio","a_y",
  "a_l1","a_l2","theta","a_r1"};
std::size_t IX(const std::string& n){
  for (std::size_t i=0;i<TRAITS.size();++i) if (TRAITS[i]==n) return i;
  return (std::size_t)-1;
}

plant::TF24_Strategy make_strategy(const Rcpp::NumericVector& pp,
                                   const std::string& over="", double v=0) {
  plant::TF24_Strategy s;
  s.control.shading_model="crown-centre"; s.control.GSS_tol_abs=1e-9;
  auto& q=s.pars;
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
  if(over=="g1_TF24")s.g1_TF24=v; else if(over=="curv_elec")q.curv_fact_elec_trans=v;
  else if(over=="curv_colim")q.curv_fact_colim=v; else if(over=="vcmax_25")q.vcmax_25=v;
  else if(over=="beta2")q.beta2=v; else if(over=="K_s")q.K_s=v; else if(over=="b")q.b=v;
  else if(over=="c")q.c=v; else if(over=="jmax_25")q.jmax_25=v; else if(over=="a")q.a=v;
  else if(over=="lma")q.lma=v; else if(over=="rho")q.rho=v; else if(over=="a_b1")q.a_b1=v;
  else if(over=="r_l")q.r_l=v; else if(over=="r_b")q.r_b=v; else if(over=="r_s")q.r_s=v;
  else if(over=="r_r")q.r_r=v; else if(over=="k_l")q.k_l=v; else if(over=="k_b")q.k_b=v;
  else if(over=="k_s")q.k_s=v; else if(over=="k_r")q.k_r=v; else if(over=="a_bio")q.a_bio=v;
  else if(over=="a_y")q.a_y=v; else if(over=="a_l1")q.a_l1=v; else if(over=="a_l2")q.a_l2=v;
  else if(over=="theta")q.theta=v; else if(over=="a_r1")q.a_r1=v;
  else if(!over.empty())Rcpp::stop("unknown TF24 trait: "+over);
  s.prepare_strategy(); return s;
}
double trait_value(const plant::TF24_Strategy& s, const std::string& t){
  if(t=="g1_TF24")return s.g1_TF24; if(t=="curv_elec")return s.pars.curv_fact_elec_trans;
  if(t=="curv_colim")return s.pars.curv_fact_colim; const auto&p=s.pars;
  if(t=="vcmax_25")return p.vcmax_25;if(t=="beta2")return p.beta2;if(t=="K_s")return p.K_s;
  if(t=="b")return p.b;if(t=="c")return p.c;if(t=="jmax_25")return p.jmax_25;if(t=="a")return p.a;
  if(t=="lma")return p.lma;if(t=="rho")return p.rho;if(t=="a_b1")return p.a_b1;if(t=="r_l")return p.r_l;
  if(t=="r_b")return p.r_b;if(t=="r_s")return p.r_s;if(t=="r_r")return p.r_r;if(t=="k_l")return p.k_l;
  if(t=="k_b")return p.k_b;if(t=="k_s")return p.k_s;if(t=="k_r")return p.k_r;if(t=="a_bio")return p.a_bio;
  if(t=="a_y")return p.a_y;if(t=="a_l1")return p.a_l1;if(t=="a_l2")return p.a_l2;if(t=="theta")return p.theta;
  if(t=="a_r1")return p.a_r1; Rcpp::stop("?"); return 0;
}
// TF24ProdPars<S> from the 27 trait scalars (cascade/area/leaf-coupled in; demographic
// fields frozen doubles from pd).
template <typename S, typename V>
plant::TF24ProdPars<S> pf_from(const V& tr, const plant::TF24ProdPars<double>& pd){
  plant::TF24ProdPars<S> p;
  p.lma=tr[IX("lma")];p.rho=tr[IX("rho")];p.theta=tr[IX("theta")];p.a_b1=tr[IX("a_b1")];
  p.a_r1=tr[IX("a_r1")];p.eta_c=S(pd.eta_c);
  p.r_l=tr[IX("r_l")];p.r_s=tr[IX("r_s")];p.r_b=tr[IX("r_b")];p.r_r=tr[IX("r_r")];
  p.k_l=tr[IX("k_l")];p.k_b=tr[IX("k_b")];p.k_s=tr[IX("k_s")];p.k_r=tr[IX("k_r")];
  p.a_bio=tr[IX("a_bio")];p.a_y=tr[IX("a_y")];p.a_l1=tr[IX("a_l1")];p.a_l2=tr[IX("a_l2")];
  p.a_f1=S(pd.a_f1);p.a_f2=S(pd.a_f2);p.hmat=S(pd.hmat);p.omega=S(pd.omega);p.a_f3=S(pd.a_f3);
  p.d_I=S(pd.d_I);p.a_dG1=S(pd.a_dG1);p.a_dG2=S(pd.a_dG2); return p;
}

// Trait-independent leaf-opt harvest at one operating point.
struct H {
  double h0_stage, profit, dprofit_dh, kmax, Eup;
  double dvcmax,dg1,dbeta2,dkmax,db,dc,djmax,da,dcelec,dccolim,dEup;
};
double profit_at(plant::TF24_Strategy& s, plant::TF24_Environment& e, double h){
  s.net_mass_production_dt(e,h,s.area_leaf(h),1.0/h); return s.leaf.profit_;
}
H harvest_at(plant::TF24_Strategy& s, plant::TF24_Environment& e, double h){
  H hh; hh.h0_stage=h; hh.profit=profit_at(s,e,h);
  const double opt=-s.leaf.root_collar_psi_;
  hh.kmax=s.leaf.leaf_specific_conductance_max_; hh.Eup=s.leaf.E_up_;
  hh.dvcmax=s.leaf.dprofit_dvcmax25(opt);hh.dg1=s.leaf.dprofit_dg1_TF24(opt);hh.dbeta2=s.leaf.dprofit_dbeta2(opt);
  hh.dkmax=s.leaf.dprofit_dkmax(opt);hh.db=s.leaf.dprofit_db(opt);hh.dc=s.leaf.dprofit_dc(opt);
  hh.djmax=s.leaf.dprofit_djmax25(opt);hh.da=s.leaf.dprofit_da(opt);hh.dcelec=s.leaf.dprofit_dcurv_elec(opt);
  hh.dccolim=s.leaf.dprofit_dcurv_colim(opt);hh.dEup=s.leaf.dprofit_dEup(opt);
  const double dd=1e-5*h;
  hh.dprofit_dh=(profit_at(s,e,h+dd)-profit_at(s,e,h-dd))/(2*dd);
  return hh;
}
// d(profit)/d(trait_k) injection coefficient from a harvest (0 for cascade-only).
double inj(std::size_t k, const H& h, const plant::TF24_Pars& p){
  const std::string& t=TRAITS[k];
  if(t=="vcmax_25")return h.dvcmax; if(t=="g1_TF24")return h.dg1; if(t=="beta2")return h.dbeta2;
  if(t=="b")return h.db; if(t=="c")return h.dc; if(t=="jmax_25")return h.djmax; if(t=="a")return h.da;
  if(t=="curv_elec")return h.dcelec; if(t=="curv_colim")return h.dccolim;
  if(t=="K_s")return h.dkmax*(h.kmax/p.K_s);
  if(t=="theta")return h.dkmax*(h.kmax/p.theta);
  if(t=="a_r1")return h.dEup*(h.Eup/p.a_r1);
  return 0.0;
}

// Pass-1 harvest of one cohort: integrate the demographic state in double (running
// the real leaf opt), recording the per-stage harvest. Mortality/offspring are not
// needed here (the leaf opt depends only on height + env).
void harvest_cohort(plant::TF24_Strategy& s, const plant::TF24ProdPars<double>& pd,
    std::vector<std::vector<plant::TF24_Environment>>& EH, std::size_t birth,
    const std::vector<double>& step_h, double h0, std::vector<H>& rec){
  rec.clear();
  auto deriv=[&](const plant::FF16State<double>& y,std::size_t n,int stage)->plant::FF16State<double>{
    plant::TF24_Environment* e=(stage==0)?((n>0)?&EH[n-1][5]:&EH[0][0]):&EH[n][stage-1];
    const double h=y.height; rec.push_back(harvest_at(s,*e,h));
    const H& hh=rec.back();
    double al=plant::tf24_area_leaf<double>(pd.a_l1,pd.a_l2,h);
    double net=plant::tf24_net_mass_production<double>(pd,h,al,hh.profit);
    plant::TF24Rates<double> r=plant::tf24_compute_rates_from_net<double>(pd,h,al,net,true);
    return plant::FF16State<double>{r.height_dt,r.mortality_dt,r.fecundity_dt,r.area_heartwood_dt,r.mass_heartwood_dt};
  };
  auto axpy=[](const plant::FF16State<double>&a,double c,const plant::FF16State<double>&k){
    return plant::FF16State<double>{a.height+c*k.height,a.mortality+c*k.mortality,a.fecundity+c*k.fecundity,a.area_heartwood+c*k.area_heartwood,a.mass_heartwood+c*k.mass_heartwood};};
  plant::FF16State<double> y{h0,0,0,0,0};
  plant::ff16_cashkarp_replay(y,step_h,birth,deriv,axpy);
}

// Pass-2: reverse-mode survival-weighted offspring of one cohort over its harvest.
// tr = the 27 active trait scalars; h0_init carries the height_0 IFT injection; Hs is
// the seedling harvest (for the establishment initial condition); decay/a_d0/ppsab fold
// to doubles. Returns the offspring (ad_t) for tw-weighting + the backward sweep.
ad_t offspring_ad(const std::vector<ad_t>& tr, const plant::TF24ProdPars<double>& pd,
    const plant::TF24_Pars& pars, const std::vector<H>& rec, const H& Hs, double decay,
    double a_d0, double ppsab, std::size_t birth, const std::vector<double>& step_h,
    ad_t h0_init, const Rcpp::NumericMatrix& ppsurv){
  using std::exp; using std::log;
  plant::TF24ProdPars<ad_t> pf = pf_from<ad_t>(tr, pd);
  // Establishment (recruitment filter), active: mort0 = -log(pr_estab), pr_estab from
  // the SEEDLING net production (leaf opt at height_0 in the birth env), so the trait
  // flows through net0 and area_leaf_0; height_0 itself is the IFT-injected h0_init.
  ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0_init);
  ad_t profit0 = Hs.profit + Hs.dprofit_dh*(h0_init - Hs.h0_stage);
  for (std::size_t k=0;k<tr.size();++k){ double c=inj(k,Hs,pars);
    if(c!=0.0) profit0 += c*(tr[k] - xad::value(tr[k])); }
  ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0_init, area0, profit0);
  ad_t uu = a_d0 * area0 / net0;
  ad_t mort0 = -log(decay / (uu*uu + 1.0));

  std::size_t idx=0;
  auto deriv=[&](const plant::FF16LifeState<ad_t>& y,std::size_t n,int stage)->plant::FF16LifeState<ad_t>{
    const H& hh=rec[idx++]; ad_t h=y.demog.height;
    ad_t profit = hh.profit + hh.dprofit_dh*(h - hh.h0_stage);
    for (std::size_t k=0;k<tr.size();++k){ double c=inj(k,hh,pars);
      if(c!=0.0) profit += c*(tr[k] - xad::value(tr[k])); }
    ad_t al=plant::tf24_area_leaf<ad_t>(pf.a_l1,pf.a_l2,h);
    ad_t net=plant::tf24_net_mass_production<ad_t>(pf,h,al,profit);
    plant::TF24Rates<ad_t> r=plant::tf24_compute_rates_from_net<ad_t>(pf,h,al,net,true);
    ad_t off_dt = r.fecundity_dt * exp(-y.demog.mortality) * ad_t(ppsurv(n,stage)/ppsab);
    return plant::FF16LifeState<ad_t>{plant::FF16State<ad_t>{r.height_dt,r.mortality_dt,
      r.fecundity_dt,r.area_heartwood_dt,r.mass_heartwood_dt}, off_dt};
  };
  auto axpy=[](const plant::FF16LifeState<ad_t>&a,double c,const plant::FF16LifeState<ad_t>&k)->plant::FF16LifeState<ad_t>{
    return plant::FF16LifeState<ad_t>{plant::FF16State<ad_t>{
      a.demog.height+c*k.demog.height,a.demog.mortality+c*k.demog.mortality,a.demog.fecundity+c*k.demog.fecundity,
      a.demog.area_heartwood+c*k.demog.area_heartwood,a.demog.mass_heartwood+c*k.demog.mass_heartwood},
      a.offspring+c*k.offspring};};
  plant::FF16LifeState<ad_t> y{plant::FF16State<ad_t>{h0_init,mort0,ad_t(0),ad_t(0),ad_t(0)}, ad_t(0)};
  return plant::ff16_cashkarp_replay(y,step_h,birth,deriv,axpy).offspring;
}

// ---------------------------------------------------------------------------
// Generic stand-gradient engine for TF24 (#472 scope B, build-order step 1 -- the
// TF24 mirror of ff16_emergent.cpp's engine). Same metric-symmetric design: a stand
// metric is a weighted reduction over the replayed cohort final states. TF24's net
// production has no adjoint tape (the hydraulic leaf optimisation), so the per-stage
// leaf OPERATING POINT is harvested in pass 1 and the trajectory replayed as a
// leaf-opt-free, profit-injected expression (as in offspring_ad). Unlike the
// offspring-only routine (which sweeps one tape per cohort), the census metrics
// (LAI/biomass/size_moment) couple cohorts through the height-trapezium, so the
// engine replays EVERY cohort onto ONE tape, then takes one reverse sweep per metric.
// ---------------------------------------------------------------------------

// Unified cohort state: 5 demog + offspring accumulator. (Census metrics would add
// log_density here, but TF24's census number density needs the SECOND-order leaf-opt
// sensitivity d(growth-rate-gradient)/d(trait): its growth gradient g' = d(height_dt)/
// d(height) is formed from the linearised leaf-opt harvest, which has no faithful
// d(g')/d(trait) -- so census metrics are deferred to a TF24 follow-up. The escape
// hatch (the 6 demographic components, all leaf-opt-free once harvested) is exact.)
template <typename S> struct Full { plant::FF16State<S> demog; S offspring; };

// Replay one cohort's Full state over its harvested operating points. profit_fn(h, hh)
// -> S supplies the leaf-opt-free profit at height h for stage harvest hh -- the TF24
// analogue of FF16's replay_cohort_full (offspring path).
template <typename S, typename ProfitFn>
Full<S> tf24_replay_full(const plant::TF24ProdPars<S>& pf, const std::vector<H>& rec,
    const H& Hs, double decay, double a_d0, double ppsab, std::size_t birth,
    const std::vector<double>& step_h, S h0_init, const Rcpp::NumericMatrix& ppsurv,
    double birth_rate, ProfitFn profit_fn) {
  using std::exp; using std::log;
  // Establishment (recruitment filter) -> initial mortality.
  S area0 = plant::tf24_area_leaf<S>(pf.a_l1, pf.a_l2, h0_init);
  S profit0 = profit_fn(h0_init, Hs);
  S net0 = plant::tf24_net_mass_production<S>(pf, h0_init, area0, profit0);
  S uu = a_d0 * area0 / net0;
  S mort0 = -log(decay / (uu * uu + 1.0));

  std::size_t idx = 0;
  auto deriv = [&](const Full<S>& y, std::size_t n, int stage) -> Full<S> {
    const H& hh = rec[idx++]; S h = y.demog.height;
    S profit = profit_fn(h, hh);
    S al = plant::tf24_area_leaf<S>(pf.a_l1, pf.a_l2, h);
    S net = plant::tf24_net_mass_production<S>(pf, h, al, profit);
    plant::TF24Rates<S> r = plant::tf24_compute_rates_from_net<S>(pf, h, al, net, true);
    S off_dt = r.fecundity_dt * exp(-y.demog.mortality) * S(ppsurv(n, stage) / ppsab);
    return Full<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
      r.area_heartwood_dt, r.mass_heartwood_dt}, off_dt};
  };
  auto axpy = [](const Full<S>& a, double c, const Full<S>& k) -> Full<S> {
    return Full<S>{plant::FF16State<S>{
      a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
      a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
      a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.offspring+c*k.offspring};
  };
  Full<S> y{plant::FF16State<S>{h0_init, mort0, S(0), S(0), S(0)}, S(0)};
  return plant::ff16_cashkarp_replay(y, step_h, birth, deriv, axpy);
}

} // namespace

// Compiled core of tf24_offspring_production_gradient(). Takes the harvested resident
// schedule (TF24_Environment per RK stage, step times sh), per-cohort birth steps /
// trapezoid-weights tw / survival ppsurv,ppsab, and the trait names. Returns
// d(offspring_production)/d(trait) (named) with attr "offspring_production" (the value
// reconstructed by the replay). Reverse mode, per cohort; establishment + height_0
// differentiated; resident light frozen.
Rcpp::NumericVector tf24_offspring_production_gradient_core(
    Rcpp::NumericVector pp,
    std::vector<std::vector<plant::TF24_Environment>> EH,   // by value: harvest_cohort needs a mutable stand
    std::vector<double> sh, std::vector<int> birth, Rcpp::NumericMatrix ppsurv,
    std::vector<double> ppsab, std::vector<double> tw, std::vector<std::string> traits) {
  plant::TF24_Strategy s0 = make_strategy(pp);
  plant::TF24ProdPars<double> pd = s0.prod_pars();
  const double h0v = s0.initial_height(), a_d0 = s0.pars.a_d0;

  const std::size_t N = EH.size();
  std::vector<double> step_h(N); for(std::size_t n=0;n<N;++n) step_h[n]=sh[n+1]-sh[n];

  const std::size_t T = TRAITS.size();
  std::vector<double> v0(T); for(std::size_t k=0;k<T;++k) v0[k]=trait_value(s0,TRAITS[k]);
  // d(height_0)/d(trait) by FD of initial_height (IFT of mass_live(h0)=omega); the
  // seedling root-find needs a 1e-4 step (a smaller step is corrupted by its noise).
  std::vector<double> dh0(T,0.0);
  for(std::size_t k=0;k<T;++k){ double dd=1e-4*std::abs(v0[k]);
    dh0[k]=(make_strategy(pp,TRAITS[k],v0[k]+dd).initial_height()
           -make_strategy(pp,TRAITS[k],v0[k]-dd).initial_height())/(2*dd); }

  std::vector<double> grad(T,0.0); double offspring=0.0;
  std::vector<H> rec;
  for (std::size_t i=0;i<birth.size();++i){
    const std::size_t b=(std::size_t)birth[i];
    plant::TF24_Environment& eb=(b>0)?EH[b-1][5]:EH[0][0];
    const double decay = std::exp(-s0.pars.recruitment_decay * sh[b]);
    // Harvest the trajectory + the seedling operating point (the leaf opts).
    harvest_cohort(s0, pd, EH, b, step_h, h0v, rec);
    H Hs = harvest_at(s0, eb, h0v);
    // ONE reverse sweep over all 27 traits for this cohort.
    ad::tape_type tape;
    std::vector<ad_t> tr(T); for(std::size_t k=0;k<T;++k) tr[k]=v0[k];
    for(auto& x:tr) tape.registerInput(x);
    tape.newRecording();
    ad_t h0_init = ad_t(h0v);
    for(std::size_t k=0;k<T;++k) if(dh0[k]!=0.0) h0_init += dh0[k]*(tr[k]-v0[k]);
    ad_t off = offspring_ad(tr, pd, s0.pars, rec, Hs, decay, a_d0, ppsab[i], b, step_h,
                            h0_init, ppsurv);
    ad_t Ji = tw[i] * off;
    tape.registerOutput(Ji); xad::derivative(Ji)=1.0; tape.computeAdjoints();
    offspring += xad::value(Ji);
    for(std::size_t k=0;k<T;++k) grad[k] += xad::derivative(tr[k]);
  }

  // Select the requested traits (the sweep covers all 27; cost is input-independent).
  Rcpp::NumericVector out(traits.size());
  for (std::size_t j=0;j<traits.size();++j){
    std::size_t k=IX(traits[j]);
    if (k==(std::size_t)-1) Rcpp::stop("unknown TF24 trait: "+traits[j]);
    out[j]=grad[k];
  }
  out.attr("names")=Rcpp::wrap(traits);
  out.attr("offspring_production")=offspring;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector tf24_offspring_production_gradient_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits) {
  const std::size_t N = eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n];
    for(R_xlen_t k=0;k<st.size();++k) EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));}
  return tf24_offspring_production_gradient_core(pp, EH, sh, birth, ppsurv, ppsab, tw, traits);
}

// Fully native TF24 offspring (invasion/selection) gradient: env + the offspring
// harvest (birth steps, trapezoid weights, per-RK-stage + at-birth survival) built
// from the live Patch -- no R harvest, mirroring FF16's offspring path.
// [[Rcpp::export]]
Rcpp::NumericVector tf24_offspring_production_gradient_native(
    SEXP scm_, Rcpp::NumericVector pp, int species, double birth_rate,
    std::vector<std::string> traits) {
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<
    plant::SCM<plant::TF24_Strategy, plant::TF24_Environment>>>(scm_);
  const auto& patch = scm->r_patch();
  double br = plant::gradient::recover_birth_rate(scm, (std::size_t)species, birth_rate);
  auto w = plant::gradient::offspring_weights(patch, (std::size_t)species, pp["S_D"], br);
  Rcpp::NumericMatrix ppsurv((int)w.N, 6);
  for (std::size_t k = 0; k < w.N; ++k)
    for (int st = 0; st < 6; ++st) ppsurv((int)k, st) = w.ppsurv[k * 6 + (std::size_t)st];
  return tf24_offspring_production_gradient_core(pp, patch.environment_history,
    patch.step_history, w.birth, ppsurv, w.ppsab, w.tw, traits);
}

namespace {
// Shared pass-1: harvest every cohort's per-stage leaf operating points + seedling
// harvests (the establishment initial condition). Fills recs/Hs/decay.
struct TF24Harvest {
  std::vector<std::vector<H>> recs; std::vector<H> Hs; std::vector<double> decay;
};
TF24Harvest tf24_pass1(plant::TF24_Strategy& s0, const plant::TF24ProdPars<double>& pd,
    std::vector<std::vector<plant::TF24_Environment>>& EH, const std::vector<int>& birth,
    const std::vector<double>& sh, const std::vector<double>& step_h, double h0v) {
  TF24Harvest H1; const std::size_t nC = birth.size();
  H1.recs.resize(nC); H1.Hs.resize(nC); H1.decay.resize(nC);
  for (std::size_t i = 0; i < nC; ++i) {
    const std::size_t b = (std::size_t)birth[i];
    plant::TF24_Environment& eb = (b > 0) ? EH[b-1][5] : EH[0][0];
    harvest_cohort(s0, pd, EH, b, step_h, h0v, H1.recs[i]);
    H1.Hs[i] = harvest_at(s0, eb, h0v);
    H1.decay[i] = std::exp(-s0.pars.recruitment_decay * sh[b]);
  }
  return H1;
}
} // namespace

// Escape hatch (TF24 mirror of ff16_state_jacobian_impl): per-cohort state x trait
// Jacobian. Each cohort is independent, so tape ONE cohort at a time, one reverse
// sweep per state component. Returns the cohort final states + the [cohort, component,
// trait] array.
// [[Rcpp::export]]
Rcpp::List tf24_state_jacobian_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits, double birth_rate) {
  plant::TF24_Strategy s0 = make_strategy(pp);
  plant::TF24ProdPars<double> pd = s0.prod_pars();
  const double h0v = s0.initial_height(), a_d0 = s0.pars.a_d0;
  const std::size_t N = eh_list.size(), nC = birth.size(), T = TRAITS.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n];
    for(R_xlen_t k=0;k<st.size();++k) EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));}
  std::vector<double> step_h(N); for(std::size_t n=0;n<N;++n) step_h[n]=sh[n+1]-sh[n];
  std::vector<double> v0(T); for(std::size_t k=0;k<T;++k) v0[k]=trait_value(s0,TRAITS[k]);
  std::vector<double> dh0(T,0.0);
  for(std::size_t k=0;k<T;++k){ double dd=1e-4*std::abs(v0[k]);
    dh0[k]=(make_strategy(pp,TRAITS[k],v0[k]+dd).initial_height()
           -make_strategy(pp,TRAITS[k],v0[k]-dd).initial_height())/(2*dd); }
  TF24Harvest H1 = tf24_pass1(s0, pd, EH, birth, sh, step_h, h0v);

  // Map requested traits to indices.
  std::vector<std::size_t> col(traits.size());
  for(std::size_t j=0;j<traits.size();++j){ col[j]=IX(traits[j]);
    if(col[j]==(std::size_t)-1) Rcpp::stop("unknown TF24 trait: "+traits[j]); }
  const std::vector<std::string> comp =
    {"height","mortality","fecundity","area_heartwood","mass_heartwood","offspring"};
  const std::size_t nS = comp.size(), nT = traits.size();
  Rcpp::NumericMatrix states(nC, nS);
  Rcpp::NumericVector jac(nC*nS*nT);
  auto JAC=[&](std::size_t i,std::size_t c,std::size_t k)->double&{ return jac[i+nC*(c+nS*k)]; };

  for (std::size_t i=0;i<nC;++i){
    ad::tape_type tape;
    std::vector<ad_t> tr(T); for(std::size_t k=0;k<T;++k) tr[k]=v0[k];
    for(auto& x:tr) tape.registerInput(x);
    tape.newRecording();
    plant::TF24ProdPars<ad_t> pf = pf_from<ad_t>(tr, pd);
    ad_t h0_init = ad_t(h0v);
    for(std::size_t k=0;k<T;++k) if(dh0[k]!=0.0) h0_init += dh0[k]*(tr[k]-v0[k]);
    auto profit_ad = [&](ad_t h, const H& hh) -> ad_t {
      ad_t p = ad_t(hh.profit) + ad_t(hh.dprofit_dh)*(h - ad_t(hh.h0_stage));
      for(std::size_t k=0;k<T;++k){ double c=inj(k,hh,s0.pars);
        if(c!=0.0) p += ad_t(c)*(tr[k]-ad_t(xad::value(tr[k]))); }
      return p; };
    Full<ad_t> y = tf24_replay_full<ad_t>(pf, H1.recs[i], H1.Hs[i], H1.decay[i], a_d0,
      ppsab[i], (std::size_t)birth[i], step_h, h0_init, ppsurv, birth_rate, profit_ad);
    ad_t out[6] = {y.demog.height, y.demog.mortality, y.demog.fecundity,
                   y.demog.area_heartwood, y.demog.mass_heartwood, y.offspring};
    for(std::size_t c=0;c<nS;++c){ states(i,c)=xad::value(out[c]); tape.registerOutput(out[c]); }
    for(std::size_t c=0;c<nS;++c){ tape.clearDerivatives(); xad::derivative(out[c])=1.0;
      tape.computeAdjoints();
      for(std::size_t j=0;j<nT;++j) JAC(i,c,j)=xad::derivative(tr[col[j]]); }
  }
  states.attr("dimnames")=Rcpp::List::create(R_NilValue, Rcpp::wrap(comp));
  jac.attr("dim")=Rcpp::IntegerVector::create((int)nC,(int)nS,(int)nT);
  jac.attr("dimnames")=Rcpp::List::create(R_NilValue, Rcpp::wrap(comp), Rcpp::wrap(traits));
  return Rcpp::List::create(Rcpp::Named("states")=states, Rcpp::Named("jacobian")=jac);
}
