// -*-c++-*-
// TF24_Model's rates: the ODE right-hand side, respiration and turnover,
// allocation, fecundity, mortality, storage, establishment and competition.
// The class is declared in tf24_model.h, which includes this file.
#ifndef PLANT_PLANT_TF24_RATES_H_
#define PLANT_PLANT_TF24_RATES_H_

#include <plant/models/tf24_model.h>

namespace plant {

// ⚠️ DO NOT HOIST the per-second -> annual factor 60*60*12*365 (seconds of
// daylight per year, 12 h day x 365 d) into one constant. It recurs in
// compute_rates and sapwood_marginal_return below and in
// net_mass_production_dt (tf24_leaf.h), each time after a double, so C++'s left-to-right evaluation makes it four rounded
// double multiplications rather than one integer; a single constant rounds
// differently, and the adaptive ODE amplifies that (offspring_production
// shifts ~0.2%).

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <typename S>
void TF24_Model<S>::compute_rates(const TF24_Environment<S>& environment,  Internals<S>& vars) {
  const S& height = vars.state(HEIGHT_INDEX);
  const S& area_leaf_ = vars.aux(aux_idx_competition_effect);
  // The Huber value's departure from theta. Read once here and threaded through
  // every sapwood quantity below, so the conducting area a plant has and the
  // one its leaf area prefers cannot disagree between two call sites.
  const S& sapwood_departure = vars.state(state_idx_log_area_sapwood_departure);
  const S& leaf_departure = vars.state(state_idx_log_area_leaf_departure);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse),
                           sapwood_departure, leaf_departure,
                           vars.state(state_idx_leaf_marginal_tracked));

  // The running mean chases the instantaneous balance at a_pl4. A first-order
  // filter, so the memory is 1/a_pl4 and nothing about it can overshoot.
  vars.set_aux(aux_idx_leaf_marginal_return, marginal_leaf_return_);
  vars.set_rate(state_idx_leaf_marginal_tracked,
                pars.a_pl4 * (marginal_leaf_return_ -
                              vars.state(state_idx_leaf_marginal_tracked)));

  // Read back the shares the budget above actually charged, rather than forming
  // them again: the carbon not spent and the tissue lost are two halves of one
  // decision, and they balance only if they come from one number.
  const S withheld         = withheld_leaf_;
  const S withheld_sapwood = withheld_sapwood_;
  const S replacement      = 1.0 - withheld;

  // Read by establishment_probability for the boundary node, so this one is a
  // rate's input rather than a reading of it.
  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);

  write_leaf_readings(vars, area_leaf_);

  // consumption rates should be emerging from net_mass_produciton_dt
  // convert evapotranspiration per leaf area per soil layer (mol H20 m^-2 s^-1) to canopy-level total 
  // yearly evapotranspiration per soil layer (m yr^-1)
  // stubbing out E_p for integration
  int soil_number_of_depths_ = environment.get_soil_number_of_depths();


  for (int i = 0; i < soil_number_of_depths_; i++) {

    // evapotranspiration (mol H20 m^-2 s^-1 layer^-1)
    // consumption rate (m yr^-1 layer ^-1)
    vars.set_consumption_rate(i, evapotranspiration_dt(area_leaf_, i)*60*60*12*365/1000*kg_per_mol_h2o);
  }

  // A storage pool buffers demography against short-term productivity swings:
  // growth is gated on having ample reserves, and mortality reads the buffered
  // relative reserves rather than instantaneous net production -- so a trough
  // draws reserves down and death is gradual, instead of the growth cutoff and
  // ~1e32 mortality spike that caused the #550 blow-up.
  // Read unclamped: a stage below empty keeps finite rates for the error
  // estimate, and non_negative_states() refuses a step that ends there.
  const S storage     = vars.state(state_idx_storage);
  const S storage_max = storage_capacity(area_leaf_, height, sapwood_departure);
  using std::exp;
  using std::sqrt;
  const S r           = storage_max > 0.0 ? storage / storage_max : S(0.0);
  // Reserve-gated growth (#517), following Daniel's intuition that a plant
  // should not grow unless it has ample carbon in storage. Growth and
  // reproduction proceed at the *production* rate, but scaled by a smooth gate
  // G(r) that is ~0 at low relative reserves and ~1 near capacity (logistic
  // centred on the growth threshold a_st2). Carbon not spent on growth refills
  // storage (dS/dt below), so: full reserves -> grow at production rate (healthy
  // dynamics preserved); low reserves -> redirect carbon to refill, growth
  // pauses; drought (net<0) -> reserves drain, growth halts and only resumes
  // once they refill (buffered growth AND mortality). Note growth is gated at
  // the production rate rather than metered out of the pool as a_st2*S: with
  // sapwood-scaled capacity, storage is tiny relative to seedling productivity,
  // so a discharge-rate limit would choke establishment.
  const S G = 1.0 / (1.0 + exp(-(r - pars.a_st2) / storage_gate_width));
  // Smooth positive part of net production (replaces the old hard net>0 cutoff).
  const S P = net_mass_production_dt_;
  const S Ppos =
    0.5 * (P + sqrt(P * P + storage_prod_eps * storage_prod_eps));
  const S growth_flux = Ppos * G;

  const S fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
  const S fraction_allocation_growth_ = fraction_allocation_growth(height);

  // --- Where the growth flux goes (#516) ------------------------------------
  // Growth moves a plant ALONG its allometry, plasticity moves it OFF, and
  // rebuilding moves it back. Keeping those three separate is what makes the
  // resting model exact and the departed one interpretable.
  //
  // So the allometric derivatives are evaluated at the leaf area the plant's
  // HEIGHT implies, area_leaf(height), not at the leaf area it actually has.
  // They are properties of the curve, and a plant's place on the curve is given
  // by its height; its leaf area is a departure from that. Evaluated this way,
  // extension is exactly departure-preserving for any departure, so the whole
  // of the departure's motion is plasticity. At rest the two coincide bitwise
  // (exp(0) == 1), so this is the same number the fixed model used.
  const S area_leaf_star = area_leaf(height);
  const S darea_leaf_dmass_live_ =
    darea_leaf_dmass_live(area_leaf_star, sapwood_departure);

  // The share of growth redirected into closing the canopy gap: it vanishes
  // exactly when the gap closes, widens as the gap does, and is gated by the
  // SAME signal that decides replacement.
  //
  // ⚠️ THE `replacement` FACTOR IS LOAD-BEARING, not a refinement. Without it
  // the gap alone drives rebuilding, so a starving plant withholds replacement
  // and spends its growth flux rebuilding at the same time -- contradictory,
  // since carbon enough to rebuild is carbon enough to maintain. Measured, the
  // unconditioned form rebuilt at ~1.5/yr against shedding's 0.42/yr, so the
  // canopy thinned by 5 per cent and stalled, and the Huber value fell instead
  // of rising because rebuilding adds leaf area without sapwood. Gating both on
  // the same signal makes a plant rebuild exactly to the degree it is willing
  // to maintain.
  //
  // The departure cannot go positive under these dynamics -- shedding only
  // lowers it and this share only raises it, stopping at zero -- so the form
  // needs no bound.
  const S rebuild_share =
    replacement * (1.0 - exp(leaf_departure / pars.a_pl2));
  const S extension_share = 1.0 - rebuild_share;

  const S area_leaf_dt =
    extension_share * growth_flux * fraction_allocation_growth_ *
      darea_leaf_dmass_live_;

  // Height follows extension alone, so rebuilding a canopy buys no height and
  // height still cannot fall: every factor here is non-negative.
  vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_star) * area_leaf_dt);
  vars.set_rate(FECUNDITY_INDEX,
    fecundity_dt(growth_flux, fraction_allocation_reproduction_));

  // Rebuilding buys leaf area at unchanged height, and it buys the WHOLE
  // package -- leaf, fine root, and the sapwood and bark that supply it -- at
  // the ratio the plant prefers rather than the one it has. Buying leaf alone
  // diluted the Huber value, and measured, that dilution dominated: plants
  // finished under-built rather than over-built (#516).
  const S rebuild_rel =
    rebuild_share * growth_flux * fraction_allocation_growth_ /
      (mass_live_per_area_leaf_at_height(height) * area_leaf_);

  // Sapwood -> heartwood conversion is turnover-driven, so it proceeds
  // regardless of carbon status. DO NOT gate it behind net production > 0.
  vars.set_rate(state_idx_area_heartwood,
                area_heartwood_dt(area_leaf_, sapwood_departure));
  const S area_sapwood_ = area_sapwood(area_leaf_, sapwood_departure);
  const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  vars.set_rate(state_idx_mass_heartwood, mass_heartwood_dt(mass_sapwood_));

  if (this->collect_all_auxiliary) {
    vars.set_aux(aux_idx_area_sapwood, area_sapwood_);
  }

  // Storage dynamics, as a charge and a drain that are each non-negative
  // without a test: sqrt(P^2 + eps^2) >= |P| makes Ppos and the drain both >= 0,
  // and G is a logistic so 1 - G >= 0. Their difference is the net flux exactly,
  // Ppos(1 - G) - (Ppos - P) = P - growth_flux, so the split itself moves
  // nothing; what it buys is somewhere to limit each direction on its own
  // (#609).
  const S charge = Ppos * (1.0 - G);      // surplus the gate withheld
  const S drain  = Ppos - P;              // shortfall met from reserves
  // Each direction is limited by the room the other has: the charge fills
  // headroom, the drain spends contents. So dS/dt = charge - (charge + drain) r,
  // the pool is a first-order filter on production, and both bounds follow from
  // the form with no scale left to choose -- at r = 0 the rate is charge >= 0,
  // at r = 1 it is -drain <= 0. It is also smooth where the net flux changes
  // sign, where the old branch stepped the slope by (r + drain_ref)/r, unbounded
  // as the reserves empty.
  //
  // Narrowing the drain to r/(r + D) is the shape this replaced and the one a
  // reader is most likely to restore, since it looks like the more faithful
  // limiter. It puts an attracting fixed point at r ~ D whose relaxation time is
  // storage_max D / drain -- under an hour for a seedling at D = 1e-3, against a
  // solver stepping in days, and measured at 26 times the accepted steps.
  //
  // The pool is capped by withholding the surplus rather than by spending it,
  // so production and the two flows no longer balance: at capacity the charge
  // the gate withheld, Ppos(1 - G) ~ 1.2e-4 of production, leaves the budget.
  //
  // The flow relaxes the pool over storage_max/(charge + drain), hours in a
  // seedling; the divisor adds storage_relaxation_offset to that time.
  const S relaxation_rate =
      storage_max > 0.0 ? (charge + drain) / storage_max : S(0.0);
  vars.set_rate(state_idx_storage,
                (charge * (1.0 - r) - drain * r) /
                    (1.0 + relaxation_rate * pars.storage_relaxation_offset));

  // --- The two departures (#516) --------------------------------------------
  // Relative rates, so both are scale-free: the same withheld fraction thins a
  // seedling's canopy and an adult's by the same proportion, which is what
  // keeps these parameters from repeating storage_prod_eps's mistake (#620) of
  // being an absolute constant against fluxes spanning six orders.
  //
  // Leaf area loses the turnover it did not replace and regains what rebuilding
  // buys. Thinning is therefore bounded by k_l -- a plant declines to replace
  // what died rather than actively shedding -- so the canopy cannot empty
  // faster than its leaves turn over.
  vars.set_rate(state_idx_log_area_leaf_departure,
                rebuild_rel - withheld * pars.k_l);

  const S sapwood_marginal_ = sapwood_marginal_return(
      area_leaf_, area_leaf_star, sapwood_departure, mass_sapwood_,
      withheld_sapwood, growth_flux, net_mass_production_dt_,
      darea_leaf_dmass_live_);
  vars.set_aux(aux_idx_sapwood_marginal_return, sapwood_marginal_);

  // The Huber value's departure is the DIFFERENCE of the two losses, so the
  // drought acclimation falls out of the arithmetic rather than being imposed:
  // leaf area thins at up to k_l = 0.457/yr while conducting area is only lost
  // to heartwood at k_s = 0.2/yr, and withholding both raises sapwood per leaf
  // area. The excess is then renewed less (withheld_sapwood > withheld), which
  // is the restoring force.
  //
  // The rebuilding term is exact geometry, not a parameter: rebuilding adds leaf
  // area and sapwood in the ratio theta, so it moves the two logs by b/A and
  // theta*b/A_s, whose difference is (b/A) * (exp(-psi) - 1). It vanishes at
  // psi = 0, which is what leaves the inward-pointing boundary flow intact.
  vars.set_rate(state_idx_log_area_sapwood_departure,
                withheld * pars.k_l - withheld_sapwood * pars.k_s +
                  rebuild_rel * (exp(-sapwood_departure) - 1.0) +
                  pars.a_sw * sapwood_marginal_);

  // Instantaneous mortality rate, now driven by relative reserves r.
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(r, vars.state(MORTALITY_INDEX)));

}

// The leaf's outputs and the root mass, written as readings: their only
// consumers are r_internals and ode_aux, both R-facing. Writing them at double alone was tried and measured at no
// difference, because a store into a slot the vector already holds pushes a
// statement and registers nothing -- it is a newly constructed active scalar
// that costs a slot, not a write to a live one.
template <typename S>
void TF24_Model<S>::write_leaf_readings(Internals<S>& vars,
                                           const S& area_leaf_) {
  vars.set_aux(aux_idx_root_mass, mass_root(area_leaf_));
  vars.set_aux(aux_idx_opt_psi_leaf, leaf.opt_psi_leaf_);
  vars.set_aux(aux_idx_opt_root_psi, leaf.opt_root_psi_);
  vars.set_aux(aux_idx_transpiration, leaf.transpiration_);
  vars.set_aux(aux_idx_E_up, leaf.E_up_);
  vars.set_aux(aux_idx_profit, leaf.profit_);
  // Written every step. At the default price the right value is zero, which is
  // also what a slot nobody writes reads, so a missing write here would show
  // only once a price is set.
  vars.set_aux(aux_idx_shadow_cost, S(leaf.shadow_cost()));
  vars.set_aux(aux_idx_stom_cond_CO2, leaf.stom_cond_CO2_);
  vars.set_aux(aux_idx_assimilation, leaf.assim_colimited_);
  // The leaf's own temperature at the operating point, not the `leaf_temp`
  // driver -- see aux_names(). Equal to the driver while pars.use_energy_balance
  // is off; solved from the leaf's transpiration when it is on.
  vars.set_aux(aux_idx_Tleaf, leaf.Tleaf_);
}

// --- Sapwood chases the ratio that maximises GROWTH -------------------------
// An integral controller on the marginal return: dpsi/dt = a_sw * R_s, which
// converges exactly onto R_s = 0 and needs no separate tracked state, because
// integrating IS averaging.
//
// R_s is the marginal return of sapwood area, in growth rather than in
// production. Growth is F * f_g * D with D = darea_leaf_dmass_live, so
//     d(growth)/dpsi  proportional to  dP/dpsi - P * (dm_sap/dA) * D
// the second term being the price of the extra stem, paid in forgone leaf
// area. Dropping it optimises production instead, which peaks near 3.3x the
// pipe-model ratio against growth's 1.4-1.8x.
template <typename S>
S TF24_Model<S>::sapwood_marginal_return(
    const S& area_leaf_, const S& area_leaf_star, const S& sapwood_departure,
    const S& mass_sapwood_, const S& withheld_sapwood, const S& growth_flux,
    const S& net_mass_production_dt_, const S& darea_leaf_dmass_live_) {
  S sapwood_marginal_ = 0.0;
  if (odelia::util::to_passive(pars.a_sw) > 0.0 &&
      maintenance_flux_ > 0.0) {
    const S mass_sapwood_now = mass_sapwood_;
    // Read, do not recompute: #617 replaced the height-linear denominator with
    // a path integral, and a recomputed copy here would still have been on the
    // old relation -- wrong by a factor 0.35 at 16 m, and silent.
    const S kmax_now = leaf_specific_conductance_max_;
    const S conv = pars.a_bio * pars.a_y;
    // psi enters kmax and sapwood mass both as exp(psi), so d/dpsi of each is
    // itself -- which is what makes these two terms so simple.
    // ⚠️ SAPWOOD TURNOVER IS GATED, so its psi-derivative is not k_s * m_s. What
    // the budget is actually charged is
    //     replacement_sapwood * k_s * m_s,  replacement_sapwood ∝ exp(-psi/a_pl2)
    // against m_s ∝ exp(psi), so the charge falls as exp((1 - 1/a_pl2) * psi) --
    // steeply, a_pl2 being well below one. Differentiating the UNGATED cost
    // instead drops a term worth +6.8 kg/yr of 13.1 at the default, and moved
    // the controller's zero from 1.35x the pipe-model ratio down to 1.05x.
    // Respiration carries no such gate, so it keeps the plain exp(psi).
    const S charged_sapwood_turnover =
      (1.0 - withheld_sapwood) * pars.k_s * mass_sapwood_now;
    const S dP_dpsi =
      conv * area_leaf_ * dprofit_dkmax_ * kmax_now * 60*60*12*365/1e6
      - conv * pars.r_s * mass_sapwood_now
      - charged_sapwood_turnover * (1.0 - 1.0 / pars.a_pl2);
    const S dmass_sap_dA =
      dmass_sapwood_darea_leaf(area_leaf_star, sapwood_departure);
    // ⚠️ SCALED BY CARBON AVAILABILITY, because re-proportioning the stem is an
    // allocation decision and has to be paid for out of the growth flux.
    // Without this the drift raises area_sapwood at fixed leaf area -- real
    // tissue -- while the budget is charged nothing, so a plant in deficit would
    // keep thickening on carbon it does not have. The reserve gate G is NOT the
    // right switch on its own: at soil 0.13 reserves are still full so G ~ 1,
    // and what has gone to zero is Ppos. growth_flux carries both.
    //
    // The factor is a share of throughput, so it is smooth, lies in [0, 1), and
    // vanishes with the growth flux. Symmetric on purpose: scaling only the
    // building direction would put a kink at R_s = 0, and the turnover route
    // out of an over-built stem stays open regardless, being a real mass flow.
    // Being strictly positive it cannot move the zero, so the controller still
    // converges on the growth optimum -- it only sets how fast, which is what
    // a_sw is for.
    //
    // Divided by the maintenance flux to make it dimensionless: strictly
    // positive, unlike P, so the sign of R_s is the sign of the growth
    // derivative everywhere including deep deficit.
    const S availability = growth_flux / (growth_flux + maintenance_flux_);
    sapwood_marginal_ = availability *
      (dP_dpsi - net_mass_production_dt_ * dmass_sap_dA *
         darea_leaf_dmass_live_) / maintenance_flux_;
  }
  return sapwood_marginal_;
}

// Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <typename S>
S TF24_Model<S>::respiration(const S& mass_leaf, const S& mass_sapwood,
                             const S& mass_bark, const S& mass_root,
                             const S& mass_coarse_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root) +
         respiration_coarse_root(mass_coarse_root);
}

template <typename S>
S TF24_Model<S>::respiration_leaf(const S& mass) const {
  return pars.r_l * mass;
}

template <typename S>
S TF24_Model<S>::respiration_bark(const S& mass) const {
  return pars.r_b * mass;
}

template <typename S>
S TF24_Model<S>::respiration_sapwood(const S& mass) const {
  return pars.r_s * mass;
}

template <typename S>
S TF24_Model<S>::respiration_root(const S& mass) const {
  return pars.r_r * mass;
}

template <typename S>
S TF24_Model<S>::respiration_coarse_root(const S& mass) const {
  return pars.r_cr * mass;
}

// Total turnover
template <typename S>
S TF24_Model<S>::turnover(const S& replacement_leaf,
                          const S& replacement_sapwood,
                          const S& mass_leaf, const S& mass_bark,
                          const S& mass_sapwood, const S& mass_root,
                          const S& mass_coarse_root) const {
   // At both fractions = 1 this is bit-identical to charging the whole flux --
   // 1.0 * x is exactly x, and the terms are summed in the same order.
   return replacement_leaf * turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          replacement_sapwood * turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root) +
          turnover_coarse_root(mass_coarse_root);
}

template <typename S>
S TF24_Model<S>::turnover_leaf(const S& mass) const {
  return pars.k_l * mass;
}

template <typename S>
S TF24_Model<S>::turnover_bark(const S& mass) const {
  return pars.k_b * mass;
}

template <typename S>
S TF24_Model<S>::turnover_sapwood(const S& mass) const {
  return pars.k_s * mass;
}

template <typename S>
S TF24_Model<S>::turnover_root(const S& mass) const {
  return pars.k_r * mass;
}

template <typename S>
S TF24_Model<S>::turnover_coarse_root(const S& mass) const {
  return pars.k_cr * mass;
}

// Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <typename S>
S TF24_Model<S>::net_mass_production_dt_A(const S& assimilation, const S& respiration,
                                const S& turnover) const {
  return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
}

// Fraction of production allocated to reproduction
template <typename S>
S TF24_Model<S>::fraction_allocation_reproduction(const S& height) const {
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <typename S>
S TF24_Model<S>::fraction_allocation_growth(const S& height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// Rate of offspring production
template <typename S>
S TF24_Model<S>::fecundity_dt(const S& net_mass_production_dt,
                               const S& fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <typename S>
S TF24_Model<S>::mortality_dt(const S& relative_reserves,
                              const S& cumulative_mortality) const {

  // Growth-dependent mortality is driven by relative NSC reserves
  // r = S/S_max (in [0,1]) rather than instantaneous productivity, so the rate
  // is bounded (see mortality_storage_dependent_dt): death is gradual as
  // reserves deplete rather than spiking under a carbon deficit.
  if (mortality_hazard_live(cumulative_mortality)) {
    return
      mortality_growth_independent_dt() +
      mortality_storage_dependent_dt(relative_reserves);
 } else {
    // Mortality probability is 1, so the rate calculations have nothing left to
    // describe and the state does not move again.
    return 0.0;
  }
}

template <typename S>
S TF24_Model<S>::mortality_growth_independent_dt() const {
  return pars.d_I;
}

// Storage-dependent growth mortality (#517), following Stefaniak et al. 2026
// (Eq 6). Bounded in [a_dG1*exp(-a_dG2), a_dG1] for r in [0,1]: full reserves
// (r=1) give near-zero excess mortality; empty reserves (r=0) give the finite
// maximum a_dG1. This boundedness is what removes the #550 ODE overflow.
template <typename S>
S TF24_Model<S>::mortality_storage_dependent_dt(const S& relative_reserves) const {
  return pars.a_dG1 * exp(-pars.a_dG2 * relative_reserves);
}

// NSC storage capacity: scales with sapwood mass (per Daniel, #517). mass_sapwood
// = area_sapwood(area_leaf, departure) * height * eta_c * rho.
template <typename S>
S TF24_Model<S>::storage_capacity(const S& area_leaf_, const S& height,
                                     const S& sapwood_departure) const {
  return pars.a_st1 *
    mass_sapwood(area_sapwood(area_leaf_, sapwood_departure), height);
}

// Initialise the storage state for a newly germinated individual at a_st3 fraction of
// its capacity (Stefaniak et al. 2026, Eq 8), so seedlings are born with
// reserves rather than starting empty (which would kill them immediately).
template <typename S>
void TF24_Model<S>::set_initial_states(const TF24_Environment<S>& environment,
                                       Internals<S>& vars) {
  // The seed's height is written here, not inherited from the height an
  // Individual was constructed with. Construction runs in plain arithmetic and an
  // active strategy receives its results rather than deriving them, so a height
  // that arrives that way carries no trait derivative and every rate the newborn
  // goes on to have inherits the loss. Declaring it by its own residual at this
  // point puts it where the parameters are already differentiable inputs, and the
  // value written is the same one either way.
  const SeedGeometry seed = seed_geometry();
  vars.set_state(HEIGHT_INDEX, seed.height);
  // A height written into vars leaves the slots it determines holding what the
  // constructor's plain height derived, which is the same number carrying no
  // derivative. Re-deriving them costs one allometry and is the only thing that
  // puts the seed's row into the leaf area every rate at birth size is scaled by.
  update_dependent_aux(HEIGHT_INDEX, vars);
  // A seed is on its allometry by definition, so both departures rest at zero.
  vars.set_state(state_idx_log_area_leaf_departure, S(0.0));
  vars.set_state(state_idx_log_area_sapwood_departure, S(0.0));
  vars.set_state(state_idx_storage,
                 pars.a_st3 * storage_capacity(seed.area_leaf, seed.height,
                                               S(0.0)));

  // Seed the tracked marginal balance from the plant's own first evaluation.
  // Left at zero it would start every seedling exactly at the gate's centre,
  // withholding half its replacement from birth; and a running mean has to
  // start somewhere its own history would have put it.
  //
  // Only where the gate is open. At a_pl0 = 0 the tracked value reaches no
  // rate, so the leaf solve this costs per birth would buy nothing -- and the
  // seed's carbon at birth is read from the first rate evaluation's aux anyway.
  S tracked = 0.0;
  if (odelia::util::to_passive(pars.a_pl0) > 0.0) {
    net_mass_production_dt(environment, vars);
    tracked = marginal_leaf_return_;
  }
  vars.set_state(state_idx_leaf_marginal_tracked, tracked);
}

// Survival of seedlings during establishment
template <typename S>
S TF24_Model<S>::establishment_probability(const TF24_Environment<S>& environment) {
  const SeedGeometry seed = seed_geometry();
  return establishment_probability(
    environment,
    net_mass_production_dt(environment, seed.height, seed.area_leaf,
                           1.0 / seed.height,
                           // A germinating seed is on its allometry by
                           // definition -- height_0 is solved from seed mass,
                           // and area_leaf_0 from height_0 -- so both of its
                           // departures are the resting ones. It has no history
                           // yet, so the tracked balance starts at break-even.
                           S(0.0), S(0.0), S(0.0)));
}

// Both forms above end here. The carbon is birth-size carbon either way, whatever
// height the caller's plant happens to be at.
template <typename S>
S TF24_Model<S>::establishment_probability(const TF24_Environment<S>& environment,
                                               const S& net_mass_production_dt_) {

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * seed_geometry().area_leaf / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

// Ratio-first hot-path overload (see header): receives the cached
// competition_effect (= area_leaf(height)) and height_inverse (= 1/height), so the
// per-call area_leaf() evaluation and z/height division are hoisted out of the
// inner competition loop.
template <typename S>
S TF24_Model<S>::compute_competition(const S& z, const S& area_leaf_,
                                          const S& height_inverse) const {
  return pars.k_I * area_leaf_ * canopy_shape.leaf_area_above(z * height_inverse);
}

}

#endif
