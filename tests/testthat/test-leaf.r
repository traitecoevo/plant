context("SCM-general")

test_that("Basic functions", {
  #first set physiological parameters
  
  FF16w_strategy <- FF16w_Strategy()
  
  vcmax_25 = FF16w_strategy$vcmax_25 #maximum carboxylation rate (umol m^-2 s^-1) 
  jmax_25 = FF16w_strategy$jmax_25 #maximum electron transport rate (umol m^-2 s^-1) 
  p_50 = FF16w_strategy$p_50 #stem water potential at 50% loss of conductivity
  c = FF16w_strategy$c #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
  b = FF16w_strategy$b #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
  psi_crit = FF16w_strategy$psi_crit #stem water potential at which conductance is 95%
  theta = FF16w_strategy$theta #huber value (m^2 sapwood area m^-2 leaf area)
  K_s = FF16w_strategy$K_s #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
  h = 5 #height or path length (m)
  beta1 = FF16w_strategy$beta1
  beta2 = FF16w_strategy$beta2
  hk_s = FF16w_strategy$hk_s
  curv_fact_elec_trans = FF16w_strategy$curv_fact_elec_trans
  a = FF16w_strategy$a
  curv_fact_colim = FF16w_strategy$curv_fact_colim
  
  GSS_tol_abs = FF16w_strategy$control$GSS_tol_abs
  vulnerability_curve_ncontrol = FF16w_strategy$control$vulnerability_curve_ncontrol
  ci_abs_tol = FF16w_strategy$control$ci_abs_tol
  ci_niter = FF16w_strategy$control$ci_niter
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  
  #without setting physiology, PPFD_, k_l_max_, and psi_soil_ should all be NA
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$leaf_specific_conductance_max_))
  expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  expect_true(is.na(l$atm_o2_kpa_))
  expect_true(is.na(l$leaf_temp_))
  
  #now set physiology, PPFD_, k_l_max_, and psi_soil_, atm_vpd_ should be not NA
  
  PPFD = 900
  sapwood_volume_per_leaf_area = theta*h
  leaf_specific_conductance_max = K_s*theta/h
  psi_soil = 2
  atm_vpd = 2
  ca = 40
  atm_o2_kpa_ = 21
  leaf_temp_ = 25
  atm_kpa_ = 101.3
  
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  expect_equal(l$PPFD_, PPFD)
  expect_equal(l$leaf_specific_conductance_max_, leaf_specific_conductance_max)
  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$ca_, ca)
  expect_equal(l$atm_vpd_, atm_vpd)
  expect_equal(l$atm_o2_kpa_, atm_o2_kpa_)
  expect_equal(l$leaf_temp_, leaf_temp_)

  #generating a new leaf object should wipe the previously stored values
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = 1000, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$leaf_specific_conductance_max_))
  expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  expect_true(is.na(l$atm_o2_kpa_))
  expect_true(is.na(l$leaf_temp_))
  
  #set physiology again for testing 
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_,atm_kpa = atm_kpa_)
  
  psi <- 1 #nominated value for water potential for testing vulnerability curve equations only (-MPa)
  
  #test conducitvity vulnerability, should be proportion value. 
  expect_equal(l$proportion_of_conductivity(psi), 0.9742527)
  
  #test calcuation of transpiration stream based on water potential of stem (-MPa)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is greater than psi_soil
  psi_stem <- psi_soil+1 #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem) > 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is less than psi_soil, creates negative value. Ordinarily an undesirable property which is typically banned (stem assumed to have minimum water potential at psi_soil)
  psi_stem <- psi_soil-1 #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem) < 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is equal to psi_soil
  psi_stem <- psi_soil #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem) == 0)
  
  upper_bound_int <- FF16w_strategy$b*((log(1/1e-5))^(1/FF16w_strategy$c))
  #for situations where psi_stem exceeds tolerance of integrator
  expect_error(l$transpiration(upper_bound_int), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #for situations where psi_soil exceeds psi_crit + tolerance
  
  psi_soil = upper_bound_int
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  psi_stem = psi_soil 
  
  expect_error(l$transpiration(psi_stem), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #test that fast E supply calculation is closely approximating full integration
  psi_soil = 0
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  psi_stem = psi_soil + 3
  
  expect_equal(l$transpiration(psi_stem), l$transpiration_full_integration(psi_stem))
  
  #test that conversion between psi and E works properly
  
  expect_equal(l$transpiration_to_psi_stem(l$transpiration(psi_stem)), psi_stem)

  c_i = 30 #intra-cellular carbon dioxide parital pressure (Pa)

  #test a function which retrieves various leaf-level states and rates from a given psi_stem value
  #for situations where psi stem is lower than psi soil

  psi_soil = 2
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = 1000, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l$set_leaf_states_rates_from_psi_stem(psi_soil - 1)
  
  #assimilation becomes dark respiration
  expect_equal(l$assim_colimited_, -l$R_d_)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_, 0)
  #transpiration becomes 0
  expect_equal(l$transpiration_, 0)
  
  #costs 0 when psi_stem == psi_soil == 0
  expect_equal(l$hydraulic_cost_Sperry(psi_soil) == 0, TRUE)
  expect_equal(l$hydraulic_cost_TF(psi_soil) > 0, TRUE)
  
  #for situations where psi stem is greater than psi soil
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_soil + 1)
  
  #assimilation becomes greater than 0 
  expect_equal(l$assim_colimited_ >0, TRUE)
  #stomatal conductance becomes greater than 0 
  expect_equal(l$stom_cond_CO2_ >0, TRUE)
  #transpiration becomes greater than 0 
  expect_equal(l$transpiration_ >0, TRUE)
  
  #calculate the hydraulic cost usign the sperry method, should be 0 when psi_soil is equivalent to psi_stem
  expect_equal(l$hydraulic_cost_Sperry(psi_soil) == 0, TRUE)
  #calculate hydraulic cost using TF method, should be greater than 0 when psi_soil is greater than 0
  expect_equal(l$hydraulic_cost_TF(psi_soil) > 0, TRUE)
  
  expect_equal(l$hydraulic_cost_Sperry(psi_soil + 1) > 0, TRUE)
  expect_equal(l$hydraulic_cost_TF(psi_soil + 1) > 0, TRUE)
  
  #ensure that hydraulic_cost returns 0 cost at 0 psi_soil/psi_stem
  psi_soil = 0
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = 1000, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l$set_leaf_states_rates_from_psi_stem(0)
  expect_equal(l$hydraulic_cost_TF(psi_soil) == 0, TRUE)
  expect_equal(l$hydraulic_cost_Sperry(psi_soil) == 0, TRUE)
  
  #psi_soil == psi_stem means A == -R_d_
  expect_equal(l$assim_colimited_, -l$R_d_)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_, 0)
  #transpiration becomes 0
  expect_equal(l$transpiration_, 0)
  
  #test behaviours related to low Ci
  expect_equal(l$assim_rubisco_limited(l$gamma_*0.1013), 0)
  expect_equal(l$assim_electron_limited(l$gamma_*0.1013), 0)
  expect_equal(l$assim_colimited(l$gamma_*0.1013), -l$R_d_)
  

  #under almost all scenarios, max ci (i.e when psi stem is set to psi crit) should be less than ca
  
  #plain version use a uniroot solving method to find ci
  l$set_leaf_states_rates_from_psi_stem(psi_crit)
  expect_equal(l$ci_< ca, TRUE)
  
  #test whether conversion between E and psi is equivalent between R and C++
  
  psi_soil = 0.5
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = 1000, ci_abs_tol = 1e-05, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l$set_leaf_states_rates_from_psi_stem(psi_crit)
  c_i = l$ci_  
  benefit_ = l$assim_colimited_
  
  kg_to_mol_h2o = 55.4939
  umol_to_mol = 1e-6
  kPa_to_Pa = 1e3
  
  g_c_ci = ((benefit_ + l$R_d_)* umol_to_mol * l$atm_kpa_ * kPa_to_Pa)/(l$ca_ - l$ci_); 
  
  E_ci = g_c_ci * 1.67 * l$atm_vpd_ / kg_to_mol_h2o / l$atm_kpa_;
  psi_stem = l$transpiration_to_psi_stem(E_ci)
  
  #conversion back and forth is not perfect
  expect_equal(psi_stem, psi_crit, tolerance = 1e-05)
  
  #let's start testing profit functions
  
  #first off- what happens when we moving psi_soil around
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

  psi_soil = 0
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)  
  #first off- what happens when we moving psi_soil around
  l$set_physiology(PPFD = PPFD, psi_soil = psi_crit + 1, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  expect_equal(l$transpiration_, 0)
  expect_equal(l$opt_psi_stem_, psi_crit+1)
  expect_equal(l$stom_cond_CO2_, 0)
  expect_equal(l$ci_, l$gamma_*0.1013)
  expect_true(l$profit_ < 0)
  
  l$set_physiology(PPFD = PPFD, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  expect_true(l$profit_ > 0)
  expect_true(l$opt_psi_stem_ > 0)
  expect_true(l$stom_cond_CO2_ > 0)
  expect_true(l$transpiration_ > 0)
  expect_true(l$hydraulic_cost_ > 0)
  
  expect_true(is.na(l$profit_psi_stem_TF(NA)))

  #test various responses to environmental gradients to check that behaviour is being conserved
  
  #light
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  high_light <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 500, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  low_light <- l$profit_
  
  expect_true(high_light > low_light)
  
  #soil moist
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  high_moist <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 2, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  low_moist <- l$profit_
  
  expect_true(high_moist > low_moist)
  
  
  #vpd'
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  low_vpd <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 3, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  high_vpd <- l$profit_
  
  expect_true(high_vpd < low_vpd)
  
  #ca
  
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 20, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  low_ca_profit <- l$profit_
  low_ca_stom_cond_co2 <- l$stom_cond_CO2_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  high_ca_profit <- l$profit_
  high_ca_stom_cond_co2 <- l$stom_cond_CO2_
  
  expect_true(low_ca_profit < high_ca_profit)
  expect_true(low_ca_stom_cond_co2 > high_ca_stom_cond_co2)
  
  
  #vcmax_25
  
  l <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  low_vcmax <- l$profit_
  
  l <- Leaf(vcmax_25 = 150, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter)
  l$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  high_vcmax <- l$profit_
  
  expect_true(high_vcmax > low_vcmax)
  
  #test effect of leaf temperature
  
  l_low_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
                            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                            ci_niter = ci_niter)
  l_low_temp$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = 20, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l_ref_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
                            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                            n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                            ci_niter = ci_niter)
  l_ref_temp$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = 25, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l_high_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
                             beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                             n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                             ci_niter = ci_niter)
  l_high_temp$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = 30, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  expect_true(l_low_temp$gamma_ < l_ref_temp$gamma_ &  l_ref_temp$gamma_ <  l_high_temp$gamma_)
  
  l_high_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, beta1 = beta1, 
                             beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                             n, GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                             ci_niter = ci_niter)
  l_high_temp$set_physiology(PPFD = 1000, psi_soil = 0, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, rho = FF16w_strategy$rho, a_bio = FF16w_strategy$a_bio, leaf_temp = 40, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  expect_equal(round(l_high_temp$ko_, 1), round(562314.4,1))
  expect_equal(round(l_high_temp$kc_, 1), round(1879.0751,1))
  expect_equal(round(l_high_temp$gamma_, 1), round(88.800391,1))
  expect_equal(round(l_high_temp$vcmax_, 1), round(35.2,1))
  
  })
