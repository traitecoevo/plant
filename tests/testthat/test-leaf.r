context("Leaf-phys")

test_that("Basic functions", {
  #first set physiological parameters
  
  # TF24_strategy <- TF24_Strategy()
  vcmax_25 = 100 #maximum carboxylation rate (umol m^-2 s^-1) 
  jmax_25 = vcmax_25*167 #maximum electron transport rate (umol m^-2 s^-1) 
  p_50 = 2 #stem water potential at 50% loss of conductivity
  c = 2.04 #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
  b = 3 #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
  psi_crit = 5 #stem water potential at which conductance is 95%
  theta = 0.000157 #huber value (m^2 sapwood area m^-2 leaf area)
  K_s = 1 #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
  h = 5 #height or path length (m)
  beta2 = 1
  hk_s = 75
  curv_fact_elec_trans = 0.7
  a = 0.3
  curv_fact_colim = 0.99
  g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8
  vulnerability_curve_ncontrol = 100
  ci_abs_tol = 1e-6
  ci_niter = 1000
  beta_R_H = 3.4e3
  beta_R_V = 9.4e4
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  
  #without setting physiology, PPFD_, k_l_max_, and psi_soil_ should all be NA
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$leaf_specific_conductance_max_))
  # expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  expect_true(is.na(l$atm_o2_kpa_))
  expect_true(is.na(l$leaf_temp_))
  expect_true(is.na(l$ci_))
  expect_true(is.na(l$stom_cond_CO2_))
  expect_true(is.na(l$assim_colimited_))
  expect_true(is.na(l$transpiration_))
  expect_true(is.na(l$profit_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  expect_true(is.na(l$hydraulic_cost_))
  expect_true(is.na(l$electron_transport_))
  expect_true(is.na(l$gamma_))
  expect_true(is.na(l$ko_))
  expect_true(is.na(l$kc_))
  expect_true(is.na(l$km_))
  expect_true(is.na(l$R_d_))
  expect_true(is.na(l$rho_))
  expect_true(is.na(l$vcmax_))
  expect_true(is.na(l$jmax_))
  expect_true(is.na(l$a_bio_))
  expect_true(is.na(l$root_collar_psi_))
  expect_true(is.na(l$opt_psi_stem_))
  expect_true(is.na(l$opt_ci_))
  expect_true(is.na(l$E_up_))

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
  
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = 1, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  
  
  expect_equal(l$c_r_V_total_, 1*(1.0/3.0))
  expect_equal(l$c_r_H_total_, 1*(2.0/3.0))
  expect_equal(l$PPFD_, PPFD)
  expect_equal(l$leaf_specific_conductance_max_, leaf_specific_conductance_max)
  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$ca_, ca)
  expect_equal(l$atm_vpd_, atm_vpd)
  expect_equal(l$atm_o2_kpa_, atm_o2_kpa_)
  expect_equal(l$leaf_temp_, leaf_temp_)
  expect_equal(l$atm_kpa_, atm_kpa_)
  expect_equal(l$rho_, 608)
  expect_equal(l$a_bio_, 0.0245)
  expect_equal(l$sapwood_volume_per_leaf_area_, sapwood_volume_per_leaf_area)
  expect_equal(l$soil_depth_, 1)
  expect_equal(l$soil_number_of_depths_, length(psi_soil))
  expect_equal(l$root_mass_, 1)



  # test that the inputs to set_physiology which take multiple values are working correctly

  psi_soil = c(1,2)
  soil_depth = c(0.5, 1)

  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$soil_number_of_depths_, length(psi_soil))
  expect_equal(l$soil_depth_, soil_depth)

  #generating a new leaf object should wipe the previously stored values
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$leaf_specific_conductance_max_))
  expect_true(is.na(l$soil_depth_))
  expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  expect_true(is.na(l$atm_o2_kpa_))
  expect_true(is.na(l$leaf_temp_))
  expect_true(is.na(l$ci_))
  expect_true(is.na(l$stom_cond_CO2_))
  expect_true(is.na(l$assim_colimited_))
  expect_true(is.na(l$transpiration_))
  expect_true(is.na(l$profit_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  expect_true(is.na(l$hydraulic_cost_))
  expect_true(is.na(l$electron_transport_))
  expect_true(is.na(l$gamma_))
  expect_true(is.na(l$ko_))
  expect_true(is.na(l$kc_))
  expect_true(is.na(l$km_))
  expect_true(is.na(l$R_d_))
  expect_true(is.na(l$rho_))
  expect_true(is.na(l$vcmax_))
  expect_true(is.na(l$jmax_))
  expect_true(is.na(l$a_bio_))
  expect_true(is.na(l$root_collar_psi_))
  expect_true(is.na(l$opt_psi_stem_))
  expect_true(is.na(l$opt_ci_))
  expect_true(is.na(l$E_up_))
  expect_true(is.na(l$soil_number_of_depths_))

  
  #set physiology again for testing 
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  psi <- 1 #nominated value for water potential for testing vulnerability curve equations only (-MPa)
  
  #test conducitvity vulnerability, should be proportion value. 
  expect_equal(l$proportion_of_conductivity(psi), 0.8991241827)
  
  #test calcuation of transpiration stream based on water potential of stem (-MPa)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is greater than psi_soil
  psi_stem <- psi_soil+1 #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem[1], psi_soil[1]) > 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is less than psi_soil, creates negative value. Ordinarily an undesirable property which is typically banned (stem assumed to have minimum water potential at psi_soil)
  psi_stem <- psi_soil-1 #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem[1], psi_soil[1]) < 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is equal to psi_soil
  psi_stem <- psi_soil #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem[1], psi_soil[1]) == 0)
  
  upper_bound_int <- 3*((log(1/1e-5))^(1/2.04))
  #for situations where psi_stem exceeds tolerance of integrator
  expect_error(l$transpiration(upper_bound_int, psi_stem[1]), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #for situations where psi_soil exceeds psi_crit + tolerance
  
  psi_soil = upper_bound_int
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  psi_stem = psi_soil 
  
  expect_error(l$transpiration(psi_stem[1], psi_soil[1]), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #test that fast E supply calculation is closely approximating full integration
  psi_soil = 0
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  psi_stem = psi_soil + 3
  
  expect_equal(l$transpiration(psi_stem[1], psi_soil[1]), l$transpiration_full_integration(psi_stem[1], psi_soil[1]))
  
  #test that conversion between psi and E works properly
  
  expect_equal(l$transpiration_to_psi_stem(l$transpiration(psi_stem[1], psi_soil[1]), psi_soil[1]), psi_stem[1], tolerance = 1e-5)

  c_i = 30 #intra-cellular carbon dioxide parital pressure (Pa)

  #test a function which retrieves various leaf-level states and rates from a given psi_stem value
  #for situations where psi stem is lower than psi soil

  psi_soil = 2
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  
  #note that this scenario should not occur in model anyway
  l$set_leaf_states_rates_from_psi_stem(psi_soil - 1, psi_soil)

  #assimilation becomes dark respiration
  expect_equal(l$assim_colimited_, -l$R_d_)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_, 0)
  #transpiration becomes 0
  expect_equal(l$transpiration_, 0)
  
  #costs 0 when psi_stem == psi_soil == 0
  expect_equal(l$hydraulic_cost_Sperry(psi_stem = psi_soil, psi_upstream = psi_soil) == 0, TRUE)
  #costs positive even when transpiration stream is 0 in hydraulic cost tf
  expect_equal(l$hydraulic_cost_TF(psi_soil) > 0, TRUE)
  
  #when psi stem is greater than psi soil
  psi_soil = 2
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  
  l$set_leaf_states_rates_from_psi_stem(psi_soil + 1, psi_soil)
  

  #assimilation becomes greater than 0 
  expect_equal(l$assim_colimited_ >0, TRUE)
  #stomatal conductance becomes greater than 0 
  expect_equal(l$stom_cond_CO2_ >0, TRUE)
  #transpiration becomes greater than 0 
  expect_equal(l$transpiration_ >0, TRUE)
  
  #when psi stem is equal to psi soil
  psi_soil = 2
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  
  l$set_leaf_states_rates_from_psi_stem(psi_soil, psi_soil)
  
  
  #calculate the hydraulic cost usign the sperry method, should be 0 when psi_soil is equivalent to psi_stem
  expect_equal(l$hydraulic_cost_Sperry(psi_soil, psi_soil) == 0, TRUE)
  #calculate hydraulic cost using TF method, should be greater than 0 when psi_soil is greater than 0
  expect_equal(l$hydraulic_cost_TF(psi_soil) > 0, TRUE)
  
  expect_equal(l$hydraulic_cost_Sperry(psi_soil + 1, psi_soil) > 0, TRUE)
  expect_equal(l$hydraulic_cost_TF(psi_soil + 1) > 0, TRUE)
  
  #ensure that hydraulic_cost returns 0 cost at 0 psi_soil/psi_stem
  psi_soil = 0
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  
  l$set_leaf_states_rates_from_psi_stem(0, 0)
  expect_equal(l$hydraulic_cost_TF(psi_soil) == 0, TRUE)
  expect_equal(l$hydraulic_cost_Sperry(psi_soil, psi_soil) == 0, TRUE)
  
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
  l$set_leaf_states_rates_from_psi_stem(psi_crit, psi_soil)
  expect_equal(l$ci_< ca, TRUE)
  
  #test whether conversion between E and psi is equivalent between R and C++
  
  psi_soil = 0.5
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  
  l$set_leaf_states_rates_from_psi_stem(psi_crit, psi_soil)
  c_i = l$ci_  
  benefit_ = l$assim_colimited_
  
  kg_to_mol_h2o = 55.4939
  umol_to_mol = 1e-6
  kPa_to_Pa = 1e3
  
  g_c_ci = ((benefit_)* umol_to_mol * l$atm_kpa_ * kPa_to_Pa)/(l$ca_ - l$ci_); 
  
  E_ci = g_c_ci * 1.67 * l$atm_vpd_ / kg_to_mol_h2o / l$atm_kpa_;
  # need to invert psi_stem to get the same value as psi_crit, which is the value we set for psi_stem in this test
  psi_stem = l$transpiration_to_psi_stem(E_ci, -psi_soil)
  
  #conversion back and forth is not perfect
  expect_equal(psi_stem, psi_crit, tolerance = 1e-05)
  
  #let's start testing profit functions
  
  #first off- what happens when we moving psi_soil around
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
    
    psi_soil = 0
    soil_depth = 0.5

  #first off- what happens when we moving psi_soil around
  # start with one soil layer. For the TF method, it will fail when there is more than one layer
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_crit+1, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
    l$optimise_psi_stem_TF()
  
  expect_equal(l$transpiration_, 0)
  expect_equal(l$opt_psi_stem_, psi_crit+1)
  expect_equal(l$stom_cond_CO2_, 0)
  expect_equal(l$ci_, l$gamma_*0.1013)
  expect_true(l$profit_ < 0)
  
  # check a more standard case, wet soil

  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  expect_true(l$profit_ > 0)
  expect_true(l$opt_psi_stem_ > 0)
  expect_true(l$stom_cond_CO2_ > 0)
  expect_true(l$transpiration_ > 0)
  expect_true(l$hydraulic_cost_ > 0)
  profit1 = l$profit_

  expect_true(is.na(l$profit_psi_stem_TF(NA, psi_soil)))

  #confirm different nubmer of alyers fails
  psi_soil = c(0, 1)
  soil_depth = c(0.5, 1)

  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  expect_error(l$optimise_psi_stem_TF(), "psi soil must have only one value to use non-root-based profit optimisation methods")

  #test various responses to environmental gradients to check that behaviour is being conserved
  
psi_soil = 1

  #light
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  high_light <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 100, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$optimise_psi_stem_TF()
  
  low_light <- l$profit_
  
  expect_true(high_light > low_light)
  
  #soil moist
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  high_moist <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 2, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  low_moist <- l$profit_
  
  expect_true(high_moist > low_moist)
  
  
  #vpd'
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  low_vpd <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 2, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  high_vpd <- l$profit_
  
  expect_true(high_vpd < low_vpd)
  
  #vcmax_25
  
  l <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  low_vcmax <- l$profit_
  
  l <- Leaf(vcmax_25 = 150, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$optimise_psi_stem_TF()
  
  high_vcmax <- l$profit_
  
  expect_true(high_vcmax > low_vcmax)
  
  #test effect of leaf temperature
  
  l_low_temp <-Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
                            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l_low_temp$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = 20, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l_ref_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
                             beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                             ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l_ref_temp$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = 25, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l_high_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
                              beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                              GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                              ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l_high_temp$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = 30, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  expect_true(l_low_temp$gamma_ < l_ref_temp$gamma_ &  l_ref_temp$gamma_ <  l_high_temp$gamma_)
  
  l_high_temp <-Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
                             beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                             ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  l_high_temp$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = 40, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  expect_equal(round(l_high_temp$ko_, 1), round(562314.4,1))
  expect_equal(round(l_high_temp$kc_, 1), round(1879.0751,1))
  expect_equal(round(l_high_temp$gamma_, 1), round(88.800391,1))
  expect_equal(round(l_high_temp$vcmax_, 1), round(35.2,1))

  # test out the root component of the leaf unit



  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24, beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  
  # test that you get error when soil_depth and psi_soil have different number of layers

  soil_depth = c(0.5,1)
  psi_soil = c(0.5)

  expect_error(l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_), "soil_depth and psi_soil must have the same number of elements")

  soil_depth = c(0.5,1)
  psi_soil = c(0.5,1)

l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

# root assignment should be equal to number of soil layers  
expect_equal(length(l$c_r_H_), length(soil_depth))
expect_equal(length(l$c_r_V_), length(soil_depth))

l$E_from_Soil_to_Root_Collar()

devtools::load_all("/Users/z3533545/Documents/GitHub/plant/")

  vcmax_25 = 100 #maximum carboxylation rate (umol m^-2 s^-1) 
  jmax_25 = vcmax_25*167 #maximum electron transport rate (umol m^-2 s^-1) 
  p_50 = 2 #stem water potential at 50% loss of conductivity
  c = 2.04 #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
  b = 3 #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
  psi_crit = 5 #stem water potential at which conductance is 95%
  theta = 0.000157 #huber value (m^2 sapwood area m^-2 leaf area)
  K_s = 1 #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
  h = 5 #height or path length (m)
  beta2 = 1
  hk_s = 75
  curv_fact_elec_trans = 0.7
  a = 0.3
  curv_fact_colim = 0.99
  g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8
  vulnerability_curve_ncontrol = 100
  ci_abs_tol = 1e-6
  ci_niter = 1000
  beta_R_H = 3.4e3
  beta_R_V = 9.4e4
  PPFD = 900
  sapwood_volume_per_leaf_area = theta*h
  leaf_specific_conductance_max = K_s*theta/h
  psi_soil = 2
  atm_vpd = 2
  ca = 40
  atm_o2_kpa_ = 21
  leaf_temp_ = 25
  atm_kpa_ = 101.3

  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, b = b, psi_crit = psi_crit, 
            beta2= beta2, hk_s = hk_s, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24/10, beta_R_H = beta_R_H, beta_R_V = beta_R_V)

  soil_depth = c(0.5,1)
  psi_soil = c(0.5,1)

l$set_physiology(root_mass = 1, rho = 608, a_bio = 0.0245, PPFD = 1000, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, sapwood_volume_per_leaf_area = sapwood_volume_per_leaf_area, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)


l$find_root_collar_psi()
l$root_collar_psi_
bound_a = l$find_root_psi(-0.5, -psi_soil, 1)
bound_b = l$find_root_psi(-0.5, -psi_soil, 0)

# l$find_psi_stem_from_psi_root((root_crit), -psi_soil)

profit <- c()
benefit <- c()
cost <- c()
transpiration_ <- c()
psi_stem <- c()
psi_root <- c()
soil_consumption <- list()
for(i in 1:length(sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100))))){
  psi_stem_c = l$find_psi_stem_from_psi_root(sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100)))[i], -psi_soil)
  l$root_collar_psi_ = -l$root_collar_psi_
  psi_root[i] = l$root_collar_psi_
  profit[i] <- l$profit_psi_stem_TF(psi_stem_c, -sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100)))[i])
  benefit[i] <- l$assim_colimited_
  cost[i] <- l$hydraulic_cost_
  transpiration_[i] <- l$transpiration_
  psi_stem[i] <- psi_stem_c
  soil_consumption[[i]] <- l$soil_consumption_
}

tibble(profit, psi_stem, psi_root) %>% View()

cols <- c("psi_stem" = "red", "root_collar_psi" = "blue", "soil_layer1" = "darkgreen", "soil_layer2" = "orange")

tibble(psi_stem, root_collar_psi2 =  sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100))), root_collar_psi = sort(c(1.102122, seq(-bound_a, -bound_b, length.out = 100)), decreasing = TRUE), soil_layer1 = 0.5, soil_layer2 = 1) %>%
  pivot_longer(cols = -root_collar_psi2) %>%
  ggplot(aes(x=  root_collar_psi2, y =-value)) +
  geom_line(aes(colour = name, group = name), size =2) +
  # geom_hline(aes(yintercept = root_collar_psi2[1])) +
  geom_vline(aes(xintercept = seq(bound_a, bound_b, length.out = 100)[profit == max(profit)])) +
  theme_bw() +
  xlab(expression(paste(psi[rc]," (MPa)"))) +
  ylab(expression(paste(psi[x]," (MPa)"))) +
  theme(text = element_text(size = 20)) +
  scale_x_reverse() + 
  scale_colour_manual(values = cols)-> a
cols <- c("total_transpiration" = "purple", "first_layer" = "darkgreen", "second_layer" = "orange")

tibble(root_collar_psi = sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100))), profit, benefit, cost) %>%
  pivot_longer(-root_collar_psi) %>%
  ggplot(aes(x= root_collar_psi, y =value)) +
  geom_line(aes(colour = name, group = name), size =2) +
  geom_vline(aes(xintercept = sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100)))[profit == max(profit)])) +
  
  theme_bw() +
  xlab(expression(paste(psi[rc]," (MPa)"))) +
  ylab("Carbon budget") +
  scale_x_reverse() + 
  theme(text = element_text(size = 20)) -> b

soil_consumption %>%
  map_dfr(., ~tibble(first_layer = .x[1], second_layer = .x[2])) %>%
  mutate(root_collar_psi = sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100)))) %>%
  mutate(total_transpiration = first_layer + second_layer) %>%
  pivot_longer(cols = -root_collar_psi) %>%
  ggplot(aes(x = root_collar_psi, y= value)) +
  geom_line(aes(colour = name, group = name), size =2) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = sort(c(-1.102122, seq(bound_a, bound_b, length.out = 100)))[profit == max(profit)])) +
  
  theme_bw() +
  xlab(expression(paste(psi[rc]," (MPa)"))) +
  ylab("Transpiration") +
  scale_x_reverse() +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values = cols) -> c_plot
library(patchwork)
a/c_plot/b
ggsave("test_leaf_root_collar.png", width = 10, height = 15)
l$find_root_collar_psi()
l$find_








max(l$psi_soil_)

E <- c()
for(i in 1:length(seq(-0.1,-5,length.out = 100))){
E[i]<-l$E_from_Soil_to_Root_Collar(P_x_r = seq(-0.1,-5,length.out = 100)[i], P_soil = c(-0.0597947, -0.1221), dz = 0.1, LA = 1)
}
E




l$E_from_Soil_to_Root_Collar(P_x_r = -5, P_soil = c(-0.0597947, -0.1221), dz = 0.1, z_soil_mid = c(0.05, 0.15), LA = 1, c_r_H = c(20,30), c_r_V = c(20,30), beta_R_H = 3.4e3, beta_R_V = 9.4e4)

  })
