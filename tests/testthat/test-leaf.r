source("tests/leaf_functions.R")
context("SCM-general")


test_that("Basic functions", {
  devtools::load_all(".")
  library(tidyverse)

  #first set physiological parameters
  
  vcmax = 100 #maximum carboxylation rate, defined by leaf nitrogen (umol m^-2 s^-1) 
  p_50 = 1.731347 #stem water potential at 50% loss of conductivity
  
  c = 2.04 #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
  b = calc_vul_b(p_50 = p_50, c = c) #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
  psi_crit = calc_psi_crit(b, c) #stem water potential at which conductance is 95%
  huber_value = 0.000157 #huber value (m^2 sapwood area m^-2 leaf area)
  K_s = 2 #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
  h = 5 #height or path length (m)
  
  l <- Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
  
  #without setting physiology, PPFD_, k_l_max_, and psi_soil_ should all be NA
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$k_l_max_))
  expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$lambda_analytical_))
  
  #now set physiology, PPFD_, k_l_max_, and psi_soil_, atm_vpd_ should be not NA
  
  PPFD = 900
  k_l_max = calc_k_l_max(K_s, huber_value, h)
  psi_soil = 2
  atm_vpd = 2
  ca = 40
  
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  
  expect_equal(l$PPFD_, PPFD)
  expect_equal(l$k_l_max_, k_l_max)
  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$ca_, ca)
  
  #generating a new leaf object should wipe the previously stored values
  
  l <- Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$k_l_max_))
  expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$electron_transport_))
  
  
  #set physiology again for testing 
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, k_l_max = k_l_max, atm_vpd = atm_vpd, ca =ca)
  
  
  psi <- 1 #nominated value for water potential for testing vulnerability curve equations only (-MPa)
  
  #test conducitvity vulnerability, should be proportion value. 
  expect_equal(l$proportion_of_conductivity(psi), calc_k_l(psi, k_l_max, b, c) / k_l_max)
  
  
  #test calcuation of transpiration stream based on water potential of stem (-MPa)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is greater than psi_soil
  psi_stem <- psi_soil+1 #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem) > 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is less than psi_soil, creates negative value. Ordinarily an undesirable property which is typically bannen (stem assumed to have minimum water potential at psi_soil)
  psi_stem <- psi_soil-1 #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem) < 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is equal to psi_soil
  psi_stem <- psi_soil #stem water potential (-MPa)
  expect_true(l$transpiration(psi_stem) == 0)
  
  #for situations where psi_stem exceeds psi_crit + tolerance
  expect_error(l$transpiration(psi_crit+psi_crit*0.1), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #for situations where psi_soil exceeds psi_crit + tolerance
  
  psi_soil = psi_crit + psi_crit*0.1
  l$set_physiology(PPFD, psi_soil = psi_soil, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  psi_stem = psi_soil 
  
  expect_error(l$transpiration(psi_stem), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #test that fast E supply calculation is closely approximating full integration
  psi_soil = 0
  l$set_physiology(PPFD, psi_soil = psi_soil, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  psi_stem = psi_soil + 3
  
  expect_equal(l$transpiration(psi_stem), l$calc_E_supply_full_integration(psi_stem))
  
  #test that conversion between psi and E works properly
  
  expect_equal(l$transpiration_to_psi_stem(l$transpiration(psi_stem)), psi_stem)

  expect_equal(l$stom_cond_CO2(psi_stem), calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b), tolerance = 1e-5)

  c_i = 30 #intra-cellular carbon dioxide parital pressure (Pa)

  
  expect_equal(l$assim_rubisco_limited(c_i), calc_A_c(c_i = c_i, vcmax = vcmax))
  
  expect_equal(l$assim_electron_limited(c_i), calc_A_j(c_i = c_i, PPFD = PPFD, vcmax = vcmax))

  expect_equal(l$assim_colimited(c_i), calc_A_lim(c_i, vcmax, PPFD))
  
  expect_equal(l$assim_colimited_analytical(c_i), calc_A_lim_one_line(c_i, vcmax = vcmax, PPFD = PPFD))
 
  #test a function which retrieves various leaf-level states and rates from a given psi_stem value
  #for situations where psi stem is lower than psi soil

  psi_soil = 2
  l <- Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_soil - 1)
  
  #assimilation becomes 0
  expect_equal(l$assim_colimited_, 0)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_, 0)
  #transpiration becomes 0
  expect_equal(l$transpiration_, 0)
  
  #for situations where psi stem is same as psi soil
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_soil)
  
  #assimilation becomes 0
  expect_equal(l$assim_colimited_, 0)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_, 0)
  #transpiration becomes 0
  expect_equal(l$transpiration_, 0)
  
  #for situations where psi stem is greater than psi soil
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_soil + 1)
  
  #assimilation becomes 0
  expect_equal(l$assim_colimited_ >0, TRUE)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_ >0, TRUE)
  #transpiration becomes 0
  expect_equal(l$transpiration_ >0, TRUE)
  
  #calculate the hydraulic cost usign the sperry method, should be 0 when psi_soil is equivalent to psi_stem
  expect_equal(l$hydraulic_cost_Sperry(psi_soil) == 0, TRUE)
  
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_crit)
  
  expect_equal(l$assim_colimited_, calc_ben_gross_one_line(psi_stem = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c, PPFD = PPFD, vcmax= vcmax, ca = ca), tolerance = 1e-03)
  
  expect_equal(l$hydraulic_cost_Sperry(psi_crit), calc_hydraulic_cost(psi_soil, psi_crit, k_l_max, b, c))
  
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_stem)

  expect_equal(l$assim_colimited_, calc_ben_gross_one_line(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c, PPFD = PPFD, vcmax= vcmax, ca = ca), tolerance = 1e-04)
  
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_crit)
  
  #this value is constant following a call to set_physiology
  expect_equal(l$lambda_analytical_, calc_ben_gross_one_line(psi_stem = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c, PPFD = PPFD, vcmax= vcmax, ca = ca) /
    calc_hydraulic_cost(psi_soil, psi_crit, k_l_max, b, c), tolerance = 1)
  
  #under almost all scenarios, max ci (i.e when psi stem is set to psi crit) should be less than ca
  
  #plain version use a uniroot solving method to find ci
  l$set_leaf_states_rates_from_psi_stem(psi_crit)
  expect_equal(l$ci< ca, TRUE)
  
  #one line method uses a more efficient, analytical based solution by using a calibrated constant
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_crit)
  expect_equal(l$ci< ca, TRUE)
  
  l$set_leaf_states_rates_from_psi_stem(psi_crit)
  ci_numerical = l$ci
  
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_crit)
  ci_analytical = l$ci
  
  #test whether max cis are roughly equivalent between numerical and analytical
  expect_equal(ci_numerical, ci_analytical, tolerance = 0.5)
  
  #test whether finding the maximum ci (i.e at psi_crit) is roughly equivalent
  expect_equal(ci_analytical - find_max_ci_one_line(psi_crit = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, vcmax = vcmax, PPFD = PPFD), 0, tolerance = 1e-03)
  
  #test whether conversion between E and psi is equivalent between R and C++
  
  E = 0.0001109062
  
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_crit)
  c_i = l$ci  
  
  benefit_ = l$assim_colimited_analytical(c_i);
  g_c_ci = (benefit_ * umol_per_mol_2_mol_per_mol * atm_kpa * kPa_2_Pa)/(ca - c_i); 
  E_ci = g_c_ci * 1.6 * atm_vpd / kg_2_mol_h20 / atm_kpa;
  psi_stem = l$transpiration_to_psi_stem(E_ci)
  
  #conversion back and forth is not perfect
  expect_equal(psi_stem, psi_crit, tolerance = 1e-06)
  
  #let's start testing profit functions
  
  psi_soil = 0
  l <- Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
  l$set_physiology(PPFD = PPFD, psi_soil = psi_soil, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  
  expect_equal(l$profit_psi_stem_Sperry_analytical(2), calc_profit_Sperry_psi_stem_one_line(psi_stem = 2, k_l_max = k_l_max, b = b, c = c), tolerance = 0.01)
  
  #test gss optimiser
  
  #first, need to find maximum ci for the c++ implementation, this is done autmoatically in R
  l$set_leaf_states_rates_from_psi_stem_analytical(psi_crit)  
  l$optimise_ci_Sperry_analytical(l$ci)  
  
  #compare values - tolerance pretty high because there is an integration method as well as a uni-root solver causing diffrences between R and c++
  expect_equal(l$profit_, opt_ci_gss(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b=b, c=c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca =ca), tolerance = 0.1)

  #test some scenarios just using the c++ implementation - if they deviate strongly we have done something to change the implementation significantly which will warn us if we don't think we've done anything
  
  #first off- what happens when we moving psi_soil around
  l$set_physiology(PPFD = 900, psi_soil = psi_crit + 1, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  l$optimise_psi_stem_Sperry_Newton_analytical(NA)
  
  expect_equal(l$profit_, 0)
  expect_equal(l$opt_psi_stem_, psi_crit+1)
  expect_equal(l$transpiration_, 0)
  
  l$set_physiology(PPFD = 900, psi_soil = 0, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  l$optimise_psi_stem_Sperry_Newton_analytical(NA)
  expect_equal(l$profit_, 14.14675, tolerance = 1e-6)
  expect_equal(l$opt_psi_stem_, 0.8451968, tolerance = 1e-6)
  
  l <- Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
  
  l$set_physiology(PPFD = 0, psi_soil = 0, k_l_max = k_l_max, atm_vpd = atm_vpd, ca = ca)
  l$optimise_psi_stem_Sperry_Newton_analytical(NA)
  
  
  #when no profit can be made at any point at the three points used for newtons method, assume profit = 0, transpiration= 0, and opt_psi_stem is set to psi_soil, representing the plant completely closing stomata. 
  expect_equal(l$profit_, 0, tolerance = 1e-6)
  expect_equal(l$opt_psi_stem_, l$psi_soil_, tolerance = 1e-6)
})
  
  
  
