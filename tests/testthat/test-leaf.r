source("tests/reference_leaf_updated.R")
context("SCM-general")


test_that("Basic functions", {
  
  #first set physiological parameters
  
  vcmax = 100 #maximum carboxylation rate, defined by leaf nitrogen (umol m^-2 s^-1) 
  p_50 = 1.731347 #stem water potential at 50% loss of conductivity
  b = 2.072101 #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
  c = 2.04 #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
  psi_crit = calc_psi_crit(b, c) #stem water potential at which conductance is 95%
  huber_value = 0.000157 #huber value (m^2 sapwood area m^-2 leaf area)
  K_s = 2 #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
  h = 5 #height or path length (m)
  
  l <- Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, beta=15000, beta_2 = 1, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
  
  #without setting physiology, PPFD_, k_l_max_, and psi_soil_ should all be NA
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$k_l_max_))
  expect_true(is.na(l$psi_soil_))
  expect_true(is.na(l$j_))
  
  
  #now set physiology, PPFD_, k_l_max_, and psi_soil_ should be not NA
  
  PPFD = 900
  k_l_max = calc_k_l_max(K_s, huber_value, h)
  psi_soil = 2
  
  l$set_physiology(PPFD, psi_soil = psi_soil, k_l_max = k_l_max)
  
  expect_equal(l$PPFD_, PPFD)
  expect_equal(l$k_l_max_, k_l_max)
  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$j_, l$calc_j())
  
  library(tidyverse)
  
  psi <- 1 #nominated value for water potential for testing vulnerability curve equations only (-MPa)
  expect_equal(l$calc_cond_vuln(psi), calc_k_l(psi, k_l_max, b, c) / k_l_max, tolerance = 1e-8)
  
  
  #for situations where psi_soil is < than psi_crit and psi_stem is greater than psi_soil
  psi_stem <- psi_soil+1 #stem water potential (-MPa)
  expect_true(l$calc_E_supply(psi_stem) > 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is less than psi_soil
  psi_stem <- psi_soil-1 #stem water potential (-MPa)
  expect_true(l$calc_E_supply(psi_stem) < 0)
  
  #for situations where psi_soil is < than psi_crit and psi_stem is equal to psi_soil
  psi_stem <- psi_soil #stem water potential (-MPa)
  expect_true(l$calc_E_supply(psi_stem) == 0)
  
  #for situations where psi_stem exceeds psi_crit + tolerance
  expect_error(l$calc_E_supply(psi_crit+psi_crit*0.1), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #for situations where psi_soil exceeds psi_crit + tolerance
  
  psi_soil = psi_crit + psi_crit*0.1
  l$set_physiology(PPFD, psi_soil = psi_soil, k_l_max = k_l_max)
  psi_stem = psi_soil 
  
  expect_error(l$calc_E_supply(psi_stem), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  expect_error(l$calc_E_supply(psi_stem + 1), "Extrapolation disabled and evaluation point outside of interpolated domain.")
  
  #test that fast E supply calculation is closely approximating full integration
  psi_soil = 0
  l$set_physiology(PPFD, psi_soil = psi_soil, k_l_max = k_l_max)
  psi_stem = psi_soil + 3
  
  expect_equal(l$calc_E_supply(psi_stem), l$calc_E_supply_full_integration(psi_stem))
  
  #test that conversion between psi and E works properly
  
  expect_equal(l$calc_psi_stem_ci(l$calc_E_supply(psi_stem)), psi_stem)

  atm_vpd <- 2 #atmopsheric vapour pressure deficit (kPa)
  
  expect_equal(l$calc_g_c(psi_stem), calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b), tolerance = 1e-5)

  c_i = 30 #intra-cellular carbon dioxide parital pressure (Pa)

  c2 = 13.13652
  
  expect_equal(l$calc_A_lim_one_line(c_i), l$calc_j() / 4 * ((c_i - 42.75 * 0.1013) /
                                                               (c_i + 13.13652)))
 
  #for situations where psi stem is lower than psi soil
  l$calc_assim_gross_one_line(psi_soil - 1)
  
  #for situations where psi stem is same as psi soil
  l$calc_assim_gross_one_line(psi_soil)
  
  #for situations where psi stem is higher than psi soil
  l$calc_assim_gross_one_line(psi_soil+ 1)
  
  ca <- 40 #atmospheric carbon dioxide partial pressure
  

  expect_equal(l$calc_A_j(c_i), calc_A_j(c_i = c_i, PPFD = PPFD, vcmax = vcmax))

  expect_equal(l$calc_A_lim(PPFD, c_i), calc_A_lim(c_i, vcmax, PPFD))

  expect_equal(l$calc_assim_gross(PPFD, psi_soil, psi_stem, k_l_max), calc_ben_gross(psi_soil = psi_soil, psi_stem = psi_stem,  b = b, c = c, k_l_max = k_l_max, PPFD = PPFD, vcmax = vcmax, ca =ca, atm_vpd = atm_vpd), tolerance = 1e-4)

  expect_equal(l$calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max), calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_stem,  b = b, c = c, k_l_max = k_l_max))

  expect_equal(l$calc_profit_Sperry(PPFD, psi_soil, psi_stem, k_l_max),calc_profit(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca), tolerance = 1e-4)

  expect_equal(l$E,calc_E_supply(psi_stem, psi_soil, k_l_max = k_l_max, b=b, c=c)$value)

  expect_equal(l$ci, solve_ci(vcmax, PPFD, calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b), ca)[[1]], tolerance = 1e-4)

  beta = 15000 #cost term for Bartlett profit equation (umol m^-3 stem s^-1)
  beta_2 = 1 #shape term for Bartlett profit equation
  
  expect_equal(l$calc_hydraulic_cost_Bartlett(psi_soil, psi_stem, k_l_max),calc_hydraulic_cost_bartlett(psi_soil, psi_stem, k_l_max, b, c, beta, beta_2, huber_value = h_v, height = h*eta_c))

  expect_equal(l$calc_profit_Bartlett(PPFD, psi_soil, psi_stem, k_l_max),
               calc_profit_bartlett(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h*eta_c), tolerance = 1e-4)

  expect_equal(l$optimise_psi_stem_Bartlett(PPFD, psi_soil, k_l_max),
  find_opt_psi_stem_bartlett(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h*eta_c)$psi_stem_root, tolerance = 1e-3)

  l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
  
  l$calc_profit_Sperry_ci(900, 0,  23.03169, k_l_max)
  plot(map_dbl(seq(0, 3.54, length.out = 100), ~l$calc_profit_Sperry(900, 0, .x, k_l_max))~seq(0, 3.54, length.out = 100))
  l$optimise_psi_stem_Sperry_Newton(900, 0, k_l_max)
  l$calc_profit_Sperry(900, 0, 0.8363347, k_l_max)
  l$optimise_ci_Sperry_Newton_recall(900, 0, k_l_max)
  l$optimise_psi_stem_Sperry_Newton_recall(900, 0, k_l_max)
  l$optimise_ci_Sperry_Newton(900, 0, k_l_max)
  l$optimise_psi_stem_Sperry(900, 0, k_l_max)
  l$optim
  
  e_max = l$calc_E_supply(k_l_max, 0, psi_crit) - l$calc_E_supply(k_l_max, 0, 0) 
  gc_max = e_max / (1.6 * atm_vpd / kg_2_mol_h20 / atm_kpa)
  map_dbl(seq(0, 40, 0.1), ~calc_A_lim(.x, 100, 900) * umol_per_mol_2_mol_per_mol - (gc_max * (ca - .x) / (atm_kpa * 1000)))
  
  calc_A_lim(31.8241, 100, 900) * umol_per_mol_2_mol_per_mol/((ca - 31.8241) / (atm_kpa * 1000))  
  0.2620655 * (1.6 * atm_vpd / kg_2_mol_h20 / atm_kpa)
  
  map_dbl(seq(0, , 0.1), ~l$calc_profit_Sperry_ci(900, 0, .x, k_l_max))
  
  ##### when psi soil is equal to psi stem
  
  expect_equal(l$calc_E_supply(k_l_max,psi_soil, psi_soil), calc_E_supply(psi_soil, psi_soil, k_l_max = k_l_max, b=b, c=c)$value)
  
  expect_equal(l$calc_E_supply(k_l_max,psi_soil, psi_soil), 0)
  
  expect_equal(l$calc_g_c(psi_soil, psi_soil, k_l_max), calc_g_c(psi_stem = psi_soil, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b))

  expect_equal(l$calc_g_c(psi_soil, psi_soil, k_l_max), 0)
  
  expect_equal(l$calc_assim_gross(PPFD, psi_soil, psi_soil, k_l_max),calc_ben_gross(psi_soil = psi_soil, psi_stem = psi_soil, atm_vpd= atm_vpd, k_l_max =k_l_max, b = b, c = c, PPFD = PPFD, vcmax = vcmax, ca = ca))
  
  expect_equal(l$calc_profit_Sperry(PPFD, psi_soil, psi_soil, k_l_max),calc_profit(psi_soil, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca))
  
  expect_equal(l$calc_hydraulic_cost_Sperry(psi_soil, psi_soil, k_l_max),calc_hydraulic_cost(psi_soil, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca))
  
  
  
  expect_equal(l$calc_hydraulic_cost_bartlett(psi_soil, psi_soil, k_l_max),calc_hydraulic_cost_bartlett(psi_soil, psi_soil, k_l_max, b, c, beta, beta_2, huber_value = h_v, height = h*eta_c))
  
  expect_equal(l$calc_profit_bartlett(PPFD, psi_soil, psi_soil, k_l_max),
               calc_profit_bartlett(psi_soil, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h*eta_c))
  
 
  ###### when psi soil is greater than psi crit
  
  psi_soil = psi_crit + 1
  
  
  expect_equal(l$optimise_psi_stem_Sperry_Newton(PPFD, psi_soil, k_l_max), psi_soil)
  
  l$calc_assim_gross(PPFD, psi_soil, psi_soil, k_l_max)
  }
)

