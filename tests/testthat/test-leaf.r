source("tests/reference_leaf.R")
context("SCM-general")


test_that("Basic functions", {
  
  l <- Leaf()
  
  library(tidyverse)
  
  K_s = 2
  h_v = 0.000157
  h = 5
  
  expect_equal(l$calc_k_l_max(K_s, h_v, h), K_s*h_v/h)
  
  p_50 = 1.731347
  c = 2.04
  
  expect_equal(l$calc_vul_b(p_50, c), p_50/(-log(1-50/100))^(1/c))
  
  k_l_max <- l$calc_k_l_max(K_s, h_v, h)
  psi <- 1
  b <- 2.072101
  c <- 2
  
  expect_equal(l$calc_cond_vuln(psi, k_l_max, b, c), k_l_max*exp(-(psi/b)^c), tolerance = 1e-6)
  
  
  control <- Control()
  e <- make_environment("FF16w") 
  
  #set rel_theta at 0.5
  
  e$set_soil_water_state(1)

  psi_aep <- 0.025
  b_ch <- 6.0
  
  expect_equal(l$calc_soil_water_pot(e, psi_aep, b_ch), psi_aep*1^-b_ch, tolerance = 1e-6)
  
  # Make it rain
  x <- seq(0, 9, 1)
  y <- rep(5, 10)
  
  atm_kpa <- 101.3 #kPa
  kg_to_g = 1000 #g kg^-1
  gh2O_to_molh20 = 1/18.02 #mol g^-1
  atm_vpd <- 1
  
  
  kg_2_mol_h20 <- kg_to_g * gh2O_to_molh20
  
  e$set_extrinsic_driver("rainfall", x, y)
  e$set_soil_water_state(1)
  psi_soil <- l$calc_soil_water_pot(e, psi_aep, b_ch)
  
  psi_stem <- psi_soil+1
  
  E_supply <- l$calc_E_supply(k_l_max, b, c, psi_soil, psi_stem)
  
  expect_equal(l$calc_g_c(e, psi_aep, b_ch, psi_stem, k_l_max, p_50, c, b, atm_kpa, kg_2_mol_h20, atm_vpd), atm_kpa*E_supply*kg_2_mol_h20/atm_vpd/1.6)
  
  vcmax = 100
  c_i = 30
  kc_25 <- 404.9
  ko_25 <- 278400
  gamma_25 <- 42.75 
  umol_per_mol_2_mol_per_mol <- 1e-6 # mol umol ^-1
  kPa_2_Pa <- 1000 #Pa/ kPa
  atm_o2_kpa <- 21 #kPa
  
  umol_per_mol_2_Pa <- atm_kpa*kPa_2_Pa*umol_per_mol_2_mol_per_mol
  km_25 <- (kc_25*umol_per_mol_2_Pa)*(1 + (atm_o2_kpa*kPa_2_Pa)/(ko_25*umol_per_mol_2_Pa))
  
  expect_equal(l$calc_A_c(c_i, vcmax, gamma_25, umol_per_mol_2_Pa, km_25), (vcmax * (c_i - gamma_25*umol_per_mol_2_Pa))/(c_i + km_25))
  
  
  PPFD <- 900
  vcmax_25_2_jmax_25 <- 1.67
  curv_fact <- 0.90
  a <- 0.3
  ca <- 40
  
  expect_equal(l$calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i), ((a*PPFD + vcmax*vcmax_25_2_jmax_25 - sqrt((a*PPFD+vcmax*vcmax_25_2_jmax_25)^2-4*curv_fact*a*PPFD*vcmax*vcmax_25_2_jmax_25))/(2*curv_fact))/4 * ((c_i - gamma_25*umol_per_mol_2_Pa)/(c_i + 2*gamma_25*umol_per_mol_2_Pa)))
  
  A_c <- l$calc_A_c(c_i, vcmax, gamma_25, umol_per_mol_2_Pa, km_25)
  A_j <- l$calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i)
  
  expect_equal(l$calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i, km_25), (A_c + A_j  - sqrt((A_c + A_j)^2-4*0.98*A_c*A_j))/(2*0.98))
  
  g_c <- atm_kpa*E_supply*kg_2_mol_h20/atm_vpd/1.6
  

  expect_equal(l$calc_hydraulic_cost(e, psi_stem, k_l_max, b, c, psi_aep, b_ch), calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_stem,  b = b, c = c, k_l_max = k_l_max), tolerance = 1e-6)
  
  expect_equal(l$calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, e, psi_aep, b_ch, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, 1), calc_ben_gross(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c, vcmax= vcmax, PPFD = PPFD), tolerance = 1e-6)
  
  expect_equal(l$calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, e, psi_aep, b_ch, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, -log(0.05)*b/c), calc_profit(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd), tolerance = 1e-6)
  
  expect_equal(l$optimise_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, e, psi_aep, b_ch, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, -log(0.05)*b/c), find_opt_psi_stem(psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)$profit_root, tolerance = 1e-6)

  expect_equal(l$optimise_profit_gss(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, e, psi_aep, b_ch, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, -log(0.05)*b/c), find_opt_psi_stem(psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)$profit_root, tolerance = 1e-6)
  
  
  }
)


psi_soil = 0
b <- 2.072101
c <- 2

K_s = 2
h_v = 0.000157
h = 5
k_l_max <- calc_k_l_max(K_s, h_v, h)

vcmax = 100
PPFD = 900

atm_vpd = 1

gr = (sqrt(5) + 1) / 2

a_bound <- psi_soil
b_bound <- -log(0.05)*b/c
c_bound = b_bound - (b_bound-a_bound)/gr
d_bound = a_bound + (b_bound-a_bound)/gr

tol <- 1e-5

while(abs(b_bound-a_bound) > tol){
  if (calc_profit(c_bound, psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd) > calc_profit(d_bound, psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)){
    b_bound = d_bound
  } else {
    a_bound = c_bound
  }
  c_bound = b_bound - (b_bound-a_bound)/gr
  d_bound = a_bound + (b_bound-a_bound)/gr
}

calc_profit(a_bound+0.01, psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)

find_opt_psi_stem(0.1, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)$profit_root 
0.02500338


psi_stem

find_opt_psi_stem(c_bound, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)$profit_root > find_opt_psi_stem(d_bound, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)$profit_root
calc_profit(0.02500338, psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)


psi_stem <- seq(psi_soil, -log(0.05)*b/c, length.out = 100)

calc_profit_func <- function(psi_stem_input){
  calc_profit(psi_stem_input, psi_soil, k_l_max, vcmax, PPFD, b, c, -log(0.05)*b/c, atm_vpd)
}  


map(psi_stem, calc_profit_func) %>%
  unlist() -> profits

plot(profits~psi_stem)
