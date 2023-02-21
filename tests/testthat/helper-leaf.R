################################
#constants
################################

kg_to_g = 1000 #mulitplicative converter from kg to g (g kg^-1)
gh2O_to_molh20 = 1/18.02 #mulitplicative converter from g h2o to mol h2o (mol g^-1)

kg_2_mol_h20 <- kg_to_g * gh2O_to_molh20 #mulitplicative converter from kg h2o to mol h2o (mol kg^-1)

R <- 8.3145 #molar gas constant J mol^-1 K^-1

kc_25 <- 404.9 #umol mol ^-1
ko_25 <- 278400 #umol mol ^-1
gamma_25 <- 42.75 #umol mol ^-1
Rd_25 <- 1 

atm_o2_kpa <- 21 #partial pressure o2 (kPa)
atm_kpa <- 101.3 #atmospheric air pressure at sea level (kPa)

kPa_2_Pa <- 1000 #multiplicative converter from kPa to Pa (Pa kPa^-1)

umol_per_mol_2_mol_per_mol <- 1e-6 #multiplicative converter from umol per mol to mol per mol (mol umol ^-1)

umol_per_mol_2_Pa <- atm_kpa*kPa_2_Pa*umol_per_mol_2_mol_per_mol  #multiplicative converter from umol per mol to Pa (Pa umol mol^-1)

a <- 0.3 #quantum yield of electron transport (mol photon mol ^-1 e)

vcmax_25_2_jmax_25 <- 1.67 # estimate vcmax from jmax (unitless)

curv_fact <- 0.9 #curvature factor for light limited transport (unitless)

km_25 <- (kc_25*umol_per_mol_2_Pa)*(1 + (atm_o2_kpa*kPa_2_Pa)/(ko_25*umol_per_mol_2_Pa)) #Pa

K_s_0 <- 2 #global central value K_s (Liu et al. 2019)


################################ functions below calculate values outside of the leaf model

#estimate sapwood-specific conductivity from p_50, estimated visually from Liu et al. 2019
p_50_2_K_s <- function(p50){
  1.638999*(p50/2)^(-1.38)
}

#estimate p_50 (-MPa) from sapwood-specific conductivity (kg m^-1 s^-1 MPa), estimated visually from Liu et al. 2019
K_s_2_p_50 <- function(K_s){
  1.731347*(K_s/K_s_0)^(-0.7246377)
}


#calculate leaf-specific conductance (kg m^-2 s^-1 MPa) from sapwood-specific conductivity (kg m^-1 s^-1 MPa), huber value (m^2 sapwood area m^-2 leaf area) and height (m)
calc_k_l_max <- function(K_s, h_v, h){
  K_s * h_v / h
}

#calculate sapwood volume per leaf area (m^3 sapwod m^-2 leaf area) 
calc_sapwood_volume_per_leaf_area <- function(h_v, h){
  h_v*h
}

#calculate b (-MPa) from c (unitless) and p_50 (-MPa)
calc_vul_b <- function(p_50, c){
  num <- p_50
  den <- (-log(1-50/100))^(1/c)
  num/den
}

#calculate stem water potentail at 95% conductivity lost
calc_psi_crit <- function(b,c) {
  b*(log(1/0.05))^(1/c)
}

#calculate PPFD throughout day
calc_PPFD_instant <- function(latitude, day, ...) {
  given_day_start <- floor(day)
  instant_through_time <- seq(given_day_start, given_day_start+1, length.out = 201)
  PPFD_instant_mol_m2_day <- PAR_given_solar_angle(solar_angle = solar_angle(latitude = (latitude), decimal_day_time = instant_through_time))
  outputs <- tibble(PPFD_umol_m2_s_instant = PPFD_instant_mol_m2_day /(60*60*24) * 1E6, instant_of_day = instant_through_time - given_day_start, day = given_day_start, latitude = latitude, ...)
  
  outputs %>%
    filter(PPFD_umol_m2_s_instant != 0) %>%
    summarise(sun_up = min(instant_of_day), sun_down = max(instant_of_day), day_length = sun_down - sun_up) -> summary_outputs
  
  bind_cols(outputs, summary_outputs)
}

#################


#calculate leaf specific conductance
integrate_E_supply <- function(x, b, c, k_l_max){
  k_l_max*exp(-(x/b)^c)
}

#calculate transpiration by intergrating across conductivity curve
calc_E_supply <- function(psi_stem, psi_soil, ...){
  integrate(integrate_E_supply, lower = psi_soil, upper = psi_stem, abs.tol = 1e-6 ,stop.on.error = TRUE,  ...)
}

#calculate stomatal conductance to h20 (g_w) and co2 (g_c) from E supply and VPD
calc_g_c <- function(psi_stem, psi_soil, atm_vpd,...){
  g_w = atm_kpa*l$transpiration(psi_stem)*kg_2_mol_h20/atm_vpd
  # g_w = atm_kpa*calc_E_supply(psi_stem = psi_stem, psi_soil = psi_soil, ...)$value*kg_2_mol_h20/atm_vpd
  g_w/1.6
}

#calculate carboxylation limited assimilation rate
calc_A_c <- function(c_i, vcmax){
  (vcmax * (c_i - gamma_25*umol_per_mol_2_Pa))/(c_i + km_25)
}

#calculate electron transport rate
calc_j <- function(PPFD, vcmax){
  jmax <- vcmax*vcmax_25_2_jmax_25
  
  (a*PPFD + jmax - sqrt((a*PPFD+jmax)^2-4*curv_fact*a*PPFD*jmax))/(2*curv_fact)
}


#calculate electron transport rate-limited assimilation
calc_A_j <- function(c_i, ...){
  calc_j(...)/4 * ((c_i - gamma_25*umol_per_mol_2_Pa)/(c_i + 2*gamma_25*umol_per_mol_2_Pa))
}

#calculate co-limited assimilation rate
calc_A_lim <- function(c_i, vcmax, PPFD){
  A_j = calc_A_j(c_i = c_i, PPFD = PPFD, vcmax = vcmax)
  A_c = calc_A_c(c_i = c_i, vcmax = vcmax)
  (A_c + A_j  - sqrt((A_c + A_j)^2-4*0.98*A_c*A_j))/(2*0.98)
}

#calculated assimilation based on just A_j, A=gc*(ca - ci) equation to be re-arranged to find ci without root-finding
calc_A_lim_one_line <- function(c_i,...) {
  # estimated parameter to smooth A_j to middle of A_j and A_c
  c2 = 13.13652;
  
  calc_j(...) / 4 * ((c_i - gamma_25 * umol_per_mol_2_Pa) /
       (c_i + c2));
}


diff_ci <- function(x, g_c, vcmax, PPFD, ca){
  calc_A_lim(c_i = x, vcmax = vcmax, PPFD = PPFD)*umol_per_mol_2_mol_per_mol - (g_c*(ca - x)/(atm_kpa*kPa_2_Pa))
}

solve_ci <- function(vcmax, PPFD, g_c, ca){
  ci_root <- uniroot(f = diff_ci,interval = c(0,ca), g_c = g_c, vcmax = vcmax, PPFD = PPFD, ca = ca, extendInt = 'no', tol = 1e-8)$root
  A_lim_root <- calc_A_lim(c_i = ci_root, vcmax = vcmax, PPFD = PPFD)
  list(ci_root, A_lim_root)
}

calc_k_l <- function(psi ,k_l_max, b, c){
  k_l_max*exp(-(psi/b)^c)
}

calc_ben_gross <- function(psi_stem, psi_soil, atm_vpd, k_l_max, b, c, PPFD, vcmax, ca) {
  
  g_c <- calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c)
  
  as.numeric(solve_ci(g_c = g_c, vcmax = vcmax, PPFD = PPFD, ca = ca)[2]) 
}

calc_ben_gross_one_line <- function(psi_stem, psi_soil, atm_vpd, k_l_max, b, c, PPFD, vcmax, ca) {
  
  g_c <- calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c)
  c2 = 13.13652;
  j_ = calc_j(PPFD, vcmax)

  first_term = 8 * j_ * umol_per_mol_2_mol_per_mol * (atm_kpa*kPa_2_Pa) * g_c * (-ca + c2 + 2 * gamma_25 * umol_per_mol_2_Pa);
  second_term = 16 * g_c^2;
  third_term = (ca + c2)^2;
  fourth_term = j_^2 * umol_per_mol_2_mol_per_mol^2 * (atm_kpa*kPa_2_Pa)^2;
  fifth_term = 4*ca*g_c;
  sixth_term = 4*c2*g_c;
  seventh_term = j_*umol_per_mol_2_mol_per_mol*(atm_kpa*kPa_2_Pa);
  eigth_term = 8*g_c;
  
  ci = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;
  A_lim = calc_A_lim_one_line(c_i =ci, vcmax = vcmax, PPFD = PPFD);
  E = g_c * 1.6 * atm_vpd / kg_2_mol_h20 / atm_kpa;
  psi = psi_stem;
  
  A_lim
}

find_max_ci_one_line <- function(psi_crit, psi_soil, atm_vpd, ...) {
  g_c = calc_g_c(psi_crit, psi_soil, atm_vpd, k_l_max = k_l_max, b = b, c = c)  
  
  c2 = 13.13652
  j = calc_j(...)
  first_term = 8 * j * umol_per_mol_2_mol_per_mol * (atm_kpa*kPa_2_Pa) * g_c * (-ca + c2 + 2 * gamma_25 * umol_per_mol_2_Pa);
  second_term = 16 * g_c^2;
  third_term = (ca + c2)^2;
  fourth_term = j^2 * umol_per_mol_2_mol_per_mol^2 * (atm_kpa*kPa_2_Pa)^2;
  fifth_term = 4*ca*g_c;
  sixth_term = 4*c2*g_c;
  seventh_term = j*umol_per_mol_2_mol_per_mol*(atm_kpa*kPa_2_Pa);
  eigth_term = 8*g_c;
  
  ci = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term    
  return(ci)
}


diff_psi <- function(x, psi_soil, E,...){
    calc_E_supply(psi_stem = x, psi_soil = 0,...)$value - calc_E_supply(psi_stem = psi_soil, psi_soil = 0,...)$value - E
}

solve_psi <- function(E, ...){
  psi_root <- uniroot(f = diff_psi,interval = c(psi_soil,psi_crit), E, ..., extendInt = 'no', tol = 1e-8)$root
  psi_root
}

calc_profit_Sperry_ci_one_line <- function(c_i, ...){                                  
  benefit_ = calc_A_lim_one_line(c_i = c_i, vcmax = vcmax, PPFD = PPFD)
  g_c_ci = (benefit_ * umol_per_mol_2_mol_per_mol * atm_kpa * kPa_2_Pa)/(ca - c_i)
  E = g_c_ci * 1.6 * atm_vpd / kg_2_mol_h20 / atm_kpa;  
  psi_stem = solve_psi(E = E, psi_soil = psi_soil, b =b, c= c, ...)
  
  lambda_ = calc_lambda(psi_soil= psi_soil, k_l_max, vcmax = vcmax, PPFD = PPFD, b= b, c= c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca= ca)

  cost_ = calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_stem, ...);
  benefit_ - lambda_*cost_;
}

calc_profit_Sperry_psi_stem_one_line <- function(psi_stem, ...){                                  
    benefit_ = calc_ben_gross_one_line(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c =c , PPFD = PPFD, vcmax = vcmax, ca = ca)

  lambda_ = calc_ben_gross_one_line(psi_stem = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c =c , PPFD = PPFD, vcmax = vcmax, ca = ca) / calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_crit, ...)
  
  cost_ = calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_stem, ...);
  benefit_ - lambda_*cost_;
}


calc_hydraulic_cost <- function(psi_soil, psi_stem, k_l_max, b, c, ...){
  (calc_k_l(psi = psi_soil, k_l_max = k_l_max, b = b, c = c) - calc_k_l(psi = psi_stem, k_l_max = k_l_max, b = b, c = c))
}

calc_lambda <-  function(psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca) {
  calc_ben_gross(psi_stem = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c, PPFD = PPFD, vcmax= vcmax, ca = ca) /
    (calc_k_l(psi = psi_soil, k_l_max = k_l_max, b = b, c = c) - calc_k_l(psi = psi_crit, k_l_max = k_l_max, b = b, c = c))
}

calc_profit <- function(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca){
  
  calc_ben_gross(psi_stem = psi_stem, psi_soil = psi_soil, vcmax = vcmax, PPFD = PPFD, k_l_max = k_l_max, b = b, c = c, atm_vpd = atm_vpd, ca = ca)- 
    calc_lambda(psi_soil = psi_soil, vcmax = vcmax, PPFD = PPFD, k_l_max = k_l_max, b = b, c = c, atm_vpd = atm_vpd, psi_crit = psi_crit, ca = ca)*
    calc_hydraulic_cost(psi_stem = psi_stem, psi_soil = psi_soil, k_l_max =  k_l_max, b = b, c = c)
}


opt_ci_gss <- function(psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca){
  gr = (sqrt(5) + 1)/2
  tol = 1e-5
  
  a_bound = gamma_25 * umol_per_mol_2_Pa
  b_bound = find_max_ci_one_line(psi_crit = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, vcmax = vcmax, PPFD = PPFD)
  
  c_bound = b_bound - (b_bound-a_bound)/gr
  d_bound = a_bound + (b_bound-a_bound)/gr
  
  while(abs(b_bound - a_bound) > tol){
    if(calc_profit_Sperry_ci_one_line(c_i = c_bound, k_l_max, b, c) >
       calc_profit_Sperry_ci_one_line(c_i = d_bound, k_l_max, b, c)) {
      b_bound = d_bound
    } else {
      a_bound = c_bound
    }
    c_bound = b_bound - (b_bound-a_bound)/gr
    d_bound = a_bound + (b_bound-a_bound)/gr
  }
  
  opt_ci = (b_bound+a_bound)/2
  
  calc_profit_Sperry_ci_one_line(c_i = opt_ci, k_l_max, b, c)
}


opt_psi_gss <- function(psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca){
  gr = (sqrt(5) + 1)/2
  tol = 1e-5
  
  a_bound = psi_soil
  b_bound = psi_crit
  
  c_bound = b_bound - (b_bound-a_bound)/gr
  d_bound = a_bound + (b_bound-a_bound)/gr
  
  while(abs(b_bound - a_bound) > tol){
    if(calc_profit_Sperry_psi_stem_one_line(c_i = c_bound, k_l_max, b, c) >
       calc_profit_Sperry_psi_stem_one_line(c_i = d_bound, k_l_max, b, c)) {
      b_bound = d_bound
    } else {
      a_bound = c_bound
    }
    c_bound = b_bound - (b_bound-a_bound)/gr
    d_bound = a_bound + (b_bound-a_bound)/gr
  }
  
  opt_ci = (b_bound+a_bound)/2
  
  calc_profit_Sperry_psi_stem_one_line(c_i = opt_ci, k_l_max, b, c)
}


