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
  instant_through_time <- seq(given_day_start, given_day_start+1, length.out = 200)
  PPFD_instant_mol_m2_day <- PAR_given_solar_angle(solar_angle = solar_angle(latitude = (latitude), decimal_day_time = instant_through_time))
  tibble(PPFD_umol_m2_s_instant = PPFD_instant_mol_m2_day /(60*60*24) * 1E6, instant_of_day = instant_through_time - given_day_start, day = given_day_start, latitude = latitude, ...)
}

#################


#calculate leaf specific conductance
integrate_E_supply <- function(x, k_l_max, b, c){
  k_l_max*exp(-(x/b)^c)
}

#calculate transpiration by intergrating across conductivity curve
calc_E_supply <- function(psi_stem, psi_soil, ...){
  integrate(integrate_E_supply, lower = psi_soil, upper = psi_stem, abs.tol = 1e-6 , ...)
}


calc_g_c <- function(psi_stem, psi_soil, atm_vpd,...){
  g_w = atm_kpa*calc_E_supply(psi_stem = psi_stem, psi_soil = psi_soil, ...)$value*kg_2_mol_h20/atm_vpd
  g_w/1.6
}

calc_A_c <- function(c_i, vcmax){
  (vcmax * (c_i - gamma_25*umol_per_mol_2_Pa))/(c_i + km_25)
}

calc_j <- function(PPFD, vcmax){
  jmax <- vcmax*vcmax_25_2_jmax_25
  
  (a*PPFD + jmax - sqrt((a*PPFD+jmax)^2-4*curv_fact*a*PPFD*jmax))/(2*curv_fact)
}

calc_A_j <- function(c_i, ...){
  calc_j(...)/4 * ((c_i - gamma_25*umol_per_mol_2_Pa)/(c_i + 2*gamma_25*umol_per_mol_2_Pa))
}

calc_A_lim <- function(c_i, vcmax, PPFD){
  A_j = calc_A_j(c_i = c_i, PPFD = PPFD, vcmax = vcmax)
  A_c = calc_A_c(c_i = c_i, vcmax = vcmax)
  (A_c + A_j  - sqrt((A_c + A_j)^2-4*0.98*A_c*A_j))/(2*0.98)
}

diff_ci <- function(x, g_c, vcmax, PPFD, ca){
  calc_A_lim(c_i = x, vcmax = vcmax, PPFD = PPFD)*umol_per_mol_2_mol_per_mol - (g_c*(ca - x)/(atm_kpa*kPa_2_Pa))
}

solve_ci <- function(vcmax, PPFD, g_c, ca){
  ci_root <- uniroot(f = diff_ci,interval = c(0,ca), g_c = g_c, vcmax = vcmax, PPFD = PPFD, ca = ca, extendInt = 'yes', tol = 1e-8)$root
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

calc_hydraulic_cost_bartlett <- function(psi_soil, psi_stem, k_l_max, b, c, beta, beta_2, huber_value, height){
  beta*huber_value*height*(1 - integrate_E_supply(psi_stem, k_l_max, b, c)/integrate_E_supply(psi_soil, k_l_max, b, c))^beta_2
}


calc_profit_bartlett <- function(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca, beta, beta_2, huber_value, height){
  
  calc_ben_gross(psi_stem = psi_stem, psi_soil = psi_soil, vcmax = vcmax, PPFD = PPFD, k_l_max = k_l_max, b = b, c = c, atm_vpd = atm_vpd, ca = ca)-
    calc_hydraulic_cost_bartlett(psi_soil = psi_soil, psi_stem = psi_stem, k_l_max = k_l_max, b = b, c = c, beta = beta, beta_2 = beta_2, huber_value = huber_value, height = height)
}

calc_upper <- function(psi_crit,...){
  psi_crit
}

calc_lower <- function(psi_soil,...){
  psi_soil
}

find_opt_psi_stem <- function(psi_soil, k_l_max, vcmax , PPFD, b, c, psi_crit, atm_vpd, ca){
  
  if(psi_soil < psi_crit){
  psi_stem_root <- ifelse(isTRUE(all.equal(PPFD, 0)), 
                          psi_soil,
                          optimise(f = calc_profit, lower= psi_soil, upper = psi_crit, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, maximum = TRUE, tol = 1e-5)$maximum)
  profit_root <- calc_profit(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)
  cost_root <- calc_hydraulic_cost(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, b = b, c = c)
  gc_root <- calc_g_c(psi_stem = psi_stem_root, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c)
  E_root <- calc_E_supply(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, b = b, c = c)$value
  ci_root <- solve_ci(vcmax = vcmax, PPFD = PPFD, g_c = gc_root, ca = ca)[[1]]
  }
  
  else {
    psi_stem_root <- psi_soil
    profit_root <- calc_profit(psi_stem = psi_soil, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)
    cost_root <- calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_soil, k_l_max = k_l_max, b = b, c = c)
    gc_root <- calc_g_c(psi_stem = psi_stem_root, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c)
    E_root <- calc_E_supply(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, b = b, c = c)$value
    ci_root <- solve_ci(vcmax = vcmax, PPFD = PPFD, g_c = gc_root, ca = ca)[[1]]
  }
  
  tibble(psi_stem_root, profit_root, cost_root, gc_root, E_root, ci_root)
}

find_opt_psi_stem_bartlett <- function(psi_soil, k_l_max, vcmax , PPFD, b, c, psi_crit, atm_vpd, ca, beta, beta_2, huber_value = huber_value, height = height){
  
  if(psi_soil < psi_crit){
    psi_stem_root <- ifelse(isTRUE(all.equal(PPFD, 0)), 
                            psi_soil,
                            optimise(f = calc_profit_bartlett, lower= psi_soil, upper = psi_crit, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, beta = beta, beta_2 = beta_2, huber_value = huber_value, height = height, maximum = TRUE, tol = 1e-5)$maximum)
    profit_root <- calc_profit_bartlett(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, beta = beta, beta_2 = beta_2, huber_value = huber_value, height = height)
  
  }
  
  else {
    psi_stem_root <- psi_soil
    profit_root <- calc_profit_bartlett(psi_stem = psi_soil, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, beta = beta, beta_2 = beta_2, huber_value = huber_value, height = height)
  }
  
  tibble(psi_stem_root, profit_root)
}


set_params <- function(K_s = 2,
                       h_v = 0.000157,
                       h = 5,
                       atm_vpd = 1,
                       vcmax = 100,
                       PPFD = 900,
                       psi_soil = 1,
                       p50 = NULL,
                       b = NULL,
                       psi_stem = 1,
                       latitude = -23, 
                       day = 1, 
                       k = 0.5, 
                       E = 1
) {
  
  c <- 2.04
  
  
  if(is.null(p50)) p50 <- K_s_2_p_50(K_s)
  
  if(is.null(b)) b <- calc_vul_b(p_50 = p50, c = c)
  
  
  k_l_max <- calc_k_l_max(K_s = K_s, h_v = h_v, h = h)
  psi_crit <- calc_psi_crit(b,c)
  
  
  tibble(K_s = K_s,
         p50 = p50,
         # p88 = p88,
         b = b,
         c = c,
         psi_crit = psi_crit,
         k_l_max = k_l_max) %>%
    expand_grid(h_v = h_v,
                h = h,
                atm_vpd = atm_vpd,
                vcmax = vcmax,
                PPFD = PPFD,
                psi_soil = psi_soil,
                psi_stem = psi_stem,
                latitude = latitude, 
                day = day, 
                k = k, 
                E = E)
}

run_models <- function(params){
  params %>%
    rowwise()%>%
    mutate(find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax , PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd),
           profit = calc_profit(k_l_max = k_l_max, b=b, c=c, psi_soil = psi_soil, atm_vpd=atm_vpd, vcmax=vcmax, PPFD = PPFD, psi_crit = psi_crit, psi_stem = psi_stem),
           raw_cost = calc_hydraulic_cost(k_l_max = k_l_max, b=b, c=c, psi_soil = psi_soil,psi_stem = psi_stem),
           lambda = calc_lambda(k_l_max = k_l_max, b=b, c=c, psi_soil = psi_soil, atm_vpd=atm_vpd, vcmax=vcmax, PPFD = PPFD, psi_crit = psi_crit),
           cost = lambda*raw_cost,
           benefit = profit+cost)
}


opt_psi_stem_gss <- function(psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca){
  gr = (sqrt(5) + 1)/2
  tol = 1e-5
  
  a_bound = psi_soil
  b_bound = psi_crit
  
  c_bound = b_bound - (b_bound-a_bound)/gr
  d_bound = a_bound + (b_bound-a_bound)/gr
  
  while(abs(b_bound - a_bound) > tol){
    if(calc_profit(c_bound, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca) >
       calc_profit(d_bound, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca)) {
      b_bound = d_bound
    } else {
      a_bound = c_bound
    }
    c_bound = b_bound - (b_bound-a_bound)/gr
    d_bound = a_bound + (b_bound-a_bound)/gr
  }
  
  opt_psi_stem = (b_bound+a_bound)/2
  
  calc_profit(k_l_max = k_l_max, b=b, c=c, psi_soil = psi_soil, atm_vpd=atm_vpd, vcmax=vcmax, PPFD = PPFD, psi_crit = psi_crit, psi_stem = opt_psi_stem, ca = ca)
  
}



