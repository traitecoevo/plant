################################
#constants
################################

sigma = 5.67e-8 #stefan boltzman constant w m^-2 k^-4
eps = 0.97 #emmisivity (unitless)
c_p = 29.3 #specific heat capacity (J mol^-1 K^-1)
t_C_2_k = 273.15 #conversion from C to K,
kg_to_g = 1000 #g kg^-1
gh2O_to_molh20 = 1/18.02 #mol g^-1

kg_2_mol_h20 <- kg_to_g * gh2O_to_molh20


R <- 8.3145 #molar gas constant J mol^-1 K^-1

kc_25 <- 404.9
ko_25 <- 278400
gamma_25 <- 42.75 
Rd_25 <- 1 

kc_dha <- 79.43 * 1000
ko_dha <- 36.38 * 1000
gamma_dha <- 37.83 * 1000 
Rd_dha <- 46.39 * 1000 

kc_c <- 38.05
ko_c <- 20.30
gamma_c <- 19.02 
Rd_c <- 18.72 


#leunig 2001
vcmax_ha <- 73637  # J mol ^-1
vcmax_hd <- 149252  # J mol ^-1
vcmax_sv <- 486   # J mol ^-1 K ^-1

jmax_ha <- 50300  # J mol ^-1
jmax_hd <- 152044  # J mol ^-1
jmax_sv <- 495   # J mol ^-1 K ^-1


atm_o2_kpa <- 21 #kPa
atm_kpa <- 101.3 #kPa


kPa_2_Pa <- 1000 #Pa/ kPa

umol_per_mol_2_mol_per_mol <- 1e-6 # mol umol ^-1

umol_per_mol_2_Pa <- atm_kpa*kPa_2_Pa*umol_per_mol_2_mol_per_mol

a <- 0.3 #quantum yield of electron transport (mol photon mol ^-1 e)

vcmax_25_2_jmax_25 <- 1.67

curv_fact <- 0.9 #curvature factor

ca <- 40 #pa atm c02

km_25 <- (kc_25*umol_per_mol_2_Pa)*(1 + (atm_o2_kpa*kPa_2_Pa)/(ko_25*umol_per_mol_2_Pa))

K_s_0 <- 2

################################


calc_k_l_max <- function(K_s, h_v, h){
  K_s * h_v / h
}

p_50_2_K_s <- function(p50){
  1.638999*(p50/2)^(-1.38)
}

K_s_2_p_50 <- function(K_s){
  1.731347*(K_s/K_s_0)^(-0.7246377)
}

calc_vul_b <- function(p_50){
  num <- p_50
  den <- (-log(1-50/100))^(1/2.04)
  num/den
}

calc_PPFD_instant <- function(latitude, day, ...) {
  given_day_start <- floor(day)
  instant_through_time <- seq(given_day_start, given_day_start+1, length.out = 200)
  PPFD_instant_mol_m2_day <- PAR_given_solar_angle(solar_angle = solar_angle(latitude = (latitude), decimal_day_time = instant_through_time))
  tibble(PPFD_umol_m2_s_instant = PPFD_instant_mol_m2_day /(60*60*24) * 1E6, instant_of_day = instant_through_time - given_day_start, day = given_day_start, latitude = latitude, ...)
}



#################

integrate_E_supply <- function(x, k_l_max, b, c){
  k_l_max*exp(-(x/b)^c)
}

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

diff_ci <- function(x, g_c, vcmax, PPFD){
  calc_A_lim(c_i = x, vcmax = vcmax, PPFD = PPFD)*umol_per_mol_2_mol_per_mol - (g_c*(ca - x)/(atm_kpa*kPa_2_Pa))
}

solve_ci <- function(vcmax, PPFD, g_c){
  ci_root <- uniroot(f = diff_ci,interval = c(0,40), g_c = g_c, vcmax = vcmax, PPFD = PPFD, extendInt = 'yes', tol = 1e-6)$root
  A_lim_root <- calc_A_lim(c_i = ci_root, vcmax = vcmax, PPFD = PPFD)
  list(ci_root, A_lim_root)
}

calc_psi_crit <- function(b,c) {
  -log(0.05)*b/c
}


calc_k_l <- function(psi ,k_l_max, b, c){
  k_l_max*exp(-(psi/b)^c)
}

calc_ben_gross <- function(psi_stem, psi_soil, atm_vpd, k_l_max, b, c, PPFD, vcmax) {
  
  g_c <- calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c)
  
  as.numeric(solve_ci(g_c = g_c, vcmax = vcmax, PPFD = PPFD)[2])
}

calc_hydraulic_cost <- function(psi_soil, psi_stem, k_l_max, b, c, ...){
  (calc_k_l(psi = psi_soil, k_l_max = k_l_max, b = b, c = c) - calc_k_l(psi = psi_stem, k_l_max = k_l_max, b = b, c = c))
}

calc_lambda <-  function(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd) {
  calc_ben_gross(psi_stem = psi_crit, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c, PPFD = PPFD, vcmax= vcmax) /
    (calc_k_l(psi = psi_soil, k_l_max = k_l_max, b = b, c = c) - k_l_max*0.05)
}

calc_profit <- function(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd){
  
  calc_ben_gross(psi_stem = psi_stem, psi_soil = psi_soil, vcmax = vcmax, PPFD = PPFD, k_l_max = k_l_max, b = b, c = c, atm_vpd = atm_vpd)- 
    calc_lambda(psi_stem = psi_stem, psi_soil = psi_soil, vcmax = vcmax, PPFD = PPFD, k_l_max = k_l_max, b = b, c = c, atm_vpd = atm_vpd, psi_crit = psi_crit)*
    calc_hydraulic_cost(psi_stem = psi_stem, psi_soil = psi_soil, k_l_max =  k_l_max, b = b, c = c)
}

calc_upper <- function(psi_crit,...){
  psi_crit
}

calc_lower <- function(psi_soil,...){
  psi_soil
}

find_opt_psi_stem <- function(psi_soil, k_l_max, vcmax , PPFD, b, c, psi_crit, atm_vpd){
  
  if(psi_soil < psi_crit){
  psi_stem_root <- ifelse(isTRUE(all.equal(PPFD, 0)), 
                          psi_soil,
                          optimise(f = calc_profit, lower= psi_soil, upper = psi_crit, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, maximum = TRUE, tol = 1e-5)$maximum)
  profit_root <- calc_profit(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd)
  cost_root <- calc_hydraulic_cost(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, b = b, c = c)
  gc_root <- calc_g_c(psi_stem = psi_stem_root, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c)
  
  E_root <- calc_E_supply(psi_stem = psi_stem_root, psi_soil = psi_soil, k_l_max = k_l_max, b = b, c = c)$value
  }
  
  else {
    psi_stem_root <- psi_soil
    profit_root <- 0
    cost_root <- 0
    gc_root <- 0
    E_root <- 0
  }
  
  tibble(psi_stem_root, profit_root, cost_root, gc_root, E_root)
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
  
  if(is.null(p50)) p50 <- K_s_2_p_50(K_s)
  
  if(is.null(b)) b <- calc_vul_b(p_50 = p50)
  
  
  k_l_max <- calc_k_l_max(K_s = K_s, h_v = h_v, h = h)
  c <- 2.04
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
