LH2o <- 2.501e6 #latent heat of water J/kg
MH2O <- 18.02 #molar mass of water g mol^-1
conv.FROM_MILI <- 1e-3 
Cp <- 29.29 #heat capacity of dry air J mol^-1 K^-1
Patm <- 101.325 #atmopsheric pressure at sea level kPa
g_2_kg <- 1/1000
kg_2_m3 <- 1/1000

PAR_2_SW = 1/(4.57 * 0.5) # conversion between PAR to SW (umol m^-2 s^-1 to J m^-2 s^-1)

  
#vpsat in kPa
calc_vpsat <- function(air_temp_c){
  0.61078 * exp(17.27 * air_temp_c / (air_temp_c + 237.3))
}

#kPa C-1
calc_slope_vpsat <- function(air_temp_c){
  (calc_vpsat(air_temp_c + 0.01) - calc_vpsat(air_temp_c))/0.01
}

#J mol-1
calc_lh <- function(ait_temp_c){
  (LH2o - 2.365e3 * ait_temp_c) * MH2O * conv.FROM_MILI
}

#kPa C^-1
calc_psych <- function(air_temp_c){
  Cp * Patm / calc_lh(air_temp_c)
}


set_params_soil <- function(R=1/365/24, #m/h
                       dsat = 1.4/24, #diffusivity m^2/h
                       nd = 6.5,#unitless
                       theta_sat = 0.482, #saturated water content
                       b_infil=8, #shape factor for infiltration (unitless)
                       ksat=12.2, #soil hydraulic conducttivity m/h^-1
                       nk=11.9, #shape factor for loss of conductivity (unitless)
                       amp=1/365/24, #amplitude for rainfall sine wave
                       sine_height=1/365/24, #height of sine wave
                       frequency=8760,
                       I_constant = 1,   
                       air_temp_c = 25,#deg C
                       theta_fc = 0.321,#field capacity of soil
                       theta_wp = 0.137,#wilting point of soil
                       swf = 1, # sensitivity power factor on empirical soil wetness curve
                       r_soil = 0.01, #soil resistance to water infil
                       LAI = 1,#m2 /m2
                       mean_PPFD = 1079.901, #umol/m2/s
                       gb = 0.01, #m/s
                       gc_max = 0.0025, #m/s
                       evap_routine = c(1,2),
                       a_psi = 206, #air entry point (Pa),
                       n_psi = 4.5, #unitless,
                       VPD = 1,
                       K_s = 2, #kg m^-1 s^-1 MPa ^-1 Liu et al. 2010 
                       h_v = 0.000157, #m^2 sapwood area m^-2 leaf area
                       h = 5, #m
                       vcmax = 100, #umol m^-2 s^-1
                       PPFD = 900, #
                       psi_soil = 1, #-MPa
                       p50 = NULL, #-MPa
                       b = NULL,
                       psi_stem = 1, #-MPa
                       latitude = -23, 
                       day = 1, 
                       k = 0.5, 
                       E = 1,
                       ca = 40,
                       fr_l = 0.1, #unitless
                       k_l = 0.5, #leaf turnover rate yr^-1
                       r_l = 148.8409, #leaf respiration rate mol y^-1 kg^-1
                       a_bio = 2.45e-2, # Biomass per mol carbon kg mol^-1
                       lma = 0.1) {
  
  expand_grid(
      R=R,
      dsat = dsat,
      nd = nd,
      theta_sat = theta_sat,
      b_infil=b_infil,
      ksat=ksat,
      nk=nk,
      amp=amp,
      sine_height=sine_height,
      frequency=frequency,
      I_constant = I_constant,   
      air_temp_c = air_temp_c,
      theta_fc = theta_fc,
      theta_wp = theta_wp,
      swf = swf, # sensitivity power factor on empirical soil wetness curve
      r_soil = r_soil, #soil resistance to water infil
      LAI = LAI,
      mean_PPFD = mean_PPFD,
      gb = gb,
      gc_max = gc_max,
      evap_routine = evap_routine,
      a_psi = a_psi,
      n_psi = n_psi,
      VPD = VPD,
      K_s = K_s, #kg m^-1 s^-1 MPa ^-1 Liu et al. 2010 
      h_v = h_v, #m^2 sapwood area m^-2 leaf area
      h = h, #m
      vcmax = vcmax, #-MPa
      p50 = NULL, #-MPa
      b = NULL, #-MPa
      latitude = -23, 
      day = day, 
      E = E,
      ca = ca,
      k = k, # light extinction coefficient
      fr_l = fr_l, #unitless
      k_l = k_l, #leaf turnover rate yr^-1
      r_l = r_l, #leaf respiration rate mol y^-1 kg^-1
      a_bio = a_bio, # Biomass per mol carbon kg mol^-1
      lma = lma) %>% #leaf mass per area  kg m^-2 leaf area 
  rowwise() %>%
  split(., 1:nrow(.))
}

calc_LAI <- function(LAI, lma, r_l, fr_l, a_bio, k_l, a_i = calc_ai(psi_s, E, LAI)){
  a_bio*(a_i*LAI-lma*LAI*r_l)*fr_l/lma-k_l*LAI
}