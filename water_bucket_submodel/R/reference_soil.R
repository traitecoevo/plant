LH2o <- 2.501e6 #latent heat of water J/kg
MH2O <- 18.02 #molar mass of water g mol^-1
conv.FROM_MILI <- 1e-3 
Cp <- 29.29 #heat capacity of dry air J mol^-1 K^-1
Patm <- 101.325 #atmopsheric pressure at sea level kPa
g_2_kg <- 1/1000
kg_2_m3 <- 1/1000
Pa_2_MPa <- 1/1e6
sec_2_hr <- 60*60
umol_2_mol <- 1/(1e6)
g <- 9.8
pw <- 1000


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

calc_LAI <- function(LAI, lma, r_l, fr_l, a_bio, k_l, a_i = calc_ai(psi_s, E, LAI)){
  a_bio*(a_i*LAI-lma*LAI*r_l)*fr_l/lma-k_l*LAI
}

calc_abs_psi <- function(theta, theta_fc, theta_wp, theta_sat){
  calc_apsi(theta_fc, theta_wp)*(theta/theta_sat)^(-calc_npsi(theta_fc, theta_wp))
}

calc_avg_theta <- function(theta){
  (theta[1]*z[1] + theta[2]*z[2])/(z[1] + z[2])
}

calc_avg_theta_r_frac <- function(theta, r_frac){
  (theta[1]*r_frac + theta[2]*(1 -r_frac))
}

calc_npsi <- function(theta_fc, theta_wp){
  -((log(1500/33))/(log(theta_wp/theta_fc)))
}

calc_apsi <- function(theta_fc, theta_wp){
  1.5e6 * theta_wp^(calc_npsi(theta_fc, theta_wp))
}

calc_D_sat <- function(theta_fc, theta_wp, k_sat, theta_sat){
  calc_apsi(theta_fc, theta_wp)*calc_npsi(theta_fc, theta_wp)*(k_sat/(g*pw))*theta_sat^(-(calc_npsi(theta_fc, theta_wp)+1)) 
}


calc_D <- function(a_psi, n_psi, K_sat, theta_sat, theta){
  calc_D_sat(a_psi, calc_npsi(theta_fc, theta_wp), K_sat, theta_sat)*(theta/theta_sat)^(calc_npsi(theta_fc, theta_wp) + 2)
}

set_params_soil <- function(R=1/365/24, #m/h
                       theta_sat = 0.485, #saturated water content
                       b_infil=8, #shape factor for infiltration (unitless)
                       ksat=3.388889e-06, #soil hydraulic conductivity m/s^-1
                       amp=1/365/24, #amplitude for rainfall sine wave
                       sine_height=1/365/24, #height of sine wave
                       frequency=8760,
                       I_constant = 1,
                       air_temp_c = 25,#deg C
                       theta_fc = 0.321,#field capacity of soil
                       theta_wp = 0.137,#wilting point of soil
                       swf = 1, # sensitivity power factor on empirical soil wetness curve
                       r_soil = 0.01, #soil resistance to water infil
                       mean_PPFD = 1079.901, #umol/m2/s
                       gb = 0.01, #m/s
                       gc_max = 0.0025, #m/s
                       evap_routine = c(1,2),
                       VPD = 1,
                       h_v = 0.000157, #m^2 sapwood area m^-2 leaf area
                       h = 5, #m
                       vcmax = 100, #umol m^-2 s^-1
                       PPFD = 900, #
                       latitude = -23, 
                       day = 1, 
                       k = 0.5, 
                       E = 1,
                       ca = 40,
                       fr_l = 0.1, #unitless
                       k_l = 0.5/365/24, #leaf turnover rate h^-1
                       r_l = 148.8409/365/24, #leaf respiration rate mol h^-1 kg^-1
                       a_bio = 2.45e-2, # Biomass per mol carbon kg mol^-1
                       lma = 0.1,
                       c = 2.04,
                       p50 = 1.731347,
                       beta1 = 15000,
                       beta2 = 1,
                       psi_soil = 1,
                       VPD_9 = 2,
                       diff_temp = 5,
                       diff_VH = 1.2) {
  
  expand_grid(
    R=R,
    theta_sat = theta_sat,
    b_infil=b_infil,
    ksat=ksat,
    amp=amp,
    sine_height=sine_height,
    frequency=frequency,
    I_constant = I_constant,   
    air_temp_c = air_temp_c,
    theta_fc = theta_fc,
    theta_wp = theta_wp,
    swf = swf, # sensitivity power factor on empirical soil wetness curve
    r_soil = r_soil, #soil resistance to water infil
    mean_PPFD = mean_PPFD,
    gb = gb,
    gc_max = gc_max,
    VPD = VPD,
    h_v = h_v, #m^2 sapwood area m^-2 leaf area
    h = h, #m
    vcmax = vcmax, #-MPa
    latitude = latitude, 
    day = day, 
    E = E,
    ca = ca,
    k = k, # light extinction coefficient
    fr_l = fr_l, #unitless
    k_l = k_l, #leaf turnover rate yr^-1
    r_l = r_l, #leaf respiration rate mol y^-1 kg^-1
    a_bio = a_bio, # Biomass per mol carbon kg mol^-1
    lma = lma,
    c = c,
    p50 = p50,    #kPa C-1    
    #kPa C^-1
    #J mol-1
    beta1 = beta1,
    beta2 = beta2,
    psi_soil = psi_soil,
    VPD_9 = VPD_9,
    diff_temp = diff_temp,
    diff_VH = diff_VH) %>%
    mutate(K_s = p_50_2_K_s(p50), #kg m^-1 s^-1 MPa ^-1 Liu et al. 2010 
           slp = calc_slope_vpsat(air_temp_c),
           gamma = calc_psych(air_temp_c),
           lh = calc_lh(air_temp_c),
           npsi = calc_npsi(theta_fc, theta_wp),
           d_sat = calc_D_sat(theta_fc, theta_wp, ksat, theta_sat),
           b = calc_vul_b(p_50 = p50, c = c),
           k_l_max = calc_k_l_max(K_s, h_v, h),
           psi_crit = calc_psi_crit(b, c),
           VPD_15 = VPD_9 *diff_VH
    )
}


