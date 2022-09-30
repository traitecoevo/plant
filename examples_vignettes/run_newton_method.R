
example_params <- tibble(mean_temp = 20, temp_diff = 17.8, latitude = 15.3, VP = 1.52, vcmax = 35.7, p_50 = 2.35, huber_value = 0.000250, day = 68.7, h = 8.95, E = 0.0376, psi_soil = 0.643, ca = 42.6) %>%
  expand_input_data() %>% 
  filter(PAR > 1e-10) %>%
  ungroup () %>% 
  pmap(calc_profit_instant) %>% 
  bind_rows

example_params <- example_params %>% slice(198)

vcmax = example_params$vcmax #maximum carboxylation rate, defined by leaf nitrogen (umol m^-2 s^-1) 
p_50 = example_params$p_50 #stem water potential at 50% loss of conductivity
c = 2.04 #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
b = calc_vul_b(p_50, c) #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
psi_crit = calc_psi_crit(b, c) #stem water potential at which conductance is 95%
huber_value = example_params$huber_value #huber value (m^2 sapwood area m^-2 leaf area)
K_s = p_50_2_K_s(p_50) #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
h = example_params$h #height or path length (m)
k_l_max = calc_k_l_max(K_s, huber_value, h)
l <- plant:::Leaf(vcmax = vcmax, p_50 = p_50, c = c, b = b, psi_crit = psi_crit, huber_value = huber_value, K_s = K_s, epsilon_leaf = 0.0001)
l$set_physiology(PPFD = example_params$PAR*example_params$E, psi_soil = example_params$psi_soil, k_l_max = k_l_max, atm_vpd = example_params$VPD_hr, example_params$ca)




diff_value = 0.001; 

opt_ci = NA;

if (is.na(opt_ci)){
  l$get_leaf_states_rates_from_psi_stem(psi_crit);
  
  opt_ci = l$ci - diff_value;
}

ci_initial = opt_ci;

x_1 = ci_initial - diff_value;
x_2 = ci_initial + diff_value;

benefit_ = l$calc_A_lim(x_2);


g_c_ci = (benefit_ * 1e-06 * 101.3 * 1000)/(l$ca_ - x_2); 


E = g_c_ci * 1.6 * l$atm_vpd_ / 55.4939 / 101.3;  

psi_stem = l$convert_E_from_ci_to_psi_stem(E);

cost_ = l$calc_hydraulic_cost_Sperry(psi_stem);
y_2 = benefit_ - l$lambda_*cost_;

y_1 = l$calc_profit_Sperry_ci(x_1);


y_0 = l$calc_profit_Sperry_ci(ci_initial);

first_dev = (y_2 - y_1)/(2*diff_value);
sec_dev = (y_2 - 2*y_0 + y_1)/diff_value^2;


opt_ci = ci_initial -  first_dev/sec_dev;

abs(opt_ci - ci_initial) < (0.0001)
