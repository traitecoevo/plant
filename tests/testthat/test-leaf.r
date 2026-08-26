
# NOTE ON ROOT CARBON. The leaf model is purely intensive: the root resistances it
# takes are **per unit leaf area**, and its soil->collar uptake is exactly
# homogeneous in (root carbon / area_leaf) -- verified bit-identical. The tests
# below therefore keep an absolute `root_carbon` and divide by `area_leaf_` at the
# call, so every pinned value keeps the meaning it had when the leaf took the two
# separately. Pass the absolute carbon by mistake and the root system is 20x too
# weak at area_leaf_ = 0.05, which moves the critical-demand collar potential from
# -0.685 to -2.57 MPa while leaving the zero-uptake collar untouched (that one is
# scale-invariant) -- a failure mode worth recognising, because only one of the two
# regression guards catches it.
#
# ⚠️ AND SINCE phylloptim #33 THE LEAF CANNOT CATCH IT AT ALL. `set_physiology`
# takes the resistances rather than the carbon, so a network built from absolute
# carbon is five vectors of plausible positive numbers and nothing in the leaf
# knows the difference. The division below is now the only thing enforcing the
# contract, which is why it stays explicit at every call site rather than being
# folded into `net()`.

# The carbon -> resistance step, which used to happen inside `set_physiology`
# (phylloptim #33). Same arithmetic and the same numbers -- what moved is who calls
# it. It carries the two architecture constants, at 10x the package defaults, which
# is what this file has always used -- they were assigned separately in each of the
# test_that() blocks that needed them and are now stated once.
#
# Namespaced, and NOT reimplemented here, on purpose: the vertical resistance
# scales with dz^2, so a second copy of this arithmetic in plant would be one slip
# away from a silent squared factor that neither package could detect. phylloptim
# is a LinkingTo dependency already; it is in Suggests for this R-level use.
test_beta_R_H <- 3.4e3
test_beta_R_V <- 9.4e4

net <- function(root_carbon_per_leaf_area, soil_depth) {
  phylloptim::root_network_from_carbon(
    root_carbon_per_leaf_area, soil_depth = soil_depth,
    beta_R_H = test_beta_R_H, beta_R_V = test_beta_R_V)
}

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
  curv_fact_elec_trans = 0.7
  a = 0.3
  curv_fact_colim = 0.99
  g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8
  vulnerability_curve_ncontrol = 100
  ci_abs_tol = 1e-6
  ci_niter = 1000
  root_c = 2.65
  root_b = 1.29
  root_psi_crit = root_b * (log(1.0 / 0.05))^(1.0 / root_c)
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  
  #without setting physiology, PPFD_, k_l_max_, and psi_soil_ should all be NA
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$leaf_specific_conductance_max_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$atm_o2_kpa_))
  expect_true(is.na(l$leaf_temp_))
  expect_true(is.na(l$ci_))
  expect_true(is.na(l$stom_cond_CO2_))
  expect_true(is.na(l$assim_colimited_))
  expect_true(is.na(l$transpiration_))
  expect_true(is.na(l$profit_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$hydraulic_cost_))
  expect_true(is.na(l$electron_transport_))
  expect_true(is.na(l$gamma_))
  expect_true(is.na(l$ko_))
  expect_true(is.na(l$kc_))
  expect_true(is.na(l$km_))
  expect_true(is.na(l$R_d_))
  expect_true(is.na(l$vcmax_))
  expect_true(is.na(l$jmax_))
  expect_true(is.na(l$opt_root_psi_))
  expect_true(is.na(l$opt_psi_stem_))
  expect_true(is.na(l$opt_ci_))
  expect_true(is.na(l$E_up_))
  expect_true(is.na(l$soil_number_of_depths_))
  
  
  expect_true(length(l$psi_soil_) == 0)
  expect_true(length(l$soil_depth_) == 0)
  expect_true(length(l$z_soil_mid_) == 0)
  expect_true(length(l$c_r_H_) == 0)
  expect_true(length(l$c_r_V_) == 0)

  
  #now set physiology, PPFD_, k_l_max_, and psi_soil_, atm_vpd_ should be not NA
  
  PPFD = 900
  leaf_specific_conductance_max = K_s*theta/h
  psi_soil = 2
  atm_vpd = 2
  ca = 40
  atm_o2_kpa_ = 21
  leaf_temp_ = 25
  atm_kpa_ = 101.3
  area_leaf_ = 0.05
  root_carbon_ = 1

  # set physiology when inputting just a single soil layer (via soil depth and psi soil)

  l$set_physiology(root_network = net((root_carbon_) / area_leaf_, 1), PPFD = PPFD, psi_soil = psi_soil, soil_depth = 1, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  # The root network is built from the per-leaf-area carbon, so these diagnostics
  # are per leaf area too: the split carbon scales up by 1/area_leaf_ and the two
  # resistances (beta/c_r) scale down by area_leaf_. The physics is unchanged --
  # uptake is exactly homogeneous in (carbon / area_leaf) -- so the old absolute
  # numbers are written here times the factor rather than replaced, which keeps
  # them checkable against the pre-package values.
  expect_equal(l$c_r_V_, (root_carbon_ / area_leaf_) * (1.0/3.0))
  expect_equal(l$c_r_H_, (root_carbon_ / area_leaf_) * (2.0/3.0))
  expect_equal(l$r_R_V, 282000 * area_leaf_)
  expect_equal(l$r_R_V_sum, 282000 * area_leaf_)
  expect_equal(l$r_R_H_min, 5100 * area_leaf_)
  expect_equal(l$PPFD_, PPFD)
  expect_equal(l$leaf_specific_conductance_max_, leaf_specific_conductance_max)
  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$ca_, ca)
  expect_equal(l$atm_vpd_, atm_vpd)
  expect_equal(l$atm_o2_kpa_, atm_o2_kpa_)
  expect_equal(l$leaf_temp_, leaf_temp_)
  expect_equal(l$atm_kpa_, atm_kpa_)
  expect_equal(l$soil_depth_, 1)
  expect_equal(l$soil_number_of_depths_, length(psi_soil))

# throw error when length of psi soil and soil depth do not match. The network is
# built for the two-layer profile so it is CONSISTENT with soil_depth: the mismatch
# under test is psi_soil against soil_depth, and since phylloptim #33 a network of
# the wrong length is a separate error with its own message.

expect_error(l$set_physiology(root_network = net(rep((root_carbon_) / area_leaf_ / 2, 2), c(1,2)), PPFD = PPFD, psi_soil = psi_soil, soil_depth = c(1,2), leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_), "soil_depth and psi_soil must have the same number of elements")

# and a network with more rooted layers than the soil profile has layers is
# rejected too -- it would be an out-of-bounds read in uptake(), not a wrong
# number. Before #33 the agreement came for free from validating root carbon
# against soil_depth.

expect_error(l$set_physiology(root_network = net(rep((root_carbon_) / area_leaf_ / 2, 2), c(1,2)), PPFD = PPFD, psi_soil = psi_soil, soil_depth = 1, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_), "rooted layers but the soil profile has only")

  # test that the inputs to set_physiology which take multiple values are working correctly

  psi_soil = c(1,2)
  soil_depth = c(0.5, 1)
  root_carbon_ = c(1,1)

  l$set_physiology(root_network = net((root_carbon_) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

  expect_equal(l$psi_soil_, psi_soil)
  expect_equal(l$soil_number_of_depths_, length(psi_soil))
  expect_equal(l$soil_depth_, soil_depth)
  expect_equal(l$c_r_V_, c(1/3,1/3) / area_leaf_)


  #generating a new leaf object should wipe the previously stored values
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  
  expect_true(is.na(l$PPFD_))
  expect_true(is.na(l$leaf_specific_conductance_max_))
  expect_true(is.na(l$atm_vpd_))
  expect_true(is.na(l$ca_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$atm_o2_kpa_))
  expect_true(is.na(l$leaf_temp_))
  expect_true(is.na(l$ci_))
  expect_true(is.na(l$stom_cond_CO2_))
  expect_true(is.na(l$assim_colimited_))
  expect_true(is.na(l$transpiration_))
  expect_true(is.na(l$profit_))
  expect_true(is.na(l$lambda_))
  expect_true(is.na(l$hydraulic_cost_))
  expect_true(is.na(l$electron_transport_))
  expect_true(is.na(l$gamma_))
  expect_true(is.na(l$ko_))
  expect_true(is.na(l$kc_))
  expect_true(is.na(l$km_))
  expect_true(is.na(l$R_d_))
  expect_true(is.na(l$vcmax_))
  expect_true(is.na(l$jmax_))
  expect_true(is.na(l$opt_root_psi_))
  expect_true(is.na(l$opt_psi_stem_))
  expect_true(is.na(l$opt_ci_))
  expect_true(is.na(l$E_up_))
  expect_true(is.na(l$soil_number_of_depths_))

  
  #set physiology again for testing 
  psi_soil = 1
  soil_depth = 0.5
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
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
  # The leaf reads the same transport spline from four places and holds a second
  # spline that is its inverse, so "outside its domain" on its own does not say
  # which lookup failed -- and that is the fact localising #576 turned on. Assert
  # the three parts that make the message a diagnosis (which spline, which end,
  # which caller) rather than the whole sentence, which is phylloptim's to word.
  err <- expect_error(l$transpiration(upper_bound_int, psi_stem[1]),
                      "evaluated outside its domain")
  expect_match(conditionMessage(err), "transpiration_from_psi")
  expect_match(conditionMessage(err), "beyond the upper end")
  expect_match(conditionMessage(err), "Leaf::transpiration, at psi_stem")

  #for situations where psi_soil exceeds psi_crit + tolerance
  
  psi_soil = upper_bound_int
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  psi_stem = psi_soil 
  
  # Same spline, same end -- and the SAME caller label, because psi_stem and
  # psi_soil are equal here and `Leaf::transpiration` checks the stem argument
  # first, so it never reaches the upstream one. Asserted rather than left
  # implicit: this case is reached from the soil side, and the label saying
  # `at psi_stem` is the thing that would otherwise read as a bug.
  err <- expect_error(l$transpiration(psi_stem[1], psi_soil[1]),
                      "evaluated outside its domain")
  expect_match(conditionMessage(err), "transpiration_from_psi")
  expect_match(conditionMessage(err), "beyond the upper end")
  expect_match(conditionMessage(err), "Leaf::transpiration, at psi_stem")

  #test that fast E supply calculation is closely approximating full integration
  psi_soil = 0
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  psi_stem = psi_soil + 3
  
  expect_equal(l$transpiration(psi_stem[1], psi_soil[1]), l$transpiration_full_integration(psi_stem[1], psi_soil[1]))
  
  #test that conversion between psi and E works properly
  
  expect_equal(l$transpiration_to_psi_stem(l$transpiration(psi_stem[1], psi_soil[1]), psi_soil[1]), psi_stem[1], tolerance = 1e-5)

  c_i = 30 #intra-cellular carbon dioxide parital pressure (Pa)

  #test a function which retrieves various leaf-level states and rates from a given psi_stem value
  #for situations where psi stem is lower than psi soil

  psi_soil = 2
  
   l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  #note that this scenario should not occur in model anyway
  l$set_leaf_states_rates_from_psi_stem(psi_soil - 1, psi_soil)

  #assimilation becomes dark respiration
  expect_equal(l$assim_colimited_, -l$R_d_)
  #stomatal conductance becomes 0
  expect_equal(l$stom_cond_CO2_, 0)
  #transpiration becomes 0
  expect_equal(l$transpiration_, 0)
  
  #costs 0 when psi_stem == psi_soil == 0
  ## phylloptim 0.6.0 removed the two-argument Sperry hydraulic cost; its
  ## successor is the normalised field hydraulic_cost_norm_, written by the
  ## ProfitMax solve rather than callable at an arbitrary psi_stem (#622).
  ## The hydraulic_cost_TF assertions beside these cover the same property
  ## for the curve plant actually runs.
  #costs positive even when transpiration stream is 0 in hydraulic cost tf
  expect_equal(l$hydraulic_cost_TF(psi_soil) > 0, TRUE)
  
  #when psi stem is greater than psi soil
  psi_soil = 2
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l$set_leaf_states_rates_from_psi_stem(psi_soil + 1, psi_soil)
  

  #assimilation becomes greater than 0 
  expect_equal(l$assim_colimited_ >0, TRUE)
  #stomatal conductance becomes greater than 0 
  expect_equal(l$stom_cond_CO2_ >0, TRUE)
  #transpiration becomes greater than 0 
  expect_equal(l$transpiration_ >0, TRUE)
  
  #when psi stem is equal to psi soil
  psi_soil = 2
  
   l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l$set_leaf_states_rates_from_psi_stem(psi_soil, psi_soil)
  
  
  #calculate the hydraulic cost usign the sperry method, should be 0 when psi_soil is equivalent to psi_stem
  #calculate hydraulic cost using TF method, should be greater than 0 when psi_soil is greater than 0
  expect_equal(l$hydraulic_cost_TF(psi_soil) > 0, TRUE)
  
  expect_equal(l$hydraulic_cost_TF(psi_soil + 1) > 0, TRUE)
  
  #ensure that hydraulic_cost returns 0 cost at 0 psi_soil/psi_stem
  psi_soil = 0
   l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
   
  l$set_leaf_states_rates_from_psi_stem(0, 0)
  expect_equal(l$hydraulic_cost_TF(psi_soil), 0)
  
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
   l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
   
  l$set_leaf_states_rates_from_psi_stem(psi_crit, psi_soil)
  c_i = l$ci_  
  benefit_ = l$assim_colimited_
  
  # Keep this R replica in step with phylloptim's constant: 0.6.0 derives it as
  # 1/molar_mass_h2o (55.509298) instead of the old hard-coded 55.4939. Leaving
  # the old value here would bias the round-trip below by exactly that 2.77e-4,
  # which is the whole content of the difference (#622).
  kg_to_mol_h2o = 1 / 0.018015
  umol_to_mol = 1e-6
  kPa_to_Pa = 1e3
  
  g_c_ci = ((benefit_)* umol_to_mol * l$atm_kpa_ * kPa_to_Pa)/(l$ca_ - l$ci_); 
  
  E_ci = g_c_ci * 1.67 * l$atm_vpd_ / kg_to_mol_h2o / l$atm_kpa_;
  # psi_upstream is a positive magnitude, same as everywhere else (phylloptim #25):
  # transpiration_to_psi_stem no longer negates it internally, so there is nothing
  # to invert here either.
  psi_stem = l$transpiration_to_psi_stem(E_ci, psi_soil)
  
  #conversion back and forth is not perfect
  # Compare against what was fed in at set_leaf_states_rates_from_psi_stem above.
  expect_equal(psi_stem, psi_crit, tolerance = 1e-05)
  
  #let's start testing profit functions
  
  #first off- what happens when we moving psi_soil around
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
   
    psi_soil = 0
    soil_depth = 0.5

  #first off- what happens when we moving psi_soil around
  # start with one soil layer. For the TF method, it will fail when there is more than one layer
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_crit + 1, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
    l$set_model("TF24", "stem"); l$optimise()
  
  expect_equal(l$transpiration_, 0)
  expect_equal(l$opt_psi_stem_, psi_crit+1)
  expect_equal(l$stom_cond_CO2_, 0)
  expect_equal(l$ci_, l$gamma_*0.1013)
  expect_true(l$profit_ < 0)
  
  # check a more standard case, wet soil

  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$set_model("TF24", "stem"); l$optimise()
  expect_true(l$profit_ > 0)
  expect_true(l$opt_psi_stem_ > 0)
  expect_true(l$stom_cond_CO2_ > 0)
  expect_true(l$transpiration_ > 0)
  expect_true(l$hydraulic_cost_ > 0)
  profit1 = l$profit_

  expect_true(is.na(l$profit_psi_stem_TF(NA, psi_soil)))

  #confirm different nubmer of layers fails for single layer methods
  psi_soil = c(0, 1)
  soil_depth = c(0.5, 1)
  root_carbon = c(0.5, 1)

  l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  expect_error({l$set_model("TF24", "stem"); l$optimise()}, "psi soil must have only one value to use non-root-based profit optimisation methods")

  #test various responses to environmental gradients to check that behaviour is being conserved
  
psi_soil = 1
soil_depth = 1
  #light
   l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
    l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  high_light <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 100, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  l$set_model("TF24", "stem"); l$optimise()
  
  low_light <- l$profit_
  
  expect_true(high_light > low_light)
  
  #soil moist
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  high_moist <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 2, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  low_moist <- l$profit_
  
  expect_true(high_moist > low_moist)
  
  
  #vpd'
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  low_vpd <- l$profit_
  
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 2, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  high_vpd <- l$profit_
  
  expect_true(high_vpd < low_vpd)
  
  #vcmax_25
  
  l <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  low_vcmax <- l$profit_
  
  l <- Leaf(vcmax_25 = 150, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)  
  l$set_model("TF24", "stem"); l$optimise()
  
  high_vcmax <- l$profit_
  
  expect_true(high_vcmax > low_vcmax)
  
  #test effect of leaf temperature
  
  l_low_temp <-Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
                            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l_low_temp$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = 20, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l_ref_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
                             beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                             ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l_ref_temp$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = 25, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  l_high_temp <- Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
                              beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                              GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                              ci_niter = ci_niter, g1_TF24 = g1_TF24)
  l_high_temp$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = 30, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  expect_true(l_low_temp$gamma_ < l_ref_temp$gamma_ &  l_ref_temp$gamma_ <  l_high_temp$gamma_)
  
    l_high_temp <-Leaf(vcmax_25 = 50, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
                             beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim,
                             GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
                             ci_niter = ci_niter, g1_TF24 = g1_TF24)
    l_high_temp$set_physiology(root_network = net((1) / area_leaf_, soil_depth), PPFD = 1000, psi_soil = 0, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = 1, ca = ca, leaf_temp = 40, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
  
  # Off the 25 C reference, phylloptim 0.6.0's derived gas constants shift these
  # in the 5th significant figure. plant runs its leaf at 25 C, where the
  # Arrhenius factors are 1 and nothing moves (#622).
  expect_equal(round(l_high_temp$ko_, 1), round(562292.4,1))
  expect_equal(round(l_high_temp$kc_, 1), round(1878.9,1))
  expect_equal(round(l_high_temp$gamma_, 1), round(88.800391,1))
  expect_equal(round(l_high_temp$vcmax_, 1), round(35.2,1))

  # test out the root component of the leaf unit

  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), 
            beta2= beta2, a = a, curv_fact_elec_trans = curv_fact_elec_trans, curv_fact_colim = curv_fact_colim, 
            GSS_tol_abs = GSS_tol_abs, vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol, 
            ci_niter = ci_niter, g1_TF24 = g1_TF24)
  
  # test that you get error when soil_depth and psi_soil have different number of layers

  soil_depth = c(0.5,1)
  psi_soil = c(0.5, 0.5)
  root_carbon = c(1, 1)

  l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

# root assignment should be equal to number of soil layers  
expect_equal(length(l$c_r_H_), length(soil_depth))
expect_equal(length(l$c_r_V_), length(soil_depth))

# test what happens when psi_root is equal to psi_soil, E should be slightly
# negative (the layer gains, because gravity still has to be paid). Collar and
# soil are both positive magnitudes now (phylloptim #25).
l$E_from_Soil_to_Root_Collar(psi_soil[1], psi_soil)
expect_true(l$E_up_ < 0)

  soil_depth = c(0.5)
  psi_soil = c(0.5)
  root_carbon = c(1)

# test what happens when psi_root is equal gravitational effect
l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
# the collar pulls exactly hard enough to lift the water and no harder, so E == 0
l$E_from_Soil_to_Root_Collar((psi_soil + l$z_soil_mid_*9.8e-3), psi_soil[1])
expect_equal(l$E_up_, 0)

# test what happens when psi_root is equal gravitational effect
l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
l$E_from_Soil_to_Root_Collar(psi_soil[1] + 0.5, psi_soil[1])
expect_true(l$E_up_ > 0)

  soil_depth = c(0.5,1)
  psi_soil = c(0.5, 0.5)
  root_carbon = c(1, 0)
l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
l$E_from_Soil_to_Root_Collar(psi_soil[1] + 0.5, psi_soil[1])

#confirm that soil_consumption is 0 when roots do not exist in that layer
expect_equal(l$soil_consumption_[2], 0)

  soil_depth = c(0.5,1)
  psi_soil = c(0.5, 0.5)
  root_carbon = c(1, 1)
l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
l$E_from_Soil_to_Root_Collar(psi_soil[1] + 0.5, psi_soil[1])

#check that soil consumption adds to E_up
expect_equal(l$E_up_, sum(l$soil_consumption_)*0.018015)

  soil_depth = c(0.5,1)
  psi_soil = c(1e6, 1e6)
  root_carbon = c(1, 1)
l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
l$find_root_collar_psi()
# R_d is a trait with its own default in phylloptim 0.6.0, not 0.015*vcmax_25
# (they coincide at vcmax_25 = 96), and psi_crit is derived. Read both off the
# leaf so this stays an identity rather than a restatement of old formulas.
expect_equal(l$profit_, -l$R_d_ - l$hydraulic_cost_TF(l$psi_crit))
l$opt_root_psi_

# assim_max_ < 0 early-exit: wet soil (so the upstream shut-down exits are NOT
# taken) but zero light, so maximum assimilation at ci = ca is below dark
# respiration. The plant operates at root_zero_E (collar where soil uptake is
# zero); at zero transpiration the stem equilibrates with the collar.
  soil_depth = c(0.5, 1)
  psi_soil = c(0, 0)
  root_carbon = c(1, 1)
l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = 0, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

# confirm we exercise the assim_max_ < 0 branch
expect_true(l$assim_max_ < 0)

l$find_root_collar_psi()

# opt_psi_stem_ is reported as a POSITIVE magnitude in every branch of
# find_root_collar_psi (the assim_max_ < 0 exit used to leak the signed
# root_zero_E here -- the sign wart fixed alongside review #7).
expect_true(l$opt_psi_stem_ > 0)
# opt_root_psi_ is a positive magnitude, like every other psi (phylloptim #25).
expect_true(l$opt_root_psi_ > 0)
# at zero transpiration the stem equilibrates with the collar: same potential, and
# now literally the same number rather than a magnitude paired with its negation.
expect_equal(l$opt_psi_stem_, l$opt_root_psi_)
})

# Medlyn stomatal-conductance model (ported/adapted from develop #450). The
# Medlyn solvers are a standalone, R-callable alternative to the root-collar
# profit optimisation and are not on the TF24 compute path; g0/g1 are exposed as
# settable fields rather than constructor args in this branch.
test_that("Medlyn stomatal model", {
  vcmax_25 = 100
  jmax_25 = vcmax_25 * 167
  c = 2.04
  b = 3
  psi_crit = 5
  theta = 0.000157
  K_s = 1
  h = 5
  beta2 = 1
  curv_fact_elec_trans = 0.7
  a = 0.3
  curv_fact_colim = 0.99
  g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8
  vulnerability_curve_ncontrol = 100
  ci_abs_tol = 1e-6
  ci_niter = 1000
  root_c = 2.65
  root_b = 1.29
  root_psi_crit = root_b * (log(1.0 / 0.05))^(1.0 / root_c)

  PPFD = 900
  leaf_specific_conductance_max = K_s * theta / h
  atm_vpd = 2
  ca = 40
  atm_o2_kpa_ = 21
  leaf_temp_ = 25
  atm_kpa_ = 101.3
  area_leaf_ = 0.05

  make_leaf <- function() {
    Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
         root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), beta2 = beta2,
         a = a, curv_fact_elec_trans = curv_fact_elec_trans,
         curv_fact_colim = curv_fact_colim, GSS_tol_abs = GSS_tol_abs,
         vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol,
         ci_niter = ci_niter, g1_TF24 = g1_TF24)
  }
  set_phys <- function(l, psi_soil = 2, atm_vpd = 2) {
    l$set_physiology(root_network = net((1) / area_leaf_, 0.5), PPFD = PPFD, psi_soil = psi_soil, soil_depth = 0.5, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
    l
  }

  # g0/g1 default to the published Medlyn (2011) values (settable fields here)
  l <- make_leaf()
  expect_equal(l$g0, 0.022)
  expect_equal(l$g1, 2.57)

  # numerical solver: reference values for this branch's build (regression guard)
  l <- set_phys(make_leaf())
  l$solve_medlyn_ci_numerical()
  # These three moved with phylloptim 0.6.0's kg_to_mol_h2o, which changed from
  # a hard-coded 55.4939 to 1/molar_mass_h2o = 55.509298 (+2.77e-4 relative).
  # It enters the E -> g_sc conversion, so g_sc carries exactly that relative
  # shift and ci and assimilation follow (#622).
  expect_equal(l$ci_, 19.6244425204, tolerance = 1e-5)
  expect_equal(l$stom_cond_CO2_, 0.120848522773, tolerance = 1e-6)
  expect_equal(l$assim_colimited_, 15.1922296035, tolerance = 1e-5)
  # operating point is physically sane: gamma* < ci < ca, positive gs and assim
  expect_true(l$ci_ > 0 && l$ci_ < ca)
  expect_true(l$stom_cond_CO2_ > 0)
  expect_true(l$assim_colimited_ > 0)

  # psi_soil has no effect on the Medlyn solution (it is a soil-water, not a
  # supply-side, model): same ci/gs/assim for a different soil potential
  l2 <- set_phys(make_leaf(), psi_soil = 3)
  l2$solve_medlyn_ci_numerical()
  expect_equal(l2$ci_, l$ci_, tolerance = 1e-5)
  expect_equal(l2$stom_cond_CO2_, l$stom_cond_CO2_, tolerance = 1e-6)

  # analytical version: stomatal conductance decreases with rising vapour-pressure
  # deficit (Medlyn 2011 1/sqrt(D) sensitivity)
  D <- seq(0.5, 5, 0.25)
  gs <- vapply(D, function(d) {
    li <- set_phys(make_leaf(), atm_vpd = d)
    li$g1 <- 3.3
    li$solve_medlyn_ci_analytical()
    li$stom_cond_CO2_
  }, numeric(1))
  expect_true(coef(stats::lm(gs ~ D))[[2]] < 0)

  # at field-capacity soil moisture (beta_ == 1) with zero residual conductance
  # (g0 == 0), the numerical optimisation and the analytical Medlyn solution must
  # coincide -- a core correctness check of the coupled solver.
  ln <- set_phys(make_leaf())
  ln$g0 <- 0; ln$g1 <- 3.3; ln$theta_ <- ln$theta_fc_
  ln$solve_medlyn_ci_numerical()
  numerical_ci <- ln$ci_

  la <- set_phys(make_leaf())
  la$g0 <- 0; la$g1 <- 3.3; la$theta_ <- la$theta_fc_
  la$solve_medlyn_ci_analytical()
  analytical_ci <- la$ci_

  expect_equal(numerical_ci, analytical_ci, tolerance = 1e-5)
})

# psi_stem_to_ci is the TF24 compute-path inner solve and the single hottest
# function in the root-water profiling (#486): given a stem water potential it
# fixes the supply-side stomatal conductance gc and returns the ci at which
# CO2 demand equals supply, i.e. the root of
#   assim_colimited(ci) * umol_to_mol  -  gc * (ca - ci) / (atm_kpa * kPa_to_Pa) = 0
# over ci in (gamma*, ca]. It is currently solved by bisection (util::uniroot)
# to an absolute tolerance of 1e-7 Pa. These tests pin the *contract* (what ci
# must satisfy) independently of the *method*, so they remain a valid safety net
# if the solver is later refactored for speed.
test_that("psi_stem_to_ci supply=demand solve", {
  vcmax_25 = 100
  jmax_25 = vcmax_25 * 167
  c = 2.04
  b = 3
  psi_crit = 5
  theta = 0.000157
  K_s = 1
  h = 5
  beta2 = 1
  curv_fact_elec_trans = 0.7
  a = 0.3
  curv_fact_colim = 0.99
  g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8
  vulnerability_curve_ncontrol = 100
  ci_abs_tol = 1e-6
  ci_niter = 1000
  root_c = 2.65
  root_b = 1.29
  root_psi_crit = root_b * (log(1.0 / 0.05))^(1.0 / root_c)

  PPFD = 900
  leaf_specific_conductance_max = K_s * theta / h
  atm_vpd = 2
  ca = 40
  atm_o2_kpa_ = 21
  leaf_temp_ = 25
  atm_kpa_ = 101.3
  area_leaf_ = 0.05

  # unit conversions used inside psi_stem_to_ci (see the leaf package's
  # leaf_model.hpp). umol_per_mol_to_Pa is derived from the pressure, not a
  # constant: the old hard-coded 0.1013 was 101.3 kPa in disguise, and the leaf
  # package now computes it per call (phylloptim #15 item 10c). Deriving it here
  # keeps the test honest if atm_kpa_ ever moves off 101.3 -- as it is on the
  # TF24 driver default of 100.5.
  umol_to_mol = 1e-6
  kPa_to_Pa = 1e3
  umol_per_mol_to_Pa = atm_kpa_ * kPa_to_Pa * umol_to_mol

  make_leaf <- function() {
    Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
         root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), beta2 = beta2,
         a = a, curv_fact_elec_trans = curv_fact_elec_trans,
         curv_fact_colim = curv_fact_colim, GSS_tol_abs = GSS_tol_abs,
         vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol,
         ci_niter = ci_niter, g1_TF24 = g1_TF24)
  }
  set_phys <- function(l, psi_soil = 2, atm_vpd = 2) {
    l$set_physiology(root_network = net((1) / area_leaf_, 0.5), PPFD = PPFD, psi_soil = psi_soil, soil_depth = 0.5, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
    l
  }

  psi_soil <- 2

  # --- 1. the returned ci is the root of the supply=demand residual ----------
  # assim_minus_stom_cond_CO2 recomputes gc from (psi_stem, psi_upstream) using
  # exactly the same expression psi_stem_to_ci solves, so the residual at the
  # returned ci must be ~0. The bisection tol is 1e-7 Pa in ci and the residual
  # slope is O(1e-6) per Pa, so the residual is well below 1e-9.
  l <- set_phys(make_leaf(), psi_soil = psi_soil)
  psi_stem <- psi_soil + 1
  ci <- l$psi_stem_to_ci(psi_stem, psi_soil)
  expect_true(abs(l$assim_minus_stom_cond_CO2(ci, psi_stem, psi_soil)) < 1e-9)

  # --- 2. ci lies strictly inside the bracket (gamma*, ca) -------------------
  gamma_star <- l$gamma_ * umol_per_mol_to_Pa
  expect_true(ci > gamma_star && ci < ca)

  # --- 3. the call writes the ci_ member (used by the compute path) ----------
  expect_equal(l$ci_, ci)

  # --- 4. method-independent reference: re-solve the identical target in pure
  # R with a high-accuracy root finder. This decouples "the correct ci" from
  # "the C++ solver", so a future method swap (Newton, Brent, analytic, ...)
  # still passes as long as it converges to the same root.
  gc_fixed <- l$stom_cond_CO2(psi_stem, psi_soil)
  target <- function(x) {
    l$assim_colimited(x) * umol_to_mol - gc_fixed * (ca - x) / (atm_kpa_ * kPa_to_Pa)
  }
  ci_ref <- stats::uniroot(target, lower = gamma_star, upper = ca, tol = 1e-10)$root
  expect_equal(ci, ci_ref, tolerance = 1e-6)

  # --- 5. monotonicity: a higher (more negative-gradient) psi_stem raises gc,
  # steepening the supply line, so the supply=demand ci increases toward ca.
  psi_stem_seq <- seq(psi_soil + 0.1, psi_crit, length.out = 25)
  ci_seq <- vapply(psi_stem_seq, function(p) l$psi_stem_to_ci(p, psi_soil), numeric(1))
  expect_true(all(diff(ci_seq) > 0))
  expect_true(all(ci_seq > gamma_star & ci_seq < ca))

  # --- 6. the identity holds across a range of soil potentials ---------------
  for (ps in c(0, 0.5, 1, 3)) {
    lp <- set_phys(make_leaf(), psi_soil = ps)
    pstem <- ps + 0.75
    ci_p <- lp$psi_stem_to_ci(pstem, ps)
    expect_true(ci_p > lp$gamma_ * umol_per_mol_to_Pa && ci_p < ca)
    expect_true(abs(lp$assim_minus_stom_cond_CO2(ci_p, pstem, ps)) < 1e-9)
  }

  # --- 7. degenerate gc=0 limit (psi_stem == psi_upstream): supply term
  # vanishes and the solve collapses to the light-compensation point where
  # net assimilation is zero (assim_colimited(ci) == 0).
  l0 <- set_phys(make_leaf(), psi_soil = psi_soil)
  ci0 <- l0$psi_stem_to_ci(psi_soil, psi_soil)
  expect_true(abs(l0$assim_colimited(ci0)) < 1e-5)

  # --- 8. regression guard: reference ci for the standard scenario on this
  # build. A solver change that alters the converged value beyond rounding is
  # expected to update this number (it is NOT bit-identical across methods).
  expect_equal(ci, 11.7559822995, tolerance = 1e-6)
})

# find_root_psi is the inner soil->root-collar continuity solve (#486). Given a
# bracket of candidate (signed, negative) root-collar potentials it returns the
# collar potential x where one of two hydraulic targets vanishes:
#   find_root_crit == 0 (root_zero_E): total soil uptake E_up_(x) == 0
#   find_root_crit == 1 (root_crit)  : E_up_(x) == stem demand at psi_crit, i.e.
#                                       E_column(x) = E_up_(x) - transpiration(psi_crit, -x) == 0
# over x in [-psi_crit, wettest_soil_layer]. It is currently solved by bisection
# (util::uniroot) to an absolute tolerance of 1e-4 in x, and sits inside the
# golden-section collar optimisation -- so it dominates the per-layer
# E_from_Soil_to_Root_Collar arithmetic that is the top hot-spot in the
# 15-layer driver. These tests pin the *contract* (what x must satisfy)
# independently of the *method*, so they remain a valid safety net if the solver
# is later swapped for a faster root finder (Brent/TOMS748, ...).
test_that("find_root_psi soil->collar continuity solve", {
  vcmax_25 = 100
  jmax_25 = vcmax_25 * 167
  c = 2.04
  b = 3
  psi_crit = 5
  theta = 0.000157
  K_s = 1
  h = 5
  beta2 = 1
  curv_fact_elec_trans = 0.7
  a = 0.3
  curv_fact_colim = 0.99
  g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8
  vulnerability_curve_ncontrol = 100
  ci_abs_tol = 1e-6
  ci_niter = 1000
  root_c = 2.65
  root_b = 1.29
  root_psi_crit = root_b * (log(1.0 / 0.05))^(1.0 / root_c)

  PPFD = 900
  leaf_specific_conductance_max = K_s * theta / h
  atm_vpd = 2
  ca = 40
  atm_o2_kpa_ = 21
  leaf_temp_ = 25
  atm_kpa_ = 101.3
  area_leaf_ = 0.05

  make_leaf <- function() {
    Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
         root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), beta2 = beta2,
         a = a, curv_fact_elec_trans = curv_fact_elec_trans,
         curv_fact_colim = curv_fact_colim, GSS_tol_abs = GSS_tol_abs,
         vulnerability_curve_ncontrol = vulnerability_curve_ncontrol, ci_abs_tol = ci_abs_tol,
         ci_niter = ci_niter, g1_TF24 = g1_TF24)
  }

  # 15-layer soil column with a mild moisture gradient (wettest at the surface),
  # mirroring the validation driver (env$set_soil_number_of_depths(15)). Distinct
  # per-layer psi_soil exercises the layer-by-layer branch switches in
  # E_from_Soil_to_Root_Collar (each x == psi_soil[i] crossing is a candidate
  # kink in the target -- the smoothness question the solver-swap hinges on).
  n_layer <- 15
  psi_soil <- seq(0.3, 0.7, length.out = n_layer)     # positive magnitudes
  soil_depth <- seq(0.1, 1.5, length.out = n_layer)   # cumulative layer depths (m)
  root_carbon <- rep(1, n_layer)

  set_phys <- function(l) {
    l$set_physiology(root_network = net((root_carbon) / area_leaf_, soil_depth), PPFD = PPFD, psi_soil = psi_soil, soil_depth = soil_depth, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)
    l
  }

  l <- set_phys(make_leaf())

  # Reconstruct the bracket exactly as find_root_collar_psi does. Every psi is a
  # positive magnitude (phylloptim #25), so nothing is flipped and the wettest layer
  # is the SMALLEST suction -- which makes it the LOWER bracket end, where the
  # signed convention had it as the upper one. The bracket runs from there to
  # psi_crit, the driest feasible collar.
  wettest <- min(psi_soil)
  lower <- wettest
  # phylloptim 0.6.0 derives psi_crit from (P50, c) rather than taking it as a
  # constructor argument, so read it off the leaf -- the block's local psi_crit
  # is no longer what the leaf is using (#622).
  psi_crit <- l$psi_crit
  upper <- psi_crit

  # Method-independent scalar targets, rebuilt in pure R from the R-exposed
  # primitives (E_from_Soil_to_Root_Collar sets E_up_; transpiration is direct).
  # These are bit-for-bit the C++ E_column_zero / E_column the solver drives.
  target0 <- function(x) {           # find_root_crit == 0
    l$E_from_Soil_to_Root_Collar(x, psi_soil)
    l$E_up_
  }
  target1 <- function(x) {           # find_root_crit == 1
    l$E_from_Soil_to_Root_Collar(x, psi_soil)
    E_up <- l$E_up_
    # E_column demands transpiration(psi_crit, x): x is the collar suction, the
    # same kind of number as psi_crit (phylloptim #25), so nothing is flipped.
    E_up - l$transpiration(psi_crit, x)
  }

  # --- 1. both targets are bracketed (opposite signs at the endpoints) --------
  expect_true(target0(lower) < 0 && target0(upper) > 0)
  expect_true(target1(lower) < 0 && target1(upper) > 0)

  # --- 2. both targets are strictly monotone over the bracket -----------------
  # (a clean single sign-change, the necessary condition for ANY bracketing
  # solver -- bisection or superlinear -- to be safe here.)
  xs <- seq(lower, upper, length.out = 400)
  t0 <- vapply(xs, target0, numeric(1))
  t1 <- vapply(xs, target1, numeric(1))
  expect_true(all(is.finite(t0)) && all(is.finite(t1)))
  # Both targets now INCREASE with x: a larger collar suction pulls more water up.
  # The sign of the slope flipped with the variable, and so did the sign pattern at
  # the endpoints -- a bracketing solver only needs them opposite, not in a
  # particular order.
  expect_true(all(diff(t0) > 0))
  expect_true(all(diff(t1) > 0))

  # --- 3. C++ root matches a tight method-independent R reference -------------
  # stats::uniroot solves the identical target to 1e-12; the C++ solver works to
  # 1e-4 in x, so agreement to ~1e-3 pins "the correct collar potential"
  # independently of the C++ method (a future Brent/TOMS748 swap still passes).
  root0 <- l$find_root_psi(wettest, psi_soil, 0L)
  root1 <- l$find_root_psi(wettest, psi_soil, 1L)
  ref0 <- stats::uniroot(target0, lower = lower, upper = upper, tol = 1e-12)$root
  ref1 <- stats::uniroot(target1, lower = lower, upper = upper, tol = 1e-12)$root
  expect_equal(root0, ref0, tolerance = 1e-3)
  expect_equal(root1, ref1, tolerance = 1e-3)

  # --- 4. residual at the returned root is ~0 ---------------------------------
  # The target slope is O(1e-4 kg/MPa); at a 1e-4 root tol the residual is well
  # below 1e-7. Also confirm the root is strictly interior to the bracket.
  expect_true(abs(target0(root0)) < 1e-7)
  expect_true(abs(target1(root1)) < 1e-7)
  expect_true(root0 > lower && root0 < upper)
  expect_true(root1 > lower && root1 < upper)

  # --- 5. ordering: the zero-uptake collar is wetter -- a SMALLER suction -- than
  # the critical-demand collar (more water is drawn at psi_crit than at zero flux).
  # The inequality reverses with the variable; the physics does not.
  expect_true(root0 < root1)

  # --- 6. find_psi_stem_from_psi_root contract --------------------------------
  # Given a collar suction it returns a finite stem suction >= it (the stem is
  # downstream, hence under at least as much tension), and the mapping is monotone
  # increasing in collar dryness. Both sides are the same kind of number now, so
  # the comparison needs no negation.
  roots <- seq(root0, root1, length.out = 8)
  psi_stems <- vapply(roots, function(r) l$find_psi_stem_from_psi_root(r, psi_soil), numeric(1))
  expect_true(all(is.finite(psi_stems)))
  # The tolerance is the continuity tol (1e-4 in x): at the zero-uptake collar
  # (root0) the flux -> 0 so psi_stem -> collar and rounding can place it
  # microscopically either side -- a 1e-8 bound would be method-dependent there.
  expect_true(all(psi_stems >= roots - 1e-3))
  expect_true(all(diff(psi_stems) > 0))          # drier collar -> drier stem

  # --- 7. NaN-input propagation: a non-finite soil potential must fail fast
  # (util::stop in E_from_Soil_to_Root_Collar), NOT return a silent NaN root.
  psi_soil_bad <- psi_soil
  psi_soil_bad[5] <- NA_real_
  expect_error(l$find_root_psi(wettest, psi_soil_bad, 0L))

  # --- 7b. and so must a SIGNED vector, which is what a pre-#25 caller has. It
  # would otherwise compile, run, and return a wrong number in silence.
  expect_error(l$find_root_psi(wettest, -psi_soil, 0L))
  expect_error(l$E_from_Soil_to_Root_Collar(2.0, -psi_soil))

  # --- 8. regression guard: hardcoded reference roots for the standard scenario.
  # These pin the physical answer independently of the build/method; both the
  # bisection and a superlinear bracketing solver converge here to the same
  # collar potentials within the 1e-4 continuity tolerance (the value is NOT
  # bit-identical across methods, hence the loose tolerance).
  # Sign flipped with the representation; the magnitudes are unchanged, which is
  # the point of #25 being a representation change and not a model change.
  expect_equal(root0, 0.4387473787, tolerance = 1e-3)
  expect_equal(root1, 0.686287940458, tolerance = 1e-3)
})

test_that("dprofit_droot_collar_psi matches a finite difference (AD/IFT gradient, #527)", {
  # Reuse the standard single-layer leaf setup from "Basic functions".
  vcmax_25 = 100; jmax_25 = vcmax_25 * 167; c = 2.04; b = 3; psi_crit = 5
  theta = 0.000157; K_s = 1; h = 5; beta2 = 1;   curv_fact_elec_trans = 0.7; a = 0.3; curv_fact_colim = 0.99; g1_TF24 = 46.32995
  GSS_tol_abs = 1e-8; vulnerability_curve_ncontrol = 100; ci_abs_tol = 1e-6
  ci_niter = 1000; root_c = 2.65; root_b = 1.29
  root_psi_crit = root_b * (log(1.0 / 0.05))^(1.0 / root_c)
  l <- Leaf(vcmax_25 = vcmax_25, jmax_25 = jmax_25, c = c, p_50 = b * (log(2))^(1 / c),
            root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c), beta2 = beta2,
            a = a, curv_fact_elec_trans = curv_fact_elec_trans,
            curv_fact_colim = curv_fact_colim, GSS_tol_abs = GSS_tol_abs,
            vulnerability_curve_ncontrol = vulnerability_curve_ncontrol,
            ci_abs_tol = ci_abs_tol, ci_niter = ci_niter, g1_TF24 = g1_TF24)
  PPFD = 900
  leaf_specific_conductance_max = K_s * theta / h
  psi_soil = 2; atm_vpd = 2; ca = 40; atm_o2_kpa_ = 21; leaf_temp_ = 25
  atm_kpa_ = 101.3; area_leaf_ = 0.05; root_carbon_ = 1
  # psi_soil was 2.0 until phylloptim #24. This leaf pairs a fragile root system
  # (root_b = 1.29, so root_psi_crit = 1.952) with a tough stem (psi_crit = 5), and
  # at psi_soil = 2.0 the soil is drier than the roots can tolerate: root_zero_E =
  # 2.005 sits BELOW root_psi_crit, so there is no operating point that both moves
  # water and stays inside the root limit, and the leaf now shuts down. The test was
  # relying on the dead clamp to let it operate anyway. 1.25 leaves the interval
  # [1.255, 1.952], which is genuinely feasible.
  psi_soil = 1.25
  l$set_physiology(root_network = net((root_carbon_) / area_leaf_, 1), PPFD = PPFD, psi_soil = psi_soil, soil_depth = 1, leaf_specific_conductance_max = leaf_specific_conductance_max, atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp_, atm_o2_kpa = atm_o2_kpa_, atm_kpa = atm_kpa_)

  # Optimise once (also sets up the soil-side caches the gradient needs).
  l$find_root_collar_psi()
  opt <- l$opt_root_psi_   # operating collar suction (positive magnitude)

  fd_grad <- function(psi, eps = 1e-5) {
    (l$evaluate_root_collar_psi(psi + eps) -
     l$evaluate_root_collar_psi(psi - eps)) / (2 * eps)
  }

  # On the feasible interior the exact AD/IFT gradient must match a fine central
  # finite difference of profit. evaluate_root_collar_psi clamps to the feasible
  # interval, so a clamped point would make the central FD straddle the boundary
  # -- skip those and test strictly-interior points.
  #
  # The points deliberately extend well past the optimum (out to opt + 0.6).
  tested <- 0
  for (psi in c(opt + 0.02, opt + 0.1, opt + 0.2, opt + 0.4, opt + 0.6)) {
    l$evaluate_root_collar_psi(psi)
    used <- l$opt_root_psi_
    if (abs(used - psi) > 1e-8) next            # clamped: not interior, skip
    ad <- l$dprofit_droot_collar_psi(psi)
    expect_true(is.finite(ad))
    # 1e-5 step here: the exact gradient tracks a fine FD to ~1e-9 across the
    # interior; 1e-6 tolerance leaves margin without admitting the 0.2%
    # extrapolation-region bug above.
    expect_equal(ad, fd_grad(psi), tolerance = 1e-6)
    tested <- tested + 1
  }
  expect_gt(tested, 2)                          # several interior points

  # --- the extrapolation region, which the clamp now puts out of the solver's reach
  #
  # The analytic dpsi_stem/dpsi term differentiates the cumulative root-vulnerability
  # spline, which CLAMP-extrapolates beyond its knot domain. Using the separate f_r
  # spline as the integrand (rather than the integral spline's own .deriv) was correct
  # on-domain and silently wrong beyond it (~0.2% at 0.4 MPa past the optimum), so
  # this property is load-bearing and must keep being tested.
  #
  # It used to be tested by the far points above. It cannot be any more: the root
  # spline's knots end at the 1%-conductivity potential (2.295 MPa here) while
  # phylloptim #24 clamps the collar to root_psi_crit, the 5% point (1.952) -- which is
  # ALWAYS the smaller. So the extrapolation region is now unreachable through
  # evaluate_root_collar_psi by construction, and a test that goes through the clamp
  # can only skip it.
  #
  # So difference the profit directly instead, on the unclamped path
  # dprofit_droot_collar_psi itself uses: psi -> psi_stem -> profit. Same identity,
  # no clamp, and it reaches past 2.295 where the two splines disagree.
  profit_at <- function(psi) {
    l$profit_psi_stem_TF(l$find_psi_stem_from_psi_root(psi, psi_soil), psi)
  }
  fd_unclamped <- function(psi, eps = 1e-5) {
    (profit_at(psi + eps) - profit_at(psi - eps)) / (2 * eps)
  }
  root_spline_end <- root_b * log(1 / 0.01)^(1 / root_c)   # 1% conductivity, 2.295
  extrapolated <- 0
  for (psi in c(root_spline_end + 0.1, root_spline_end + 0.3, root_spline_end + 0.5)) {
    ad <- l$dprofit_droot_collar_psi(psi)
    expect_true(is.finite(ad))
    expect_equal(ad, fd_unclamped(psi), tolerance = 1e-6)
    extrapolated <- extrapolated + 1
  }
  expect_gt(extrapolated, 2)   # all three, beyond the root spline's knot domain
})

test_that("Leaf() errors on misspelled argument names (issue #377)", {
  # R's partial matching silently accepts abbreviations; the wrapper should
  # catch any name that does not exactly match a declared parameter.
  vcmax_25 <- 100; jmax_25 <- 100 * 167; c <- 2.04; b <- 3; psi_crit <- 5
  root_c <- 2.65; root_b <- 1.29
  root_psi_crit <- root_b * (log(1 / 0.05))^(1 / root_c)
  common_args <- list(
    c = c, p_50 = b * (log(2))^(1 / c), root_c = root_c, root_p_50 = root_b * (log(2))^(1 / root_c),
    beta2 = 1, jmax_25 = jmax_25,
    a = 0.3, curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
    GSS_tol_abs = 1e-8, vulnerability_curve_ncontrol = 100,
    ci_abs_tol = 1e-6, ci_niter = 1000, g1_TF24 = 46.33
  )

  # exact name is accepted
  expect_silent(do.call(Leaf, c(list(vcmax_25 = vcmax_25), common_args)))

  # misspelled name should error
  expect_error(
    do.call(Leaf, c(list(vcma = vcmax_25), common_args)),
    regexp = "vcma"
  )
})
