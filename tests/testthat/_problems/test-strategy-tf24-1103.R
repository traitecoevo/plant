# Extracted from test-strategy-tf24.R:1103

# prequel ----------------------------------------------------------------------
tf24_storage_capacity_r <- function(height, s = TF24_Strategy()) {
  z <- rep(0, length(height))
  s$pars$a_st1 * TF24_strategy_expand_allometry(s, height, z, z)$mass_sapwood
}
tf24_storage_at <- function(storage, height = 5, light = 1, theta = 0.4,
                            s = TF24_Strategy()) {
  env <- Environment("TF24")
  env$set_fixed_environment(light, height_max = 150)
  env$set_soil_water_state(rep(theta, env$get_soil_number_of_depths()))
  env$time <- 5
  ind <- Individual("TF24", "TF24_Env")(s)
  ind$set_state("height", height)
  ind$set_state("storage", storage)
  ind$compute_rates(env)
  list(dS = ind$internals$rates[[match("storage", ind$ode_names)]],
       P = ind$aux("net_mass_production_dt"))
}
tf24_departure_env <- function() {
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.2, 5))
  env$set_fixed_environment(1.0, 40)
  env
}
tf24_sapwood_probe <- function(psi, theta_soil, a_sw = 1.0, height = 10) {
  s <- TF24_Strategy(collect_all_auxiliary = TRUE)
  s$pars$a_sw <- a_sw
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5L)
  env$set_soil_water_state(rep(theta_soil, 5L))
  env$set_fixed_environment(1.0, 40)
  ind <- TF24_Individual(s)
  ind$set_state("height", height)
  ind$set_state("log_area_sapwood_departure", psi)
  ind$set_initial_states(env)
  ind$compute_rates(env)
  list(R_s = ind$aux("sapwood_marginal_return"),
       dpsi = ind$rate("log_area_sapwood_departure"),
       dheight = ind$rate("height"))
}

# test -------------------------------------------------------------------------
s <- TF24_Strategy()
s$pars$a_sw <- 1.0
s$control$shading_model <- "deep_crown"
expect_error(TF24_Individual(s), "deep-crown")
