# Extracted from test-strategy-tf24.R:812

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

# test -------------------------------------------------------------------------
s <- TF24_Strategy(collect_all_auxiliary = TRUE)
s$pars$a_pl0 <- 1.0
wet <- Environment("TF24")
wet$set_soil_number_of_depths(5)
wet$set_soil_water_state(rep(0.30, 5))
wet$set_fixed_environment(1.0, 40)
dry <- Environment("TF24")
dry$set_soil_number_of_depths(5)
dry$set_soil_water_state(rep(0.11, 5))
dry$set_fixed_environment(1.0, 40)
ind <- TF24_Individual(s)
ind$set_state("height", 5.0)
ind$set_initial_states(wet)
ind$compute_rates(dry)
expect_lt(ind$aux("leaf_marginal_return", "sapwood_marginal_return"), 0)
