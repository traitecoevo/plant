# Extracted from test-strategy-tf24.R:769

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
probe <- function(theta_soil) {
    env <- Environment("TF24")
    env$set_soil_number_of_depths(5)
    env$set_soil_water_state(rep(theta_soil, 5))
    env$set_fixed_environment(1.0, 40)
    ind <- TF24_Individual(s)
    ind$set_state("height", 5.0)
    ind$set_initial_states(env)          # reserves at a_st3 = 0.8 of capacity
    ind$compute_rates(env)
    list(dphi = ind$rate("log_area_leaf_departure"),
         dpsi = ind$rate("log_area_sapwood_departure"),
         margin = ind$aux("leaf_marginal_return", "sapwood_marginal_return"))
  }
