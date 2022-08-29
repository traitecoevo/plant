library(purrr)
library(tidyverse)

devtools::load_all(".")

p_50_2_K_s <- function(p50){
  1.638999*(p50/2)^(-1.38)
}

calc_vul_b <- function(p_50, c){
  num <- p_50
  den <- (-log(1-50/100))^(1/c)
  num/den
}


calc_psi_crit <- function(b,c) {
  b*(log(1/0.05))^(1/c)
}

run_FF16w_model <- function(inputs){
  
  
  ctrl = scm_base_control()
  
  p0 <- scm_base_parameters("FF16w")
  p0$max_patch_lifetime <- inputs$max_patch_lifetime
  p50 <- inputs$p50
  
  p1 <- expand_parameters(trait_matrix(c(p50), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1))
  
  p1$strategies[[1]]$K_s <- p_50_2_K_s(p50)
  
  p1$strategies[[1]]$b <- calc_vul_b(p50, p1$strategies[[1]]$c)
  
  p1$strategies[[1]]$psi_crit <- calc_psi_crit(p1$strategies[[1]]$b, p1$strategies[[1]]$c)

  p1$strategies[[1]]$epsilon_leaf <- 0.00001
  
  env <- make_environment("FF16w", soil_initial_state = rep(0.5, 1), rainfall = inputs$rainfall)
  system.time(result_low <- build_schedule(p1, env, ctrl))
  system.time(result_low <- run_scm_collect(p1, env, ctrl))
  
  result_low
}

posslm1 = safely(.f = run_FF16w_model, otherwise = "Error")

expand_grid(rainfall = seq(0.1, 3, length.out = 2), p50 = seq(0.1, 12, length.out = 2), max_patch_lifetime = 10) %>%
  split(1:nrow(.)) %>%
  map(., ~posslm1(.)) -> results
