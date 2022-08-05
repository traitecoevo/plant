rm(list=ls())
devtools::load_all(".")

p0 <- scm_base_parameters("FF16w")
p0$disturbance_mean_interval <- 2
p0$max_patch_lifetime <- 5

p50_1 <- 2
p50_2 <- 4
p50_3 <- 5

p1 <- expand_parameters(trait_matrix(c(p50_1,p50_2,p50_3), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1,1,1))
env <- make_environment("FF16w", soil_initial_state = rep(0.2, 1), rainfall = 0.5)


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


p1$strategies[[1]]$K_s <- p_50_2_K_s(p50_1)
p1$strategies[[2]]$K_s <- p_50_2_K_s(p50_2)
p1$strategies[[3]]$K_s <- p_50_2_K_s(p50_3)

p1$strategies[[1]]$b <- calc_vul_b(p50_1, p1$strategies[[1]]$c)
p1$strategies[[2]]$b <- calc_vul_b(p50_2, p1$strategies[[2]]$c)
p1$strategies[[3]]$b <- calc_vul_b(p50_3, p1$strategies[[3]]$c)

p1$strategies[[1]]$psi_crit <- calc_psi_crit(p1$strategies[[1]]$b, p1$strategies[[1]]$c)
p1$strategies[[2]]$psi_crit <- calc_psi_crit(p1$strategies[[2]]$b, p1$strategies[[2]]$c)
p1$strategies[[3]]$psi_crit <- calc_psi_crit(p1$strategies[[3]]$b, p1$strategies[[3]]$c)

# 
# test_low5 <- c(1,2)
# saveRDS(test_low5, "plant/test_low5.RDS")


ctrl = scm_base_control()
# # 
# ctrl$ode_tol_abs <- ctrl$ode_tol_rel <- 1e-3
# ctrl$ode_step_size_initial <- ctrl$ode_step_size_min <- 1e-5
# ctrl$schedule_eps <- 0.01

system.time(result <- run_scm_collect(p1, env, ctrl))
saveRDS(result, "plant/result_low_2.5.RDS")
ctrl$schedule_eps


results_tidy <- 
  result2 %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded

data_species_tot_med <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_med %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()

results_tidy$species %>%
  plot_size_distribution()

results_tidy$time



