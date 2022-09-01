# 
# test_start <- c(1,2)
# saveRDS(test_start, "plant/test_start.RDS")

install.packages("devtools")
library(devtools)
# 
# test_low3 <- c(1,2)
# saveRDS(test_low3, "plant/test_low3.RDS")


packages <- c("Rcpp", "R6", "crayon", "nleqslv", "BB" ,"BH")
install.packages(packages)
# 
# test_low2 <- c(1,2)
# saveRDS(test_low2, "plant/test_low2.RDS")

devtools::install_deps("plant", upgrade = TRUE)
# 
# test_low1 <- c(1,2)
# saveRDS(test_low1, "plant/test_low1.RDS")

devtools::load_all("plant")
# 
# test_low <- c(1,2)
# saveRDS(test_low, "plant/test_low.RDS")

devtools::load_all(".")


p0 <- scm_base_parameters("FF16w")
p0$disturbance_mean_interval <- 2
p0$max_patch_lifetime <- 30

p50_1 <- 2
p50_2 <- 4
p50_3 <- 6


p1 <- expand_parameters(trait_matrix(c(p50_1,p50_2, p50_3), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1,1,1))
env <- make_environment("FF16w", soil_initial_state = rep(0.5, 1), rainfall = 0.5)
# 
# 
# test_low4 <- c(1,2)
# saveRDS(test_low4, "plant/test_low4.RDS")


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

p1$strategies[[1]]$epsilon_leaf <- 0.0002
p1$strategies[[2]]$epsilon_leaf <- 0.0002
p1$strategies[[3]]$epsilon_leaf <- 0.0002

# 
# test_low5 <- c(1,2)
# saveRDS(test_low5, "plant/test_low5.RDS")


ctrl = scm_base_control()
# #
ctrl$ode_tol_abs <- 1e-4
ctrl$ode_step_size_initial <- 1e-04
ctrl$ode_step_size_min <- 1e-5
ctrl$schedule_eps <- 0.001
ctrl$ode_step_size_max <- 10
ctrl$equilibrium_eps <- 0.01
ctrl$ode_tol_rel <- 1e-3

system.time(result <- build_schedule(p1, env, ctrl))
system.time(result <- run_scm_collect(result, env))
saveRDS(result, "plant/result_low_2.5.RDS")



results_tidy <- 
  result %>% 
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

