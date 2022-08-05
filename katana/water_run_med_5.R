# 
# test_start <- c(1,2)
# saveRDS(test_start, "plant/test_start.RDS")

install.packages("devtools")
library(devtools)
# 
# test_med3 <- c(1,2)
# saveRDS(test_med3, "plant/test_med3.RDS")


packages <- c("Rcpp", "R6", "crayon", "nleqslv", "BB" ,"BH")
install.packages(packages)
# 
# test_med2 <- c(1,2)
# saveRDS(test_med2, "plant/test_med2.RDS")

devtools::install_deps("plant", upgrade = TRUE)
# 
# test_med1 <- c(1,2)
# saveRDS(test_med1, "plant/test_med1.RDS")

devtools::load_all("plant")
# 
# test_med <- c(1,2)
# saveRDS(test_med, "plant/test_med.RDS")

p0 <- scm_base_parameters("FF16w")
p0$disturbance_mean_interval <- 2
p0$max_patch_lifetime <- 2.5

p50_1 <- 2
p50_2 <- 4

p1 <- expand_parameters(trait_matrix(c(p50_1, p50_2), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1,1))
env <- make_environment("FF16w", soil_initial_state = rep(0.2, 1), rainfall = 1)
# 
# 
# test_med4 <- c(1,2)
# saveRDS(test_med4, "plant/test_med4.RDS")
# 

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

p1$strategies[[1]]$b <- calc_vul_b(p50_1, p1$strategies[[1]]$c)
p1$strategies[[2]]$b <- calc_vul_b(p50_2, p1$strategies[[2]]$c)

p1$strategies[[1]]$psi_crit <- calc_psi_crit(p1$strategies[[1]]$b, p1$strategies[[1]]$c)
p1$strategies[[2]]$psi_crit <- calc_psi_crit(p1$strategies[[2]]$b, p1$strategies[[2]]$c)

# 
# test_med5 <- c(1,2)
# saveRDS(test_med5, "plant/test_med5.RDS")


ctrl = scm_base_control()
# 
ctrl$ode_tol_abs <- ctrl$ode_tol_rel <- 1e-3
ctrl$ode_step_size_initial <- ctrl$ode_step_size_min <- 1e-5
ctrl$schedule_eps <- 0.01

result <- run_scm_collect(p1, env, ctrl)
saveRDS(result, "plant/result_med_2.5.RDS")