test <- c(1,2)
saveRDS(test, "test.RDS")

saveRDS(test, "plant/test.RDS")

install.packages("devtools")
library(devtools)

test2 <- c(1,2)

saveRDS(test2, "plant/test2.RDS")
packages <- c("Rcpp", "R6", "crayon", "nleqslv", "BB" ,"BH")
install.packages(packages)


test3 <- c(1,2)

saveRDS(test3, "plant/test3.RDS")

test4 <- rep(NA, length(packages))

for (i in 1:length(packages)){
test4[i] <- print(nzchar(system.file(package = packages[i])))
}

saveRDS(test4, "plant/test4.RDS")

devtools::install_deps("plant", upgrade = TRUE)

test5 <- c(1,2)

saveRDS(test5, "plant/test5.RDS")

devtools::load_all("plant")

test6 <- c(1,2)

saveRDS(test6, "plant/test6.RDS")

p0 <- scm_base_parameters("FF16w")
p0$disturbance_mean_interval <- 2
p0$max_patch_lifetime <- 0.1

p1 <- expand_parameters(trait_matrix(c(3,4), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1,1))
env <- make_environment("FF16w", soil_initial_state = rep(0.2, 1), rainfall = 1.5)

p1$strategies[[1]]$K_s <- 3
p1$strategies[[2]]$K_s <- 1



ctrl = scm_base_control()
# 
ctrl$ode_tol_abs <- ctrl$ode_tol_rel <- 1e-3
ctrl$ode_step_size_initial <- ctrl$ode_step_size_min <- 1e-5
ctrl$schedule_eps <- 0.01

# remember to `load_all` if you change the C++ code
# types <- extract_RcppR6_template_types(p1, "Parameters")
# scm <- do.call('SCM', types)(p1, env, ctrl)
# 
# scm$run()




test7 <- c(1,2)

saveRDS(test7, "plant/test7.RDS")

# test5 <- c(1,2)
# 
# saveRDS(test4, "plant/test5.RDS")

result <- run_scm_collect(p1, env, ctrl)
saveRDS(result, "plant/result.RDS")


library(tidyverse)
result_high <- readRDS("result_high.RDS")

results_tidy_high <- 
  result_high %>% 
  tidy_patch()


results_tidy_high%>% 
  FF16_expand_state() -> results_tidy_expanded_high



data_species_tot_high <- 
  results_tidy_expanded_high$species %>% 
  integrate_over_size_distribution()


  