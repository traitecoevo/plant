devtools::load_all(".")
devtools::install(".")
library(plant)
library(tidyverse)

p0 <- scm_base_parameters("FF16w")


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


# set default values 
run_FF16w_model <- function(inputs){
  
  
  ctrl = scm_base_control()
  
  for(v in names(ctrl)){
    if(length(inputs[[v]]) > 0){
      if(!is.na(inputs[[v]])){
        ctrl[[v]] <- inputs[[v]]
      }
    }
  }
  
  p0 <- scm_base_parameters("FF16w")
  p0$max_patch_lifetime <- inputs$max_patch_lifetime
  
  p50_1 <- inputs$p50_1
  p50_2 <- inputs$p50_2
  
  if(is.na(inputs$epsilon_leaf)){
    epsilon_leaf = 0.00001
  } else{
    epsilon_leaf = inputs$epsilon_leaf
  }
  
  ## increase birth rates to 10 or 50
  
  p1 <- expand_parameters(trait_matrix(c(p50_1, p50_2), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1,1))
  
  p1$strategies[[1]]$K_s <- p_50_2_K_s(p50_1)
  p1$strategies[[1]]$b <- calc_vul_b(p50_1, p1$strategies[[1]]$c)
  p1$strategies[[1]]$psi_crit <- calc_psi_crit(p1$strategies[[1]]$b, p1$strategies[[1]]$c)
  p1$strategies[[1]]$epsilon_leaf <- epsilon_leaf
  p1$strategies[[1]]$recruitment_decay <- inputs$recruitment_decay
  
  
  p1$strategies[[2]]$K_s <- p_50_2_K_s(p50_2)
  p1$strategies[[2]]$b <- calc_vul_b(p50_2, p1$strategies[[2]]$c)
  p1$strategies[[2]]$psi_crit <- calc_psi_crit(p1$strategies[[2]]$b, p1$strategies[[2]]$c)
  p1$strategies[[2]]$epsilon_leaf <- epsilon_leaf
  p1$strategies[[2]]$epsilon_leaf <- inputs$recruitment_decay
  
  env <- make_environment("FF16w", soil_initial_state = rep(inputs$initial_theta, 1), rainfall = inputs$rainfall)
  time1 <- system.time(result <- build_schedule(p1, env, ctrl))
  time2 <- system.time(result <- run_scm_collect(result, env, ctrl, collect_auxiliary_variables = TRUE))
  return(list(result, time1, time2))
}


poss_run_plant = safely(.f = run_FF16w_model, otherwise = "Error")


inputs <- expand_grid(rainfall = 2.93, p50_1 = 5.99, p50_2 = 8.46, initial_theta = c(0.1), max_patch_lifetime = 100, schedule_eps = NA, ode_tol_abs = NA, ode_tol_rel = NA, epsilon_leaf = 0.001, recruitment_decay = c(1)) %>%
  slice(sample(1:n()))

run_and_save_plant_water_model <- function(..., save = FALSE, dir = NULL){
  input_parameters <- tibble(...)

  if(save == TRUE){
  if(!file.exists(paste("examples_vignettes/outputs/", dir, sep = ""))){
  browser()

    dir.create(paste("examples_vignettes/outputs/", dir, sep = ""))
    
  }
  
  dir <- paste("examples_vignettes/outputs/", dir, "/", sep = "")
  
  filename <-paste("rainfall_", round(input_parameters$rainfall,2), 
                   "p501_", round(input_parameters$p50_1, 2), 
                   "p502_", round(input_parameters$p50_2, 2), 
                   "lifetime_", input_parameters$max_patch_lifetime,
                   "schedule_eps", round(input_parameters$schedule_eps, 6),
                   "ode_tol_abs", round(input_parameters$ode_tol_abs, 6),
                   "ode_tol_rel", round(input_parameters$ode_tol_rel, 6),
                   "epsilon_leaf", round(input_parameters$epsilon_leaf, 6),
                   "theta", input_parameters$initial_theta,
                   "recruitment_decay", input_parameters$recruitment_decay,
                   ".RDS", sep="")
                   

  path <- paste(dir, filename, sep="")
  
  if(!file.exists(path)){
  result <- poss_run_plant(input_parameters)
  result<-list(result, input_parameters)
  saveRDS(result, path)
  }
  else {
    print(paste("skipping", path, "already done", sep = " "))
  }
  } else {
  
  result <- poss_run_plant(input_parameters)
  result<-list(result, input_parameters)
  
  return(result)
  }
  
  
  return
}


system.time(result <- pmap(inputs, run_and_save_plant_water_model, save = FALSE, dir = "recruitment_dynamics"))

result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("species") %>%  
  ggplot(aes(x = time, y = opt_psi_stem_)) +
  geom_line(aes(group = interaction(species, node), colour = species)) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab("Years") +
  ylab(expression(paste(psi[opt], " (-MPa)", sep = ""))) -> a


result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("species") %>%
  dplyr::filter(!is.na(.data$density))  %>%
  dplyr::mutate(relative_log_density = rel(.data$log_density)) %>%
  ggplot(aes(x = time, y = opt_psi_stem_, alpha = relative_log_density)) +
  geom_line(aes(group = interaction(species, node), colour = species)) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab("Years") +
  ylab(expression(paste(psi[opt], " (-MPa)", sep = ""))) -> a

result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("species") %>%
  dplyr::filter(!is.na(.data$density))  %>%
  dplyr::mutate(relative_log_density = rel(.data$log_density)) %>%
  ggplot(aes(x = time, y = opt_psi_stem_, alpha = relative_log_density)) +
  geom_line(aes(group = interaction(species, node), colour = species)) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab("Years") +
  ylab(expression(paste(psi[opt], " (-MPa)", sep = ""))) + 
  ggplot2::labs(x = "Time (years)", color = "Species", alpha = "Relative log(Density [/m/m2])") + 
  ggplot2::scale_alpha(range = c(0.01, 1)) + ggplot2::theme_classic() -> a
  
png("examples_vignettes/outputs/psi_distribution.png", height = 800, width =  1200, res = 200)
a
dev.off()


result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("species") %>%
  dplyr::filter(!is.na(.data$density))  %>%
  dplyr::mutate(relative_log_density = rel(.data$log_density)) %>%
  ggplot(aes(x = time, y = height, alpha = relative_log_density)) +
  geom_line(aes(group = interaction(species, node), colour = species)) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab("Years") +
  ylab(expression(paste(Height, " (m)", sep = ""))) + 
  ggplot2::labs(x = "Time (years)", color = "Species", alpha = "Relative log(Density [/m/m2])") + 
  ggplot2::scale_alpha(range = c(0.01, 1)) + ggplot2::theme_classic() -> b

png("examples_vignettes/outputs/size_distribution.png", height = 800, width = 1200, res = 200)
b
dev.off()








result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("species") %>%
  plot_size_distribution() -> b

result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("env") %>%
  ggplot() +
  geom_line(aes(x=time, y = theta))+
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab("Years") +
  ylab(expression(paste(Soil~moisture~content, " (",m^{3}~H[2],"O",~m^{-3}~Soil, ")",sep = ""))) -> c

png("examples_vignettes/outputs/theta.png", height = 800, width = 800, res = 100)
c
dev.off()

result[[1]][[1]]$result[[1]] %>%
  tidy_patch() %>%
  FF16_expand_state() %>%
  pluck("species") %>%
  select(-c(opt_ci_,count)) %>%
  integrate_over_size_distribution() %>%
  ggplot(aes(x = time, y = area_leaf)) +
  geom_line(aes(group = species, colour = species)) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab("Years") +
  ylab(expression(paste(Total~leaf~area, " (", m^{-2}, ")",sep = ""))) -> b

png("opt_psi_stem.png", height = 500, width = 1000, res = 100)
a
dev.off()


png("leaf_area_theta.png", height = 500, width = 1000, res = 100)
cowplot::plot_grid(b,c, ncol=2, nrow =1)
dev.off()

