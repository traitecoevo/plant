devtools::load_all()
library(tidyverse)


# FF16
base_mulga_parameters <- function() {
  p0 <- scm_base_parameters("FF16", "FF16_Env")
  
  p0$strategy_default$lma <- 0.0645
  p0$strategy_default$hmat <- 5
  p0$strategy_default$rho <- 1100
  p0$strategy_default$narea <- 0.0032
  p0$strategy_default$a_l1 <- 2.17
  p0$strategy_default$a_l2 <- 0.5
  p0$strategy_default$omega <- 0.00000926
  p0
}



set_patch <- function(traits =  trait_matrix(c(0.0645), c("lma")),
                      B_lf1 = 1,
                      seed_rain = 100,
                      max_patch_lifetime = 250,
                      p0 = base_mulga_parameters(),
                      latitude = 28.182,
                      schedule = NULL,
                      optimise_schedule = FALSE,
                      schedule_reduction_factor = 5) {
  
  hyper_par_fn = make_FF16_hyperpar(B_lf1 = B_lf1, latitude = latitude)
  p1 <- expand_parameters(traits, p0, hyper_par_fn, mutant = FALSE)
  p1$seed_rain <- seed_rain
  p1$max_patch_lifetime <- max_patch_lifetime
  
  if(!is.null(schedule)) {
    p1$cohort_schedule_times = schedule
  }
  
  if(optimise_schedule) {
    # use built-in schedule
    nodes <- cohort_schedule_times_default(max_patch_lifetime)
    
    # then downsample, taking every nth integration node 
    nth_element <- function(vector, n = 1, starting_position = 1) { 
      vector[c(1:starting_position, seq(starting_position, length(vector), n))] 
    }
    
    p1$cohort_schedule_times[[1]] <- nth_element(nodes, schedule_reduction_factor)
    
    result <- build_schedule(p1)
  } else
    result <- p1

  return(result)
}

net_biomass <- function(schedule, B_lf1) {
  
  p <- set_patch(schedule = schedule, 
                 B_lf1 = B_lf1)
  
  # gather outputs at each time step
  results <- run_scm_collect(p) %>% 
    tidy_patch() %>%
    FF16_expand_state()

  v <- c("mass_leaf", "mass_bark", "mass_sapwood", "mass_heartwood")

  tissues <- results$species %>% 
    integrate_over_size_distribution() %>%
    slice(-n()) %>% # drop last node
    select(time,species, one_of(v)) %>%
    pivot_longer(cols=starts_with("mass"), names_to = "tissue") %>%
    mutate(across(tissue, factor, levels = v)) %>%

  return(tissues)
}


# Build schedule and visualise against tree yield formula
p <- set_patch(optimise_schedule = T)
tissues <- net_biomass(schedule = p$cohort_schedule_times,
                       B_lf1 = 0.6931932)

tyf <- function(t, m = 50, g=12.25, r=1){
  k = 2 * g - 1.25
  agb = r * m * exp(-k / t)
  # bgb = agb * .67
  
  return(agb / 10) # t.ha-1 -> kg.m-2
}

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  xlim(c(0,250)) +
  facet_wrap(~species)


# Optimising routines for B_lf1 search
source("scripts/optimise-bayesopt.R")

calibrate <- function(x, biomass_fn = tyf, 
                      schedule = p$cohort_schedule_times,
                      n_control_pts = 50,
                      match_tissue = "mass_heartwood",
                      show_error = FALSE,
                      ...) {
  
  tissues <- net_biomass(schedule, B_lf1 = x[1])
  
  # generate equidistant control points
  s = schedule[[1]]
  pts <- s[findInterval(seq(0, max(s), len = n_control_pts), s)]
  
  totals <- filter(tissues, 
                   tissue == match_tissue, 
                   time %in% pts) %>%
    mutate(total_biomass = biomass_fn(time)) %>%
    mutate(error = total_biomass - value)
  
  if(show_error)
    plot(totals$time, totals$error)
  
  # loss fn
  return(sum(abs(totals$error)))
}


# optimise using expected improvement acquisition fn
fit <- bayesian_optimize(calibrate, bounds = c(0.5, 1))

# re-run with current optima
tissues <- net_biomass(schedule = p$cohort_schedule_times,
                       B_lf1 =  fit$minima)

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  xlim(c(0,250)) +
  facet_wrap(~species)
 

# continue, re-using previous results
fit <- bayesian_optimize(calibrate, 
                         evals = fit$evaluations,
                         iter = 5,
                         bounds = c(0.5, 1))
