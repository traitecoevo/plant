devtools::load_all()
library(tidyverse)
library(tictoc)

source("scripts/optimise-bayesopt.R")


# strategy defaults
mulga <- function() {
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


create_environment <- function() FF16_make_environment()
create_site <- make_FF16_hyperpar

create_species <- function(lma = 0.0645, hmat = 5, eta = 12,
                           birth_rate = 100,
                           site) {
  
  traits = trait_matrix(c(lma, hmat, eta), c("lma", "hmat", "eta"))
  sp = expand_parameters(traits, mulga(), site, mutant = FALSE)
  
  sp$birth_rate <- birth_rate
  
  return(sp)
}

# This is clunky - assumes one schedule will be suitable even while optimising
# traits and site characteristics
create_schedule <- function(max_patch_lifetime = 250,
                            optimise_schedule = FALSE,
                            schedule_reduction_factor = 5) {
  
  site <- create_site()
  parameters <- create_species(site = site)
  
  # update patch longevity
  parameters$max_patch_lifetime <- max_patch_lifetime
  
  if(optimise_schedule) {

    # start with built-in schedule
    nodes <- node_schedule_times_default(max_patch_lifetime)
    
    # then downsample, taking every nth integration node 
    nth_element <- function(vector, n = 1, starting_position = 1) { 
      vector[c(1:starting_position, seq(starting_position, length(vector), n))] 
    }
    
    times <- nth_element(nodes, schedule_reduction_factor)
    
    parameters$node_schedule_times[[1]] <- times #1spp
    
    # and rebuild
    parameters <- build_schedule(parameters)
  }

  return(parameters)
}


calculate_net_biomass <- function(parameters, environment) {

  # gather outputs at each time step
  results <- run_scm_collect(parameters, environment) %>% 
    tidy_patch() %>%
    FF16_expand_state()
  
  v <- c("mass_leaf", "mass_bark", "mass_sapwood", "mass_heartwood")
  
  tissues <- results$species %>% 
    integrate_over_size_distribution() %>%
    select(time,species, one_of(v)) %>%
    pivot_longer(cols=starts_with("mass"), names_to = "tissue") %>%
    mutate(across(tissue, factor, levels = v)) %>%
    filter(time != max(time)) # drop last node
    
  return(tissues)
}





# Example -----------------------------------------------------------------


p <- create_schedule(max_patch_lifetime = 50, optimise_schedule = T)

tyf <- function(t, m = 50, g=12.25, r=1){
  k = 2 * g - 1.25
  agb = r * m * exp(-k / t)
  # bgb = agb * .67
  
  return(agb / 10) # t.ha-1 -> kg.m-2
}

# bounds = list(B_lf1 = list(min = 0.5, max = 1.5))

bounds = list(B_lf1 = list(min = log(0.5), max = log(1.5)),
              hmat = list(min = log(2), max = log(8)),
              eta = list(min = log(2), max = log(15)))


seed <- generate_seeds(calibrate, bounds, n = 6, target = tyf, parameters = p)

fit <- optim.EI(calibrate, tyf, bounds, n_iter = 4, parameters = p, evals = seed)

tic()
fit <- optim.EI(calibrate, tyf, bounds, n_iter = 50, parameters = p, 
                evals = fit$evaluations,
                gp = fit$gp)
toc()

plot_preds(fit$gp$pointer, bounds, fit$evaluations, par = 2,)

ggplot(data.frame(fit$evaluations),
       aes(x.B_lf1, x.eta, size = y, color = y)) +
  geom_point() +
  scale_color_viridis_c()


plot(diff(fit$evaluations$x[, 1]), fit$evaluations$y[-1])
plot(diff(fit$evaluations$x[, 2]), fit$evaluations$y[-1])

# deleteGPseps()

# # optimise using expected improvement acquisition fn
# fit <- bayesian_optimize(calibrate, bounds = c(0.5, 1))


# Test results
m = which.min(fit$evaluations$y)
pred <- exp(fit$evaluations$x)

best <- pred[m, ]


site <- create_site(B_lf1 = best[1], latitude = 28.182)
sp <- create_species(hmat = best[2],
                     eta = best[3],
                     site = site)

# sp <- create_species(site = site)

sp$node_schedule_times <- p$node_schedule_times

e <- create_environment()

# re-run with current optima
tissues <- calculate_net_biomass(sp, e) %>%
  filter(time != max(p$node_schedule_times[[1]]))

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  facet_wrap(~species)



# # continue, re-using previous results
# fit <- bayesian_optimize(calibrate, 
#                          evals = fit$evaluations,
#                          iter = 1,
#                          bounds = c(0.5, 1))
