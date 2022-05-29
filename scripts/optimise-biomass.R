devtools::load_all()
library(tidyverse)
library(tictoc)

source("scripts/optimise-bayesopt.R")

# Example -----------------------------------------------------------------


p <- create_schedule(max_patch_lifetime = 200, optimise_schedule = T)

# fixed parameters, only vary belowground effects
site <- create_site(latitude = 28.182)
sp <- create_species(#k_2 = 0.05,
                     site = site)

run_scm_collect(sp) %>%
  tidy_patch() %>%
  pluck("species") %>%
  drop_na() %>%
  plot_size_distribution()


# calibrate against tree yield formula for a given max. bio
tyf <- function(t, m = 50, g=12.25, r=1){
  k = 2 * g - 1.25
  agb = r * m * exp(-k / t)
  # bgb = agb * .67
  
  return(agb / 10) # t.ha-1 -> kg.m-2
}


# parameter bounds on logscale
bounds = list(B_lf1 = list(min = 0.5, max = 1),
              hmat = list(min = 4, max = 8),
              eta = list(min = 5, max = 12),
              k_2 = list(min = 0, max = 0.5),
              birth_rate = list(min = 0.01, max = 2))

# generate seed
seed <- generate_seeds(calibrate, bounds, n = 30, target = tyf, parameters = p)

# run smoke test - few iterations
fit <- bayesopt(calibrate, tyf, bounds, n_iter = 4, parameters = p, evals = seed)

# check for evidence of minima
fit

# update with EI which searches broadly, but doesn't like boundaries
tic()
fit <- bayesopt(calibrate, tyf, bounds, n_iter = 20, parameters = p, 
                evals = fit$evaluations,
                gp = fit$gp, search_fn = "EI")
toc()


# IECI finds optima quickly, but explores conservatively
fit <- bayesopt(calibrate, tyf, bounds, n_iter = 5, parameters = p, 
                evals = fit$evaluations,
                gp = fit$gp, search_fn = "IECI")


ggplot(data.frame(fit$evaluations),
       aes(x_raw.B_lf1, x_raw.birth_rate, 
           size = -y, color = -y)) +
  geom_point() +
  scale_color_viridis_c()

# clean up
# deleteGPseps()
# deleteGPs()


# Test results
m = which.min(fit$evaluations$y)
pred <- unscale(fit$evaluations$x_raw, bounds)

best <- pred[m, ]

site <- create_site(B_lf1 = best[1],
                    latitude = 28.182)

sp <- create_species(hmat = best[2], 
                     eta = best[3],
                     k_2 = best[4],
                     birth_rate = best[5],
                     site = site)

sp$max_patch_lifetime <- 50
sp$node_schedule_times <- p$node_schedule_times

# re-run with current optima
tissues <- calculate_net_biomass(sp) 

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  facet_wrap(~species)

res <- run_scm_collect(sp) %>%
  tidy_patch() %>%
  pluck("species") %>%
  drop_na() 

plot_size_distribution(res)

res %>%
  integrate_over_size_distribution() %>%
  select(time, individuals) %>%
  ggplot(., aes(time, individuals)) +
  geom_line() +
  theme_classic()


# increase birthrate ------------------------------------------------------
site <- create_site(B_lf1 = best[1],
                    latitude = 28.182)

locked_sp <- create_species(hmat = best[2], 
                     eta = best[3],
                     k_2 = best[4],
                     birth_rate = 10 * best[5],
                     site = site)

locked_sp$max_patch_lifetime <- 50
locked_sp$node_schedule_times <- p$node_schedule_times

tissues <- calculate_net_biomass(locked_sp) 

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  facet_wrap(~species)


res <- run_scm_collect(locked_sp) %>%
  tidy_patch() %>%
  pluck("species") %>%
  drop_na() 


plot_size_distribution(res)

res %>%
    integrate_over_size_distribution() %>%
    select(time, individuals) %>%
    ggplot(., aes(time, individuals)) +
      geom_line() +
      theme_classic()


# add thinning treatment --------------------------------------------------
unlocked_sp <- locked_sp
unlocked_sp$strategies[[1]]$is_variable_mortality_rate = T
unlocked_sp$strategies[[1]]$mortality_rate_x = seq(0, 200, len = 200)
unlocked_sp$strategies[[1]]$mortality_rate_y = 0.75 * exp(-0.1 * unlocked_sp$strategies[[1]]$mortality_rate_x)

tissues <- calculate_net_biomass(unlocked_sp) 

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  facet_wrap(~species)


res <- run_scm_collect(unlocked_sp) %>%
  tidy_patch() %>%
  pluck("species") %>%
  drop_na() 

plot_size_distribution(res)

res %>%
  integrate_over_size_distribution() %>%
  select(time, individuals) %>%
  ggplot(., aes(time, individuals)) +
  geom_line() +
  theme_classic()
