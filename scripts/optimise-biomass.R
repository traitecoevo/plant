devtools::load_all()
library(tidyverse)
library(tictoc)

source("scripts/optimise-bayesopt.R")

# Example -----------------------------------------------------------------


p <- create_schedule(max_patch_lifetime = 50, optimise_schedule = T)

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
              hmat = list(min = 2, max = 8),
              eta = list(min = 5, max = 12),
              k_2 = list(min = 0, max = 0.01))

# generate seed
seed <- generate_seeds(calibrate, bounds, n = 6, target = tyf, parameters = p)

# run smoke test - few iterations
fit <- bayesopt(calibrate, tyf, bounds, n_iter = 4, parameters = p, evals = seed)

# check for evidence of minima
fit

# IECI is fast but conservative
fit <- bayesopt(calibrate, tyf, bounds, n_iter = 10, parameters = p, 
                evals = fit$evaluations,
                gp = fit$gp, search_fn = "IECI")

# this is a bit janky
plot_preds(fit$gp$pointer, bounds, fit$evaluations, par = 1)


# update with EI which searches more broadly
tic()
fit <- bayesopt(calibrate, tyf, bounds, n_iter = 100, parameters = p, 
                evals = fit$evaluations,
                gp = fit$gp, search_fn = "EI")
toc()


ggplot(data.frame(fit$evaluations),
       aes(x_raw.B_lf1, x_raw.k_2, 
           size = -y, color = y)) +
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
                     site = site)


node_schedule_times_default <- p$node_schedule_times
sp$node_schedule_times <- p$node_schedule_times

# re-run with current optima
tissues <- calculate_net_biomass(sp) 

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  facet_wrap(~species)



