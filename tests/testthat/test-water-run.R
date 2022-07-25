p0 <- scm_base_parameters("FF16w")
p0$disturbance_mean_interval <- 2
p0$max_patch_lifetime <- 0.1

p1 <- expand_parameters(trait_matrix(c(0.000157, 0.000157*2), "huber_value"), p0, mutant=FALSE, birth_rate_list = c(1,1))
env <- make_environment("FF16w", soil_initial_state = rep(0.2, 1), rainfall = 1.5)

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

result <- run_scm_collect(p1, env, ctrl)
saveRDS(result, "result2species.RDS")


saveRDS(result, "result2.RDS")
# result <-readRDS("result.RDS")
result <- readRDS("result2.RDS")
library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded



data_species_tot <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()

results_tidy_expanded$species %>%
  drop_na() %>%
  View()

results_tidy_expanded$species %>%
  drop_na() %>% View()
  group_by(step, time) %>%
  summarise(area_leaf = sum(area_leaf)) %>%
  pull(area_leaf) %>% diff() %>% min()
  ggplot(aes(x = time, y = area_leaf)) +
  geom_line()

results_tidy_expanded$species %>%
  drop_na() %>%
plot_size_distribution()

  results_tidy <- 
    result %>% 
    tidy_patch() %>%
    FF16_expand_state() 

  results_tidy$species %>%
    drop_na() %>% View()

  
  hist(results_tidy$env$canopy_openness)
  