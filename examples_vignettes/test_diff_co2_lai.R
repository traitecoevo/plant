p0 <- scm_base_parameters("FF16w")
p0$max_patch_lifetime <- 50
ctrl = scm_base_control()


p1 <- expand_parameters(trait_matrix(c(0.0825), "lma"), p0, mutant=FALSE, birth_rate_list = c(1))
env <- make_environment("FF16w", soil_initial_state = rep(0.5, 1), rainfall = 2, co2 = 50)
time1 <- system.time(result <- build_schedule(p1, env, ctrl))
time2 <- system.time(result_high <- run_scm_collect(result, env, ctrl, collect_auxiliary_variables = TRUE))

devtools::install()
3
library(plant)
library(tidyverse)
result_high %>%
  tidy_patch %>% 
  FF16_expand_state() %>%
  pluck("species") %>% 
  select(-opt_ci_, -count) %>%
  integrate_over_size_distribution() %>%
  mutate(ca = 50) -> high_ca


result_low %>%
  tidy_patch %>% 
  FF16_expand_state() %>%
  pluck("species") %>% 
  select(-opt_ci_, -count) %>%
  integrate_over_size_distribution() %>%
  mutate(ca = 35) -> low_ca

low_ca %>% left_join(result_low$env %>%
                       tidy_env() %>% pluck("soil_moist")) %>%
  ggplot(aes(x=time,y=soil_moist)) +
  geom_line() -> low_ca_soil_moist
  
high_ca %>%
  bind_rows(low_ca)%>% left_join(result_high$env %>%
                       tidy_env() %>% pluck("soil_moist") %>% mutate(ca = 50) %>%
                        bind_rows(result_low$env %>%
                                    tidy_env() %>% pluck("soil_moist") %>% mutate(ca = 35))) %>%
  mutate(ca = factor(ca)) %>%
  ggplot(aes(x=time,y=soil_moist)) +
  geom_line(aes(colour = ca, group = ca)) +
  theme_classic()+
  theme(text = element_text(size=16)) +
  labs(colour = expression(paste(CO[2]," (", Pa, ")"))) +
  xlab("Time (yrs)") +
  ylab(expression(paste(theta, " (",m^3~H[2],"O ",m^{-3}~soil,")"))) -> soil_moist_plot


high_ca %>%
  bind_rows(low_ca) %>%
  mutate(ca  = factor(ca)) %>%
  ggplot(aes(x = time, y = area_leaf)) +
  geom_line(aes(col = ca, group = ca)) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(colour = expression(paste(CO[2]," (", Pa, ")"))) +
  xlab("Time (yrs)") +
  ylab(expression(paste("Leaf area index"))) -> co2_plot

cowplot::plot_grid(co2_plot, low_ca_soil_moist, ncol=1)


  
  png("examples_vignettes/outputs/co2_lai_plot_w_soil_moist.png", height = 2400, width = 2000, res = 300)
  cowplot::plot_grid(co2_plot, soil_moist_plot, ncol=1)
  dev.off()
                  