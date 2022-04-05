# minimum reproducible example of initial size density problem

env <- make_environment("FF16")
ctrl <- scm_base_control()
load("scripts/initial_conditions.Rdata")

out <- run_scm_collect(scenario$optimised_parameters$parameters, env, ctrl)

results <- tidy_patch(out)

results %>%
  {.$species} %>%
  tidyr::drop_na() %>%
  plot_size_distribution() +
  coord_cartesian(expand = F) +
  scale_colour_manual(values = "black")


totals <- results %>%
  FF16_expand_state() %>%
  {.$species} %>%
  integrate_over_size_distribution() %>%
  slice(-n()) # final node looks weird

v <- c("mass_heartwood", "mass_sapwood", "mass_bark", "mass_leaf")

totals %>% 
  select(time, species, one_of(v)) %>%
  pivot_longer(cols=starts_with("mass"), names_to = "tissue") %>%
  mutate(across(tissue, factor, levels = v)) %>%
  ggplot(., aes(time, value, fill=tissue)) +
  geom_area() +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  coord_cartesian(expand = F)

totals %>% 
  select(time, area_leaf, individuals) %>%
  gather(metric, value, -time) %>%
  ggplot(., aes(time, value)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ metric, scales = "free")
