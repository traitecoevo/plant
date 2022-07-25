result_high <- readRDS("result_high.RDS")


result_high$







result_mid <- readRDS("result_med.RDS")
result_low <- readRDS("result_low.RDS")

result_high$env


library(tidyverse)

results_tidy <- 
  result_high %>% 
  tidy_patch()

results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded

data_species_tot_high <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution() %>%
  mutate(rain = "high")


library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded



data_species_tot_high <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_high %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()


result <- readRDS("result_med_2.5.RDS")

library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded



data_species_tot_med <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_med %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()


result <- readRDS("result_low_2.5.RDS")

library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded



data_species_tot_low <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_low %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()

data_species_tot_high %>%
  mutate(prec = "3 m") -> data_species_tot_high


data_species_tot_med %>%
  mutate(prec = "1 m ") -> data_species_tot_med


data_species_tot_low %>%
  mutate(prec = "0.5 m") -> data_species_tot_low

data_species_tot_high %>%
  bind_rows(data_species_tot_med) %>%
  bind_rows(data_species_tot_low) %>% 
  ggplot(aes(time, area_leaf, colour = prec, linetype = species))+
  geom_line()



data_species_tot_med %>%
  bind_rows(data_species_tot_low) %>%
  mutate(p_50 = if_else(species == 1, "2", "4")) %>%
  group_by(time, prec) %>%
  mutate(total_area_leaf = sum(area_leaf)) %>%
  ungroup() %>%
  mutate(frac_area_leaf = area_leaf/total_area_leaf) %>%
  ggplot(aes(time, frac_area_leaf, colour = prec, linetype = p_50))+
  geom_line()


data_species_tot_high %>%
  bind_rows(data_species_tot_med) %>%
  bind_rows(data_species_tot_low) %>%
  group_by(time, prec) %>%
  mutate(total_area_leaf = sum(area_leaf)) %>%
  ungroup() %>%
  ggplot(aes(time, total_area_leaf, colour = prec))+
  geom_line()



result <- readRDS("result_med.RDS")

library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded

data_species_tot_med <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_med %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()


result <- readRDS("result_low.RDS")

library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded

data_species_tot_med <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_med %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()


result <- readRDS("result_low_2.5.RDS")

library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded

data_species_tot_med <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_med %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()

result <- readRDS("result_med_2.5.RDS")

library(tidyverse)

results_tidy <- 
  result %>% 
  tidy_patch()


results_tidy%>% 
  FF16_expand_state() -> results_tidy_expanded

data_species_tot_med <- 
  results_tidy_expanded$species %>% 
  integrate_over_size_distribution()

library(ggplot2)

data_species_tot_med %>% 
  ggplot(aes(time, area_leaf, colour = species)) +
  geom_line()



