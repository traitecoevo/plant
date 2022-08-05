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
  result_update %>% 
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


diff_value = 0.07

psi_stem_initial = 3 + diff_value

psi_stem_initial = psi_stem_next

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_0 <- l$calc_profit_Sperry_one_line(PPFD = 266.93, psi_soil = 2.31371, psi_stem = psi_stem_initial, k_l_max = 0.000324148)

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_1 <- l$calc_profit_Sperry_one_line(PPFD = 266.93, psi_soil = 2.31371, psi_stem = psi_stem_initial - diff_value, k_l_max = 0.000324148)

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_2 <- l$calc_profit_Sperry_one_line(PPFD = 266.93, psi_soil = 2.31371, psi_stem = psi_stem_initial + diff_value, k_l_max = 0.000324148)

first_dev = (y_2 - y_1)/(2*diff_value)
sec_dev = (y_2 - 2*y_0 + y_1)/(diff_value^2)

psi_stem_next = psi_stem_initial -  first_dev/sec_dev;

library(purrr)
library(tidyverse)

tibble(psi_stem = seq(2.31371, 5, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = l$calc_profit_Sperry_one_line(PPFD = 266.93, psi_soil = 2.31371, psi_stem = psi_stem, k_l_max= 0.000324148)) %>%
  ggplot(aes(x= psi_stem, y = profit)) +
  geom_line() +
  geom_point(aes(x= psi_stem_initial, y = y_0)) +
  geom_point(aes(x= psi_stem_initial - diff_value, y = y_1)) +
  geom_point(aes(x= psi_stem_initial + diff_value, y = y_2))
l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
l$psi_stem_next <- 3
l$optimise_psi_stem_Sperry_Newton_recall_one_line(266.93,2.31371, 0.000324148)




diff_value = 0.01

psi_stem_initial = 3.33733

psi_stem_initial = psi_stem_next

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_0 <- l$calc_profit_Sperry_one_line(PPFD = 266.838, psi_soil = 2.31554, psi_stem = psi_stem_initial, k_l_max = 0.000324148)

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_1 <- l$calc_profit_Sperry_one_line(PPFD = 266.838, psi_soil = 2.31554, psi_stem = psi_stem_initial - diff_value, k_l_max = 0.000324148)

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_2 <- l$calc_profit_Sperry_one_line(PPFD = 266.838, psi_soil = 2.31554, psi_stem = psi_stem_initial + diff_value, k_l_max = 0.000324148)

first_dev = (y_2 - y_1)/(2*diff_value)
sec_dev = (y_2 - 2*y_0 + y_1)/(diff_value^2)

psi_stem_next = psi_stem_initial -  first_dev/sec_dev;


diff_value = 0.01

psi_stem_initial = 3.33733

psi_stem_initial = psi_stem_next

l <- Leaf(vcmax = 100, p_50 = p50_2, c = 2.04, b = p1$strategies[[2]]$b, psi_crit = p1$strategies[[2]]$psi_crit , beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 0.6297342)
y_0 <- l$calc_profit_Sperry_one_line(PPFD = 266.838, psi_soil = 2.31554, psi_stem = psi_stem_initial, k_l_max = 0.000324148)

l <- Leaf(vcmax = 100, p_50 = p50_2, c = 2.04, b = p1$strategies[[2]]$b, psi_crit = p1$strategies[[2]]$psi_crit , beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 0.6297342)
y_1 <- l$calc_profit_Sperry_one_line(PPFD = 266.838, psi_soil = 2.31554, psi_stem = psi_stem_initial - diff_value, k_l_max = 0.000324148)

l <- Leaf(vcmax = 100, p_50 = p50_2, c = 2.04, b = p1$strategies[[2]]$b, psi_crit = p1$strategies[[2]]$psi_crit, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 0.6297342)
y_2 <- l$calc_profit_Sperry_one_line(PPFD = 266.838, psi_soil = 2.31554, psi_stem = psi_stem_initial + diff_value, k_l_max = 0.000324148)

first_dev = (y_2 - y_1)/(2*diff_value)
sec_dev = (y_2 - 2*y_0 + y_1)/(diff_value^2)

psi_stem_next = psi_stem_initial -  first_dev/sec_dev;



l <- Leaf(vcmax = 100, p_50 = p50_2, c = 2.04, b = p1$strategies[[3]]$b, psi_crit = p1$strategies[[3]]$psi_crit, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = p1$strategies[[3]]$K_s)


library(purrr)
library(tidyverse)

tibble(psi_stem = seq(2.31371, p1$strategies[[3]]$psi_crit, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = l$calc_profit_Sperry_one_line(PPFD = 266.93, psi_soil = 2.31371, psi_stem = psi_stem, k_l_max= 0.000324148)) %>%
  ggplot(aes(x= psi_stem, y = profit)) +
  geom_line() +
  geom_point(aes(x= psi_stem_initial, y = y_0)) +
  geom_point(aes(x= psi_stem_initial - diff_value, y = y_1)) +
  geom_point(aes(x= psi_stem_initial + diff_value, y = y_2))
l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
l$psi_stem_next <- 3
l$optimise_psi_stem_Sperry_Newton_recall_one_line(266.93,2.31371, 0.000324148)















psi_stem_initial = 3.5
diff_value = 0.01
psi_soil = 2.924

psi_stem_initial = psi_stem_next


l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_0 <- l$calc_profit_Sperry_one_line(PPFD = 500, psi_soil = psi_soil, psi_stem = psi_stem_initial, k_l_max = 0.000524148)

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_1 <- l$calc_profit_Sperry_one_line(PPFD = 500, psi_soil = psi_soil, psi_stem = psi_stem_initial - diff_value, k_l_max = 0.000524148)

l <- Leaf(vcmax = 100, p_50 = 1.731347, c = 2.04, b = 2.072101, psi_crit = 3.548059, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = 2)
y_2 <- l$calc_profit_Sperry_one_line(PPFD = 500, psi_soil = psi_soil, psi_stem = psi_stem_initial + diff_value, k_l_max = 0.000524148)


first_dev = (y_2 - y_1)/(2*diff_value)
sec_dev = (y_2 - 2*y_0 + y_1)/(diff_value^2)

psi_stem_next = psi_stem_initial -  first_dev/sec_dev;

library(purrr)
library(tidyverse)

tibble(psi_stem = seq(psi_soil, 4, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = l$calc_profit_Sperry_one_line(PPFD = 500, psi_soil = psi_soil, psi_stem = psi_stem, k_l_max= 0.000524148)) %>%
  ggplot(aes(x= psi_stem, y = profit)) +
  geom_line() +
  geom_point(aes(x= psi_stem_initial, y = y_0)) +
  geom_point(aes(x= psi_stem_initial - diff_value, y = y_1)) +
  geom_point(aes(x= psi_stem_initial + diff_value, y = y_2))



