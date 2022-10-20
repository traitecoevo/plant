library(purrr)

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

calc_k_l_max <- function(K_s, h_v, h){
  K_s * h_v / h
}

c <- 2.04

h_v = 0.000157 #huber value (m^2 sapwood area m^-2 leaf area)
h = 5 #height/path length (m)
eta <- 12.0
eta_c <- 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta)

p50 <- 2

K_s = p_50_2_K_s(p50)
b = calc_vul_b(p50, c)
psi_crit = calc_psi_crit(b, c)
k_l_max = calc_k_l_max(K_s, h_v, h*eta_c)

l <- Leaf(vcmax = 30, p_50 = p50, c = 2.04, b = b, psi_crit = psi_crit, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = K_s)
l$set_physiology(PPFD = 900, psi_soil = 0, k_l_max = k_l_max)

library(tidyverse)
tibble(psi_soil = seq(0, 4.1, 0.01)) %>%
mutate(profit = map_dbl(.x = psi_soil, .f = ~l$calc_profit_Sperry_one_line(.x))) %>% View()
  ggplot(aes(x = psi_soil, y= profit)) +
  geom_line()
  

l <- Leaf(vcmax = 100, p_50 = 2, c = 2.04, b = p1$strategies[[1]]$b, psi_crit = p1$strategies[[1]]$psi_crit, beta=15000, beta_2 = 1, huber_value = 0.000157, K_s = p_50_2_K_s(2))
l$set_physiology(PPFD = 1.21886, psi_soil = 1.63328, k_l_max = 0.000843655-05)
l$PPFD_

psi_soil = 0
diff_value = 0.001
psi_stem_initial = 1.64477
x_1 = 1.64377
x_2 = 1.64577

y_0 = l$calc_profit_Sperry_one_line(psi_stem_initial)
y_1 = l$calc_profit_Sperry_one_line(x_1)
y_2 = l$calc_profit_Sperry_one_line(x_2)

psi_stem_initial - ((y_2 - y_1)/(2*diff_value))/((y_2 - 2*y_0 + y_1)/(diff_value^2))

l$optimise_psi_stem_Sperry_Newton_recall_one_line(NA)


