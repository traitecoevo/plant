rm(list=ls())

source("tests/reference_leaf.R")

library(ggplot2)

l <- Leaf()

K_s = 2
h_v = 0.000157
h = 5

p_50 = 1.731347

k_l_max <- l$calc_k_l_max(K_s, h_v, h)

kg_to_g = 1000 #g kg^-1
gh2O_to_molh20 = 1/18.02 #mol g^-1

kg_2_mol_h20 <- kg_to_g * gh2O_to_molh20


kc_25 <- 404.9
ko_25 <- 278400
gamma_25 <- 42.75 
umol_per_mol_2_mol_per_mol <- 1e-6 # mol umol ^-1
atm_o2_kpa <- 21 #kPa
atm_kpa = 101.3 #kPa
kPa_2_Pa = 1000 #Pa/ kPa
umol_per_mol_2_Pa <- atm_kpa*kPa_2_Pa*umol_per_mol_2_mol_per_mol
km_25 <- (kc_25*umol_per_mol_2_Pa)*(1 + (atm_o2_kpa*kPa_2_Pa)/(ko_25*umol_per_mol_2_Pa))

inputs <- expand_grid(
  PPFD = c(300,600,900), 
  vcmax = c(30,100),
  vcmax_25_2_jmax_25 = 1.67,
  curv_fact = 0.90,
  a = 0.3,
  gamma_25 = 42.75,
  umol_per_mol_2_Pa = umol_per_mol_2_Pa,
  km_25 = km_25,
  psi_soil = seq(0, 3.5, 0.1), 
  k_l_max = c(0.5*6.28e-05,6.28e-05,2*6.28e-05),
  p_50 = 1.731347,
  c = 2,
  b = 2.072101,
  kg_2_mol_h20 = kg_2_mol_h20,
  umol_per_mol_2_mol_per_mol = 1e-6,
  atm_vpd = c(2), 
  ca =  c(30,40,50), 
  atm_kpa = 101.3, #kPa
  kPa_2_Pa = 1000, #Pa/ kPa
  psi_crit = b*(log(1/0.05))^(1/c)
) 

inputs %>%
  rowwise() %>%
  mutate(outputs = tibble(profit = l$optimise_profit_gss(PPFD = PPFD, vcmax = vcmax, vcmax_25_2_jmax_25 = vcmax_25_2_jmax_25, curv_fact = curv_fact, a = a, gamma_25 = gamma_25, umol_per_mol_2_Pa = umol_per_mol_2_Pa, km_25 = km_25, 
                                                         psi_soil = psi_soil, k_l_max = k_l_max, p_50 = p_50, c = c, b = b, kg_2_mol_h2o = kg_2_mol_h20, umol_per_mol_2_mol_per_mol = umol_per_mol_2_mol_per_mol, 
                                                         atm_vpd = atm_vpd, ca = ca, atm_kpa = atm_kpa, kPa_2_Pa = kPa_2_Pa, psi_crit = psi_crit), E = l$E, psi_stem = l$psi, g_c = l$g_c, ci = l$ci)) %>%
  unnest(outputs)-> outputs_cpp



outputs_cpp %>%
  left_join(outputs_R) %>%
  rowwise() %>%
  mutate(gc_replugged = l$calc_g_c(psi_soil = psi_soil, psi_stem = psi_stem_root, k_l_max = k_l_max, p_50 = p_50, c = c, b = b, atm_kpa = atm_kpa, kg_2_mol_h2o = kg_2_mol_h20, atm_vpd = atm_vpd)) -> outputs_cpp


plot(outputs_cpp$gc_replugged/outputs_cpp$g_c~outputs_cpp$psi_soil)

inputs %>%
  rowwise() %>%
  mutate(find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)) %>%
  mutate(calc_g_c(psi_stem= psi_stem_root, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, b = b, c = c))-> outputs_R

plot(outputs_cpp$gc_replugged/outputs_R$`calc_g_c(...)`~outputs_R$psi_soil)

plot(outputs_R/outputs_cpp$gc_replugged~outputs_R$psi_soil)



inputs %>%
  rowwise() %>%
  mutate(psi_stem = psi_soil+runif(1, 0, 2)) %>%
  mutate(tibble(gc_cpp = l$calc_g_c(psi_soil = psi_soil, psi_stem = psi_stem, k_l_max = k_l_max, p_50 = p_50, c = c, b = b, atm_kpa = atm_kpa, kg_2_mol_h2o = kg_2_mol_h20, atm_vpd = atm_vpd))) -> outputs_cpp

outputs_cpp %>%
  rowwise() %>%
  mutate(calc_g_c(psi_soil = psi_soil, psi_stem = psi_stem, k_l_max = k_l_max, b = b, c = c, atm_vpd = atm_vpd)) -> outputs_R


plot(outputs_R$gc_cpp/outputs_R$`calc_g_c(...)`~outputs_R$psi_soil)

inputs %>%
  rowwise() %>%
  mutate(outputs = tibble(profit = l$optimise_profit_gss(PPFD = PPFD, vcmax = vcmax, vcmax_25_2_jmax_25 = vcmax_25_2_jmax_25, curv_fact = curv_fact, a = a, gamma_25 = gamma_25, umol_per_mol_2_Pa = umol_per_mol_2_Pa, km_25 = km_25, 
                        psi_soil = psi_soil, k_l_max = k_l_max, p_50 = p_50, c = c, b = b, kg_2_mol_h2o = kg_2_mol_h20, umol_per_mol_2_mol_per_mol = umol_per_mol_2_mol_per_mol, 
                        atm_vpd = atm_vpd, ca = ca, atm_kpa = atm_kpa, kPa_2_Pa = kPa_2_Pa, psi_crit = psi_crit), E = l$E, psi_stem = l$psi, g_c = l$g_c, ci = l$ci)) %>%
  unnest(outputs) %>%
  rename(profit_cpp = profit,
         psi_stem_cpp = psi_stem)-> outputs_cpp

inputs %>%
  rowwise() %>%
  mutate(find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)) %>%
  rename(profit_R = profit_root,
         psi_stem_R = psi_stem_root,
         gc_R = gc_root) -> outputs_R

inputs %>%
  rowwise() %>%
  mutate(opt_psi_stem_gss(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)) -> outputs_R_manual

plot(outputs_cpp$E/outputs_R$E_root~outputs_cpp$psi_soil)
plot(outputs_cpp$g_c/outputs_R$gc_R~outputs_cpp$psi_soil)



plot(outputs_cpp$profit_cpp/outputs_R$profit_R~outputs_cpp$psi_soil)

plot(outputs_cpp$psi_stem_cpp/outputs_R$psi_stem_R~outputs_cpp$psi_soil)


plot(outputs_cpp$psi_stem_cpp,outputs_R$psi_stem_R)
plot(outputs_cpp$profit_cpp,outputs_R$profit_R)


left_join(outputs_cpp, outputs_R) %>%
  select(PPFD, vcmax, psi_soil, k_l_max, atm_vpd, profit_cpp, psi_stem_cpp, psi_stem_R, profit_R) %>%
  pivot_longer(cols = c(profit_cpp, psi_stem_cpp, psi_stem_R, profit_R)) %>%
  mutate(metric = ifelse(grepl("profit", name), "profit", "psi_stem")) %>%
  mutate(program = ifelse(grepl("cpp", name), "cpp", "R"))-> data


left_join(outputs_cpp, outputs_R) %>%
  select(PPFD, vcmax, psi_soil, k_l_max, atm_vpd, profit_cpp, psi_stem_cpp, psi_stem_R, profit_R, ca) -> data_wide

a_plot <- data_wide %>%
  ggplot() +
  geom_point(aes(x = psi_stem_cpp , y = psi_stem_R)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", col = "red") +
  theme_classic()

b_plot <- data_wide %>%
  ggplot() +
  geom_point(aes(x = profit_cpp , y = profit_R)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", col = "red") +
  theme_classic()

c_plot <- data_wide %>%
  ggplot() +
  geom_point(aes(x = psi_soil , y = psi_stem_R/psi_stem_cpp)) + 
  theme_classic()

d_plot <- data_wide %>%
  ggplot() +
  geom_point(aes(x = psi_soil , y = profit_R/profit_cpp)) + 
  theme_classic()

d_plot <- (data_wide$psi_stem_cpp / data_wide$psi_stem_R) %>% hist(main = "Psi_stem (cpp/R)")
e_plot <-(data_wide$profit_cpp / data_wide$profit_R) %>% hist(main = "Profit (cpp/R)")


png("R_vs_cpp_comparison.png")
cowplot::plot_grid(plotlist = list(a_plot, b_plot, c_plot, d_plot))
dev.off()


ggplot()+
  geom_point(aes(x = data_wide$profit_cpp, y = data_wide$profit_R)) + 
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", col = "red")
  
ggplot()+
  geom_point(aes(x = data_wide$psi_stem_cpp, y = data_wide$psi_stem_R)) + 
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", col = "red")

(data_wide$psi_stem_cpp / data_wide$psi_stem_R) %>% hist()
(data_wide$profit_cpp / data_wide$profit_R) %>% hist()

data_wide %>%
  mutate(profit_ratio = profit_cpp / profit_R) %>%
  ggplot(aes(x = psi_stem_cpp, y = psi_stem_R)) +
  geom_point() +
  facet_wrap(~profit_ratio < 0.8)

data_wide %>%
  mutate(profit_ratio = abs(log10(profit_cpp / profit_R))) %>%
  ggplot(aes(x = profit_cpp, y = profit_R)) +
  geom_point() +
  facet_wrap(~profit_ratio > 0.025)

                        
inputs <- expand_grid(
  PPFD = c(300,600,900), 
  vcmax = c(30,100),
  vcmax_25_2_jmax_25 = 1.67,
  curv_fact = 0.90,
  a = 0.3,
  gamma_25 = 42.75,
  umol_per_mol_2_Pa = umol_per_mol_2_Pa,
  km_25 = km_25,
  psi_soil = seq(0, 3.5, 0.1),
  k_l_max = c(0.5*6.28e-05,6.28e-05,2*6.28e-05),
  p_50 = 1.731347,
  c = 2,
  b = 2.072101,
  kg_2_mol_h20 = kg_2_mol_h20,
  umol_per_mol_2_mol_per_mol = 1e-6,
  atm_vpd = c(2), 
  ca =  c(30,40,50), 
  atm_kpa = 101.3, #kPa
  kPa_2_Pa = 1000, #Pa/ kPa
  psi_crit = b*(log(1/0.05))^(1/c)
)                 
  


inputs %>%
  rowwise() %>%
  mutate(psi_stem = psi_soil+runif(1, 0, 2)) %>%
  filter(psi_stem < psi_crit) %>%
  mutate(outputs = tibble(profit = l$calc_profit(PPFD = PPFD, vcmax = vcmax, vcmax_25_2_jmax_25 = vcmax_25_2_jmax_25, curv_fact = curv_fact, a = a, gamma_25 = gamma_25, umol_per_mol_2_Pa = umol_per_mol_2_Pa, km_25 = km_25, 
                                                         psi_soil = psi_soil, psi_stem = psi_stem, k_l_max = k_l_max, p_50 = p_50, c = c, b = b, kg_2_mol_h2o = kg_2_mol_h20, umol_per_mol_2_mol_per_mol = umol_per_mol_2_mol_per_mol, 
                                                         atm_vpd = atm_vpd, ca = ca, atm_kpa = atm_kpa, kPa_2_Pa = kPa_2_Pa, psi_crit = psi_crit))) -> outputs_cpp
  
outputs_cpp %>%
  rowwise() %>%
  mutate(calc_profit(psi_stem = psi_stem, psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)) -> outputs



plot(outputs$outputs$profit/outputs$`calc_profit(...)`~outputs$psi_soil)
