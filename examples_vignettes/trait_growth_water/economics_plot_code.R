p0 <- scm_base_parameters("FF16w")
leaf <- Leaf(p0$strategy_default$vcmax, p0$strategy_default$p_50, p0$strategy_default$c, p0$strategy_default$b, p0$strategy_default$psi_crit, p0$strategy_default$K_s,
     p0$strategy_default$epsilon_leaf, p0$strategy_default$beta1, p0$strategy_default$beta2)

k_l_max = p0$strategy_default$K_s * p0$strategy_default$theta / 10
sapwood_volume_per_leaf_area = p0$strategy_default$theta * 10
leaf$set_physiology(1000, 0.0, k_l_max, 1, 40,sapwood_volume_per_leaf_area)
library(tidyverse)
tibble(psi_stem = seq(0, p0$strategy_default$psi_crit, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = leaf$profit_psi_stem_Bartlett(psi_stem),
         cost = leaf$hydraulic_cost_Bartlett(psi_stem),
         benefit = profit + cost) %>%
  pivot_longer(cols = c(profit, cost, benefit)) %>%
  ggplot(aes(x = psi_stem, y = value)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_line(aes(colour = name, group = name), size = 2) +
  theme_classic() +
  xlab(expression(paste(psi[soil]," (-MPa)"))) +
  ylab(expression(paste(Assimilated~CO[2], " (",mu,"mol",~m^{-2}~s^{-1},")"))) +
  theme(text = element_text(size = 18), legend.position = "none") 

tibble(psi_stem = seq(0, 3, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = leaf$profit_psi_stem_Bartlett(psi_stem),
         cost = leaf$hydraulic_cost_Bartlett(psi_stem),
         benefit = profit + cost) %>%
  pivot_longer(cols = c(profit, cost, benefit)) %>%
  filter(name %in% c("benefit")) %>%
  ggplot(aes(x = psi_stem, y = value)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_line(aes(colour = name, group = name), size = 2) +
  theme_classic() +
  xlab(expression(paste(psi[leaf]," (-MPa)"))) +
  ylab(expression(paste(Assimilated~CO[2], " (",mu,"mol",~m^{-2}~s^{-1},")"))) +
  theme(text = element_text(size = 25), legend.position = "none") +
  ylim(-10, 30) +
  scale_color_discrete() -> a


tibble(psi_stem = seq(0, 3, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = leaf$profit_psi_stem_Bartlett(psi_stem),
         cost = leaf$hydraulic_cost_Bartlett(psi_stem),
         benefit = profit + cost) %>%
  pivot_longer(cols = c(profit, cost, benefit)) %>%
  filter(name %in% c("benefit","cost")) %>%
  mutate(name = recode_factor(name, Benefit = "benefit", Cost = "cost")) %>%
  ggplot(aes(x = psi_stem, y = value)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_line(aes(colour = name, group = name), size = 2) +
  theme_classic() +
  xlab(expression(paste(psi[leaf]," (-MPa)"))) +
  ylab(expression(paste(Assimilated~CO[2], " (",mu,"mol",~m^{-2}~s^{-1},")"))) +
  theme(text = element_text(size = 25), legend.position = "none") +
  ylim(-10, 30) +
  scale_color_discrete() -> b


tibble(psi_stem = seq(0, 3, length.out = 100)) %>%
  rowwise() %>%
  mutate(profit = leaf$profit_psi_stem_Bartlett(psi_stem),
         cost = leaf$hydraulic_cost_Bartlett(psi_stem),
         benefit = profit + cost) %>%
  pivot_longer(cols = c(profit, cost, benefit)) %>%
  filter(name %in% c("benefit","profit","cost")) %>%
  mutate(name = recode_factor(name, Benefit = "benefit", Profit = "profit", Cost = "cost")) %>%
  ggplot(aes(x = psi_stem, y = value)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_line(aes(colour = name, group = name), size = 2) +
  theme_classic() +
  xlab(expression(paste(psi[leaf]," (-MPa)"))) +
  ylab(expression(paste(Assimilated~CO[2], " (",mu,"mol",~m^{-2}~s^{-1},")"))) +
  theme(text = element_text(size = 25), legend.position = "none") +
  ylim(-10, 30) -> c

png("examples_vignettes/outputs/ESA_2022/one.png", height = 1000, width = 1200, res = 180)
a
dev.off()
png("examples_vignettes/outputs/ESA_2022/two.png", height = 1000, width = 1200, res = 180)
b
dev.off()

png("examples_vignettes/outputs/ESA_2022/three.png", height = 1000, width = 1200, res = 180)
c
dev.off()


f <- function(x){
  profit = leaf$profit_psi_stem_Bartlett
  cost = leaf$hydraulic_cost_Bartlett
}
tibble(psi_stem = seq(0, p0$strategy_default$psi_crit, length.out = 100))  %>%
  mutate(benefit = profitz)


