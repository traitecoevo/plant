
f <- function(x){
  leaf <- Leaf(vcmax = 300, p_50 = p0$strategy_default$p_50, c = p0$strategy_default$c, b = p0$strategy_default$b, psi_crit = p0$strategy_default$psi_crit, K_s = p0$strategy_default$K_s, epsilon_leaf = p0$strategy_default$epsilon_leaf, beta1 = p0$strategy_default$beta1,beta2=  p0$strategy_default$beta2)
  leaf$set_physiology(PPFD = 500, psi_soil = 1, leaf_specific_conductance_max = k_l_max, atm_vpd = 2, ca = 40, sapwood_volume_per_leaf_area = sapwood_volume)
  cost = leaf$hydraulic_cost_Bartlett(x)
  leaf$set_leaf_states_rates_from_psi_stem(x)
  benefit = leaf$assim_colimited_
  leaf$optimise_psi_stem_Bartlett()
  opt_psi_stem_ = leaf$opt_psi_stem_
  opt_profit = leaf$profit_
  tibble(benefit, cost, opt_psi_stem_,opt_profit)
}


p0 = scm_base_parameters("FF16w")

k_l_max = p0$strategy_default$K_s*p0$strategy_default$theta/1.75923
sapwood_volume = p0$strategy_default$theta*1.75923

tibble(psi_stem = seq(1, p0$strategy_default$psi_crit, length.out = 100)) %>%
  rowwise() %>%
  mutate(economics = f(psi_stem)) %>%
  unnest(economics) %>%
  mutate(profit = benefit - cost) %>%
  pivot_longer(cols = c(benefit, cost, profit)) %>%
  ggplot(aes(x = psi_stem, y = value)) +
  geom_line(aes(group = name, colour =name)) +
  geom_point(aes(x = opt_psi_stem_, y = opt_profit))


rm(list=ls())



out <- optim(rowMeans(bounds), f, method="L-BFGS-B",
             lower=bounds[, "lower"], upper=bounds[, "upper"],
             control=list(fnscale=-1, factr=1e10))
p0 = scm_base_parameters("FF16w")

f <- function(x){
  leaf <- Leaf(1000, 1.0, 2.04, 2.0, 3.42, 2, 0.001, 20000, 1.5)
  leaf$set_physiology(1000, 1, k_l_max, 2, 40, sapwood_volume)
  cost = leaf$hydraulic_cost_Bartlett(x)
  leaf$set_leaf_states_rates_from_psi_stem(x)
  benefit = leaf$assim_colimited_
  leaf$optimise_psi_stem_Bartlett()
  opt_psi_stem_ = leaf$opt_psi_stem_
  opt_profit = leaf$profit_
  tibble(benefit, cost, opt_psi_stem_,opt_profit)
}


trait_by_height <- function(...){

inputs <- tibble(...)

env <- FF16w_make_environment(vpd = inputs$vpd)
env$set_soil_water_state(env$soil_moist_from_psi(inputs$psi*1e06)) 
env$set_fixed_environment(inputs$E, 150)
  
p0 <- scm_base_parameters("FF16w")

f <- function(x){
s <- strategy(trait_matrix(x,  inputs$trait_name), p0, birth_rate_list = 1)
s$hmat <- 50
indv <- FF16w_Individual(s)
res <- plant::grow_individual_to_height(indv, heights = inputs$height, env, time_max =100)
res <- res$rate[colnames(res$rate) == "height"]
if(is.na(res)){
res <-  0
}
return(res)
}

ret <- optim(par= 0.000157, f, method="L-BFGS-B",
      control=list(fnscale=-1, trace = 1), lower = 1e-07, upper = 10)

tibble(trait = ret$par, growth = ret$value)
}

future::plan(multisession, workers = 6)

expand_grid(trait_name = "theta", E = c(1), psi = c(0.1,0.5,1,2), height = c(0.5,1,2), vpd = c(1,2,3,4,5)) %>%
  mutate(results = pmap(., trait_by_height)) -> outputs

expand_grid(trait_name = "theta", E = c(1), psi = c(0.1, 0.5, 1, 2), height = c(0.5, 1, 2, 5), vpd = c(1,2,3,5)) %>%
  mutate(results = pmap(., trait_by_height)) -> outputs

outputs %>%
  unnest(results)

outputs %>%
  unnest(results) %>%
  mutate(trait = ifelse(is.na(growth), NA, trait)) %>% 
  ggplot(aes(x = vpd, y = height, fill = trait)) +
  geom_tile() +
  facet_wrap(~psi) +
  scale_fill_continuous(trans = "log")



env <- FF16w_make_environment(vpd = c(1))
env$set_soil_water_state(env$soil_moist_from_psi(0.5*1e06)) 
env$set_fixed_environment(1, 150)

f <- function(x){
  s <- strategy(trait_matrix(x,  "theta"), p0, birth_rate_list = 1)
  indv <- FF16w_Individual(s)
  res <- plant::grow_individual_to_height(indv, heights = 0.5, env, time_max =100)
  res$rate[colnames(res$rate) == "height"]
}

tibble(K_s = c(0.000157,seq(0.000001, 0.01, length.out =50))) %>%
  rowwise() %>%
  mutate(height_growth = f(K_s)) %>%
  ggplot(aes(x = K_s, y = height_growth)) +
  geom_line() +
  geom_point(aes(x= 0.000157, y= 1.86))
  
