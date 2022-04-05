devtools::load_all()
library(dfoptim)
library(tidyverse)

# Tree yield formula
tyf <- function(t, m = 500, g=12.25, r=1){
  k = 2 * g - 1.25
  agb = r * m * exp(-k / t)
  # bgb = agb * .67
  
  return(agb / 10) # t.ha-1 -> kg.m-2
}

plot(1:1000, tyf(1:1000), xlab = "Years", ylab = "Biomass")


# FF16
base_mulga_parameters <- function() {
  p0 <- scm_base_parameters("FF16", "FF16_Env")
  
  p0$strategy_default$lma <- 0.0645
  p0$strategy_default$hmat <- 5
  p0$strategy_default$rho <- 1100
  p0$strategy_default$narea <- 0.0032
  p0$strategy_default$a_l1 <- 2.17
  p0$strategy_default$a_l2 <- 0.5
  p0$strategy_default$omega <- 0.00000926
  p0
}

nth_element <- function(vector, n = 1, starting_position = 1) { 
  vector[c(1:starting_position, seq(starting_position, length(vector), n))] 
}

set_patch <- function(
    traits =  trait_matrix(c(0.0645), c("lma")),
    B_lf1 = 1,
    seed_rain = 100,
    max_patch_lifetime = 250,
    p0 = base_mulga_parameters(),
    find_equilbrium = FALSE,
    optimise_schedule = FALSE,
    latitude = 28.182,
    nth_reduction = 5,
    schedule = NULL
) {
  
  hyper_par_fn = make_FF16_hyperpar(B_lf1 = B_lf1, latitude = latitude)
  p1 <- expand_parameters(traits, p0, hyper_par_fn, mutant = FALSE)
  p1$seed_rain <- seed_rain
  p1$max_patch_lifetime <- max_patch_lifetime
  
  if(!is.null(schedule)) {
    p1$cohort_schedule_times = schedule
  }
  
  if(find_equilbrium)
    result <- equilibrium_seed_rain(p1)
  else if(optimise_schedule) {
    integration_points <- cohort_schedule_times_default(max_patch_lifetime)
    p1$cohort_schedule_times[[1]] <- nth_element(integration_points, nth_reduction)
    
    result <- build_schedule(p1)
  } else
    result <- p1

  return(result)
}

p <- set_patch(optimise_schedule = T)
  
net_biomass <- function(schedule, lma) {
  
  p <- set_patch(schedule = schedule, 
                 traits = trait_matrix(c(lma), c("lma")))
  
  # gather outputs at each time step
  results <- run_scm_collect(p) %>% 
    tidy_patch() %>%
    FF16_expand_state()

  v <- c("mass_leaf", "mass_bark", "mass_sapwood", "mass_heartwood")

  tissues <- results$species %>% 
    integrate_over_size_distribution() %>%
    slice(-n()) %>% # drop last node
    select(time,species, one_of(v)) %>%
    pivot_longer(cols=starts_with("mass"), names_to = "tissue") %>%
    mutate(across(tissue, factor, levels = v)) %>%

  return(tissues)
}

tissues <- net_biomass(schedule = p$cohort_schedule_times,
                       lma = 0.0645)

ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  xlim(c(0,250)) +
  facet_wrap(~species)


calibrate <- function(x, biomass_fn = tyf, 
                      schedule = p$cohort_schedule_times,
                      match_tissue = "mass_heartwood",
                      ...) {
  
  tissues <- net_biomass(schedule, lma= x[1])
  
  totals <- filter(tissues, tissue == match_tissue) %>%
    mutate(total_biomass = biomass_fn(time, m = x[2])) %>%
    mutate(error = abs(total_biomass - value) * seq(0, 1, len = n()))
  
  return(-log(sum(totals$error)))
}

init = c(0.065, 250)
calibrate(init)
  
optimised <- nmkb(init, calibrate, 
                  control = list(maximize = TRUE),
                  lower=c(0.06, 0), upper = c(0.07, 1e4),
                  biomass_fn = tyf)

tissues <- net_biomass(schedule = p$cohort_schedule_times,
                       lma = optimised$par[1])


ggplot(tissues, aes(time, value)) +
  geom_area(aes(fill = tissue)) +
  geom_function(fun=tyf, linetype = "dashed") +
  geom_function(fun= ~ tyf(., m = 450), 
                color = "red", linetype = "dashed") +
  labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
  theme_classic() + 
  xlim(c(0,100)) +
  facet_wrap(~species)
 
total %>% 
  ggplot(aes(time, individuals, colour=species)) +
  geom_line()

np <- 10
set.seed(123)

p0 <- rnorm(np)
xm1 <- nmk(fn=rosbkext, par=p0) # maximum `fevals' is not sufficient to find correct minimum
xm1b <- nmkb(fn=rosbkext, par=p0, lower=-2, upper=2)
