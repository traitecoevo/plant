library(tidyverse)
devtools::load_all(".")



## growth 
sp = expand_parameters(trait_matrix(c(0.0825), "lma"), 
                      scm_base_parameters("FF16drivers"))

time_pts <- seq(0, 150)

sp$strategies[[1]]$growth_rate_x = time_pts
sp$strategies[[1]]$growth_rate_y = 1 + sin(time_pts)
sp$strategies[[1]]$is_variable_growth_rate = T

sp_smooth <- build_schedule(sp)

result <-   run_scm_collect(sp) %>%
  tidy_patch() %>% 
  purrr::pluck("species") %>% 
  tidyr::drop_na()

ggplot(result, aes(x = time, y = height, group = node, alpha = density)) +
  geom_line() +
  labs(
    x = "Time (yrs)",
    y = "Height (m)",
    alpha = "Density (individuals.m^-2)",
    subtitle = "Time-varying growth: 1 + sin(t)"
  ) +
  scale_alpha(na.value = 0) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 0.7)



## herbivory - leaf turnover 
sp = expand_parameters(trait_matrix(c(0.0825), "lma"), 
                       scm_base_parameters("FF16drivers"))

time_pts <- seq(0, 150)

sp$strategies[[1]]$herbivory_rate_x = time_pts
sp$strategies[[1]]$herbivory_rate_y = rep(2, 151)
sp$strategies[[1]]$is_variable_herbivory_rate = T
sp$strategies[[1]]$herbivory_size_threshold = 1.2

#sp_smooth <- build_schedule(sp)

result <-   run_scm_collect(sp) %>%
  tidy_patch() %>% 
  purrr::pluck("species") %>% 
  tidyr::drop_na()

ggplot(result, aes(x = time, y = height, group = node, alpha = density)) +
  geom_line() +
  labs(
    x = "Time (yrs)",
    y = "Height (m)",
    alpha = "Density (individuals.m^-2)",
    subtitle = "Double leaf turnover for cohorts < 1.2m"
  ) +
  scale_alpha(na.value = 0) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 0.7)


## mortality
# strategy defaults
mulga <- function() {
  p0 <- plant::scm_base_parameters("FF16drivers", "FF16_Env")
  
  p0$strategy_default$lma <- 0.0645
  p0$strategy_default$hmat <- 5
  p0$strategy_default$rho <- 1100
  p0$strategy_default$narea <- 0.0032
  p0$strategy_default$a_l1 <- 2.17
  p0$strategy_default$a_l2 <- 0.5
  p0$strategy_default$omega <- 0.00000926
  p0
}

sp = expand_parameters(trait_matrix(c(0.0645), "lma"), mulga())

time_pts <- seq(0, 150)

sp$strategies[[1]]$mortality_rate_x = time_pts
sp$strategies[[1]]$mortality_rate_y = c(rep(100, 60), rep(1, 91))
sp$strategies[[1]]$is_variable_mortality_rate = T
sp$strategies[[1]]$herbivory_size_threshold = 1.2

#sp_smooth <- build_schedule(sp)

result <-   run_scm_collect(sp) %>%
  tidy_patch() %>% 
  purrr::pluck("species") %>% 
  tidyr::drop_na() 

times_in <- function(df, times) {
  df %>%
    filter(floor(time) %in% times)
}

result %>%
  times_in(., c(10, 20, 40, 50, 80, 100)) %>%
  ggplot(., aes(x = height, y = log_density)) +
    geom_line() +
    facet_wrap(~ floor(time)) +
    theme(aspect.ratio = 1)


ggplot(result, aes(x = time, y = height, group = node, alpha = density, colour = density)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (yrs)",
    y = "Height (m)",
    colour = "Density (individuals.m^-2)",
    alpha = "Density (individuals.m^-2)",
    subtitle = "Time-varying mortality affecting < 1.2 m trees: "
  ) +
  scale_alpha(na.value = 0) +
  scale_colour_viridis_c() +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 0.7)

