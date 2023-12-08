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


## mortality
sp = expand_parameters(trait_matrix(c(0.0825), "lma"), 
                       scm_base_parameters("FF16drivers"))

time_pts <- seq(0, 150)

sp$strategies[[1]]$mortality_rate_x = time_pts
sp$strategies[[1]]$mortality_rate_y = 1 + sin(time_pts)
sp$strategies[[1]]$is_variable_mortality_rate = T

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
    subtitle = "Time-varying mortality: 1 + sin(t)"
  ) +
  scale_alpha(na.value = 0) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 0.7)
