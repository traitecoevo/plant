library(tidyverse)
devtools::load_all(".")

# base parameters
p0 <- scm_base_parameters("FF16drivers")

## Define sin function to vary through time from 0-amplitude 
f_sin <- function(time_yrs, period_yrs = 1.0, amplitude = 1.0) {
  amplitude * (1 + sin(2 * pi * time_yrs / period_yrs)) / 2.0
}

## Add species
sp <- expand_parameters(trait_matrix(c(0.0825), "lma"), p0)

# Varying growth rate
## Create time-varying driver 
time_pts <- seq(0, 150, length.out=1000)
## Multiplier on assimilation, varying between g_min and 1.0 of it's potential value
g_min <- 0.7

sp$strategies[[1]]$growth_rate_x <- time_pts
sp$strategies[[1]]$growth_rate_y = g_min + (1-g_min) * f_sin(time_pts, period_yrs = 10)
sp$strategies[[1]]$is_variable_growth_rate <- TRUE

## Check, plot through time
plot(sp$strategies[[1]]$growth_rate_x, sp$strategies[[1]]$growth_rate_y, type = "l")

## Force high resolution run, by adding more node introduction times to the default. This saves running Build_schedule

n_extra_nodes <- 500
sp$node_schedule_times[[1]] <- 
  c(p0$node_schedule_times_default, seq(0, p0$max_patch_lifetime, length.out = n_extra_nodes)) %>% sort() %>% unique()

## Run the stand
system.time({sp_smooth <- run_scm_collect(sp)})

## Process results
result <-
  sp_smooth %>%
  tidy_patch() %>% 
  purrr::pluck("species") %>% 
  tidyr::drop_na()

plot_size_distribution(result) +
  labs(title = "Growth 0.7 to 1")

ggsave(filename = "growth 0.75-1.png")

FF16_generate_stand_report(sp_smooth, overwrite = TRUE)


# Varying mortality - as above but with mortality rate
sp = expand_parameters(trait_matrix(c(0.0825), "lma"), 
                       scm_base_parameters("FF16drivers"))

## Multiplier on assimilation, varying between m_min and 1.0 of it's potential value
time_pts <- seq(0, 150, length.out = 1000)
m_min <- 0.7

sp$strategies[[1]]$mortality_rate_x <- time_pts
sp$strategies[[1]]$mortality_rate_y <- m_min + (1 - m_min) * f_sin(time_pts, period_yrs = 10)
sp$strategies[[1]]$is_variable_mortality_rate <- TRUE

system.time({sp_smooth <- run_scm_collect(sp)})

result <-  sp_smooth %>%
  tidy_patch() %>% 
  purrr::pluck("species") %>% 
  tidyr::drop_na()

plot_size_distribution(result) +
  labs(title = "Mortality 0.7 to 1")

