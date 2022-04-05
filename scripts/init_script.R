devtools::load_all()
library(tidyverse)
library(patchwork)

env <- make_environment("FF16")
ctrl <- scm_base_control()
p0 <- scm_base_parameters("FF16")

names <- c("height", "mortality", "fecundity", "area_heartwood",
           "mass_heartwood", "offspring_produced_survival_weighted",
           "log_density")

p1 <- expand_parameters(trait_matrix(c(0.0825), "lma"), p0, FF16_hyperpar,FALSE)
s = p1$strategies[[1]]


# log-normal
ln_init <- function(mean_height = 6, log_sd = 0.8,
                    density = 1,
                    min_h = 0.3920458,
                    max_h = 10,
                    n = 100) {
  
    log_mu <- log(mean_height)
  
    init_height <- seq(min_h, max_h, len = n)
    init_density <- dlnorm(init_height, log_mu, log_sd)
    
    # normalise to user density
    area_density <-   trapezium(init_height, init_density)
    init_density_scaled <- init_density / area_density * density
  
    # sense check initial leaf area
    la = (init_height / s$a_l1)^(1.0 / s$a_l2)
    leaf_area <- trapezium(init_height, init_density_scaled * la)
    
    # fit splines
    init_log_normal <- c(init_height, 
                         rep(0, n),
                         rep(0, n),
                         rep(0, n), 
                         rep(0, n),
                         rep(0, n),
                         log(init_density_scaled)) %>%
      matrix(., ncol = n, byrow = T)
    
    rownames(init_log_normal) <- names
    
    splines <- init_spline(list(init_log_normal), size_idx = 1)
    

    return(list(leaf_area = leaf_area,
                splines = splines))
}

# Lower mean density
init_df <- data.frame(mu = c(2, 2, 2, 2, 2, 2),
                      sd = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                      density = c(1.0, 0.5, 0.1, 0.05, 0.02, 0.01)) %>%
  mutate(init = pmap(list(mu, sd, density), ~ ln_init(..1, ..2, density = ..3))) %>%
  unnest_wider(init)


init_df
  
plant_log_console()

res <- map(init_df$splines, ~ build_schedule(p1, env, ctrl, splines = .x))

out <- map(res, ~ run_scm_collect(.$parameters, env, ctrl))

results <- map(out, tidy_patch)

fn <- function(init_df, tidy_patch, build_schedule_results) {
  scenario = sprintf("Mode: %d, density: %.3f, LA: %.2f", 
                 init_df$mu, init_df$density, init_df$leaf_area)
  
  n_nodes = length(build_schedule_results$parameters$cohort_schedule_times[[1]])
  diagnostics = sprintf("Complete: %s, # steps: %d, # nodes: %d",
                        build_schedule_results$complete,
                        build_schedule_results$n_steps,
                        n_nodes)
  
  tidy_patch %>%
    {.$species} %>%
    tidyr::drop_na() %>%
    plot_size_distribution() +
    coord_cartesian(expand = F) +
    labs(title = scenario,
         subtitle = diagnostics) +
    ggdark::dark_theme_bw()
}


y <- transpose(init_df)
#fn(y[[1]], results[[1]], res[[1]])

p <- pmap(list(y, results, res), ~ fn(..1, ..2, ..3))

wrap_plots(p) + 
  plot_layout(guides = 'collect')

# disable darkmode
ggdark::invert_geom_defaults() 



fn2 <- function(tidy_results) {
  tidy_results %>%
    FF16_expand_state() %>%
    {.$species} %>%
    integrate_over_size_distribution()
}

fn3 <- function(totals) {
  totals %>% 
    select(time, area_leaf, individuals) %>%
    gather(metric, value, -time) %>%
    ggplot(., aes(time, value)) +
      geom_line() +
      coord_cartesian(expand = F) +
      theme_classic() +
      facet_wrap(~ metric, scales = "free")
}

totals <- map(results, fn2)
p <- map(totals, fn3)

wrap_plots(p) + 
  plot_layout(guides = 'collect')

v <- c("mass_heartwood", "mass_sapwood", "mass_bark", "mass_leaf")

fn4 <- function(totals) {
  totals %>% 
    select(time, species, one_of(v)) %>%
    pivot_longer(cols=starts_with("mass"), names_to = "tissue") %>%
    mutate(across(tissue, factor, levels = v)) %>%
    ggplot(., aes(time, value, fill=tissue)) +
      geom_area() +
      labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
      theme_classic() + 
      coord_cartesian(expand = F)
}

p <- map(totals, fn4)

wrap_plots(p) + 
  plot_layout(guides = 'collect')


# Pathological scenarios
init_df <- data.frame(mu = c(10, 10, 10, 10, 10, 10),
                      sd = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                      density = c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001)) %>%
  mutate(init = pmap(list(mu, sd, density), ~ ln_init(..1, ..2, density = ..3))) %>%
  unnest_wider(init)


## Focus on initial cohorts
p1$strategies[[1]]$recruitment_decay <- 2
p1$strategies[[1]]$a_d0 = 0

## Adjust times
# set crude starting point
nth_element <- function(vec, interval){
  vec[seq(1, length(vec), interval)]
}
times <- p1$cohort_schedule_times_default
times_half <- nth_element(times, 2)
p1$cohort_schedule_times[[1]] <- times
