devtools::load_all()

# TODO: fix for mutants
collect_state <- function(scm) {
  res <- list(scm$state)
  while (!scm$complete) {
    scm$run_next()
    res <- c(res, list(scm$state))
  }
  return(res)
}

run_one_patch <- function(birth_rate = 1, lifetime = 105.32) {

  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime

  p1 <- expand_parameters(trait_matrix(0.8, "lma"), p, mutant = F,
                          birth_rate_list = birth_rate)

  e <- make_environment("FF16")

  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = T

  types <- extract_RcppR6_template_types(p1, "Parameters")
  scm <- do.call('SCM', types)(p1, e, ctrl)
  
  resident <- collect_state(scm)

  time <- purrr::map_dbl(resident, `[[`, "time") %>%
    data.frame(time = .) %>%
    dplyr::mutate(step = seq_along(time))
  
  state <- purrr::map(resident, `[[`, "species") %>%
    purrr::map_df(.,~ purrr::pluck(., 1) %>% t() %>%
                    tibble::as_tibble(), .id="step") %>%
    dplyr::mutate(dplyr::across(step, readr::parse_number)) %>%
    dplyr::group_by(step) %>%
    dplyr::mutate(node = seq_len(dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(time, by = "step")
    
  resident_rr <- scm$net_reproduction_ratios

  # use stateful SCM for efficient mutant runs
  scm$run_mutant(p1$strategies, append = F, update_schedule = F)
  
  mutant_rr <- scm$net_reproduction_ratios
  
  mutant_endpoint <- scm$state$species %>% 
    purrr::pluck(., 1) %>%
    t() %>% tibble::as_tibble() %>%
    dplyr::mutate(time = scm$time, 
                  node = seq_along(time))
#  sink()
  return(list(fitnesses = c(resident_rr, mutant_rr),
              resident_state = state,
              mutant_endpoint = mutant_endpoint))
}

result <- run_one_patch(1, 1)

# [1] 7.066999e-24 7.066999e-24
result$fitnesses

result <- run_one_patch(1)

# [1] 47.66136 47.66136
result$fitnesses

# test over grid
birth_rates <- c(100, 50, 10, 5, 1, 0.1, 0.01, 0.001)

grid_search <- purrr::map(birth_rates, run_one_patch)


fitnesses <- purrr::map_df(grid_search,
                           ~ purrr::pluck(., "fitnesses") %>% 
                             tibble::tibble(., rownames = c("resident_rr", "mutant_rr"),
                                            .name_repair = ~ c("net_reproduction_ratio", "invasion_type")),
                           .id = "grid_step") %>%
  tidyr::spread(invasion_type, net_reproduction_ratio) %>%
  dplyr::mutate(birth_rate = birth_rates,
                unit = "seeds_per_m2_yr") %>%
  dplyr::group_by(grid_step, birth_rate, unit) %>%
  dplyr::mutate(prop = scales::percent(mutant_rr / resident_rr))

library(ggplot2)
theme_set(theme_classic() +
            theme(aspect.ratio = 1))

ggplot(fitnesses, aes(birth_rate, resident_rr)) +
  geom_point(col="red") +
  geom_line(col="red") +
  geom_line(aes(birth_rate, y=1)) +
  scale_x_log10() +
  scale_y_log10()


q <- purrr::pluck(grid_search, 8)

x <- dplyr::group_by(q$resident_state, node) %>%
  dplyr::mutate(age = time - min(time),
                birth_year = min(time),
                density = exp(log_density)) %>%
  dplyr::ungroup()

ggplot(x, aes(time, density, group = node, color = birth_year)) +
  geom_line() +
  scale_colour_viridis_c()

ggplot(x, aes(time, height, group = node, colour = birth_year)) +
 geom_line() +
 scale_colour_viridis_c()

y <- x  %>%
 dplyr::select(-node, -log_density, -age, -birth_year) %>% stats::na.omit() %>%
 dplyr::filter(step > 1) %>% 
 dplyr::group_by(step, time) %>% 
 dplyr::summarise(
   density_integrated = -trapezium(height, density), 
   min_height = min(height),
   dplyr::across(where(is.double) & !c(density, density_integrated, min_height), 
                 ~-trapezium(height, density * .x)), 
   .groups = "drop"
 ) %>% 
 dplyr::rename(density = density_integrated)

ggplot(y, aes(time, density)) +
 geom_line()

# helper function - predicts to new values with spline
# only needed to ensure predictions for xout outside the ange of x are set to NA
f <- function(x, y, xout) {
 y_pred <- stats::spline(x, y, xout=xout)$y
 y_pred[!dplyr::between(xout, min(x), max(x))] <- NA
 y_pred
}

heights = c(1, 2, 8, 10, 16)

z <- x %>%
 dplyr::group_by(time, step) %>%
 dplyr::summarise(
   dplyr::across(where(is.double), ~f(height, .x, xout=heights)),
   .groups = "keep") %>%
 dplyr::mutate(density = exp(log_density),
               node = seq_len(dplyr::n())) %>%
 dplyr::ungroup()

ggplot(z, aes(time, density, group = node, colour = factor(height))) +
 geom_line()

