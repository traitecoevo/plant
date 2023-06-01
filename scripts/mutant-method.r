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

run_one_patch <- function(birth_rate, logfile = "logs.txt") {

  sink(logfile)
  p <- scm_base_parameters("FF16")

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
    dplyr::left_join(time)
    
  resident_rr <- scm$net_reproduction_ratios

  # use stateful SCM for efficient mutant runs
  scm$run_mutant()
  
  mutant_rr <- scm$net_reproduction_ratios
  
  mutant_endpoint <- scm$state$species %>% 
    purrr::pluck(., 1) %>%
    t() %>% tibble::as_tibble() %>%
    dplyr::mutate(time = scm$time, 
                  node = seq_along(time))
  sink()
  return(list(fitnesses = c(resident_rr, mutant_rr),
              resident_state = state,
              mutant_endpoint = mutant_endpoint))
}

birth_rates <- c(100, 50, 10, 5, 1, 0.1, 0.01, 0.001)

grid_search <- purrr::map(birth_rates[1], run_one_patch)


fitnesses <- purrr::map_df(grid_search,
                           ~ purrr::pluck(., "fitnesses") %>% 
                             tibble::tibble(., rownames = c("resident_rr", "mutant_rr"),
                                            .name_repair = ~ c("net_reproduction_ratio", "invasion_type")),
                           .id = "grid_step")

tidyr::spread(fitnesses, invasion_type, net_reproduction_ratio) %>%
  dplyr::mutate(birth_rate = birth_rates,
                unit = "seeds_per_m2_yr") %>%
  dplyr::group_by(grid_step, birth_rate, unit) %>%
  dplyr::mutate(prop = scales::percent(mutant_rr / resident_rr))

library(ggplot2)
theme_set(theme_classic() +
            theme(aspect.ratio = 1))

purrr::pluck(grid_search, 1) %>%
  {
    ggplot(.$resident_state, aes(time, exp(log_density), group = node)) +
    geom_line() +
    geom_point(data = .$mutant_endpoint, color = "red")
  
  ggplot(.$resident_state, aes(time, height, group = node)) +
    geom_line() +
    geom_point(data = .$mutant_endpoint, color = "red")
  }
