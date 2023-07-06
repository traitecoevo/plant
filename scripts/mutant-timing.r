library(dplyr)
library(purrr)
devtools::load_all()

# basic setup 
p <- scm_base_parameters("FF16")
p$max_patch_lifetime <- 100

create_resident <- function(p1, cache = F) {
  
  e <- make_environment("FF16")
  
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = cache
  
  types <- extract_RcppR6_template_types(p1, "Parameters")
  scm <- do.call('SCM', types)(p1, e, ctrl)
  scm$run()
  
  return(scm)
}

p1 <- expand_parameters(trait_matrix(0.1978791, "lma"), p,
  mutant = F,
  birth_rate_list = 40
)

# timing for resident only without caching
system.time(create_resident(p1))

# timing for resident only with caching
system.time(scm <- create_resident(p1, cache = T))

# timing for mutant re-using node schedule
system.time(scm$run_mutant(p1$strategies, append = F, update_schedule = F))
  
# timing for mutant with one node per ode step
# note - overwrites nodeschedule, so recreate to go back to default
system.time(scm$run_mutant(p1$strategies, append = F, update_schedule = T))  

# iterating one mutant at a time with original schedule
scm <- create_resident(p1, cache = T)
n = 10
lmas <- seq(0.17, 0.22, len = n)

lma_strategies_single <- purrr::map(lmas, 
  ~expand_parameters(trait_matrix(.x, "lma"), p, mutant = F, birth_rate_list = 1)$strategies)

system.time({
  for(i in 1:n) {
    scm$run_mutant(lma_strategies_single[[i]], append = F, update_schedule = F)
  }
})

# iterating one mutant at a time with one node per ode step
scm <- create_resident(p1, cache = T)
system.time({
  out_single <- purrr::map_df(lma_strategies_single, 
    ~ {
    scm$run_mutant(.x, append = F, update_schedule = T)
    tibble(lma = .x[[1]]$lma,
      rr = scm$net_reproduction_ratios)
    }, .id="lma_grid")
  })

# running many mutants at once
# note - mutants appear to interact, cause unknown, so only demonstrating runtime
lma_strategies_many <- 
  expand_parameters(trait_matrix(lmas, "lma"), p, mutant = F, birth_rate_list = rep(1, n))$strategies

# must update schedule to run many at once
system.time({
  scm$run_mutant(lma_strategies_many, append = F, update_schedule = T)
  out_many <- tibble(
    lma = map_dbl(lma_strategies_many, ~.x$lma),
      rr_many = scm$net_reproduction_ratios)
  })

# do we get same results using 1 at a time vs many
left_join(out_single, out_many) %>%
mutate(diff = rr - rr_many)

# but we never set the birth_rates for the mutants! I think that may be causing the issues with many mutants case

# What is speed cost of caching? Run same resident without

# why is mutant as slow as resident? Loading cache??? Use pointers
n <- 5
birth_total <- 40

p1 <- expand_parameters(trait_matrix(0.1978791, "lma"), p,
  mutant = F,
  birth_rate_list = birth_total
)

p1_diff <- expand_parameters(trait_matrix(0.25, "lma"), p,
  mutant = F,
  birth_rate_list = 1
)

scm1 <- create_resident(p1, cache = T)
scm1$net_reproduction_ratios
scm1$run_mutant(p1$strategies, append = F, update_schedule = F)
scm1$net_reproduction_ratios

scm1$run_mutant(p1_diff$strategies, append = F, update_schedule = F)
scm1$net_reproduction_ratios

pn <- expand_parameters(trait_matrix(rep(0.1978791, n), "lma"), p,
  mutant = F,
  birth_rate_list = rep(birth_total, n)/n
)

scmn <- create_resident(pn, cache = T)
scmn$net_reproduction_ratios

## Works!!!

scmn$run_mutant(pn$strategies, append = F, update_schedule = F)
scmn$net_reproduction_ratios

## Works!!!

scmn$run_mutant(p1_diff$strategies, append = F, update_schedule = F)
scmn$net_reproduction_ratios

## Fails-  memory allocation fault
scmn$run_mutant(pn$strategies, append = F, update_schedule = F)
scmn$net_reproduction_ratios


## Same number of strategies --
pn_diff <- expand_parameters(trait_matrix(c(rep(0.25, n)), "lma"), p,
  mutant = F,
  birth_rate_list = rep(1, n)
)

# Works
scmn$run_mutant(pn_diff$strategies, append = F, update_schedule = F)
scmn$net_reproduction_ratios

## Extra number of strategies -- Try adding 1 more strategy, predict this fails
## Fails -- Error: Need at least two points for the trapezium rule
pn_diff <- expand_parameters(trait_matrix(c(rep(0.25, n+1)), "lma"), p,
  mutant = F,
  birth_rate_list = rep(1, n+1)
)

scmn$run_mutant(pn_diff$strategies, append = F, update_schedule = F)
scmn$net_reproduction_ratios

