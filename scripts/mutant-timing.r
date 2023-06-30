devtools::load_all()

# basic setup 
create_resident <- function(cache = F) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 100
  
  p1 <- expand_parameters(trait_matrix(0.1978791, "lma"), p, mutant = F,
                          birth_rate_list = 40)
  
  e <- make_environment("FF16")
  
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = cache
  
  types <- extract_RcppR6_template_types(p1, "Parameters")
  scm <- do.call('SCM', types)(p1, e, ctrl)
  scm$run()
  
  return(scm)
}

# timing for resident only without caching
system.time(create_resident())

# timing for resident only with caching
system.time(scm <- create_resident(cache = T))

# timing for mutant re-using node schedule
s = FF16_Strategy()
system.time(scm$run_mutant(list(s), append = F, update_schedule = F))
  
# timing for mutant with one node per ode step
# note - overwrites nodeschedule, so recreate to go back to default
system.time(scm$run_mutant(list(s), append = F, update_schedule = T))  

# iterating one mutant at a time with original schedule
scm <- create_resident(cache = T)
n = 10
lmas <- seq(0.17, 0.22, len = n)

system.time({
  for(i in 1:n) {
    s$lma = lmas[i]
    scm$run_mutant(list(s), append = F, update_schedule = F)
  }
})

# iterating one mutant at a time with one node per ode step
scm <- create_resident(cache = T)
system.time({
  for(i in 1:n) {
    s$lma = lmas[i]
    scm$run_mutant(list(s), append = F, update_schedule = T)
  }
})

# running many mutants at once
# note - mutants appear to interact, cause unknown, so only demonstrating runtime
m = lapply(lmas, function(lma) {x = FF16_Strategy(); x$lma = lma; return(x)})

# must update schedule to run many at once
scm <- create_resident(cache = T)
system.time(scm$run_mutant(m, append = F, update_schedule = T))
