## Build a single schedule:
plant_log_console()
p <- scm_base_parameters()
p$strategies <- list(FF16_Strategy())
p$birth_rate <- 1.1
p <- build_schedule(p)
