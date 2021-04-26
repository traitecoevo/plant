## Build a single schedule:
plant_log_console()
p <- scm_base_parameters()
p$strategies <- list(FF16_Strategy())
p$net_reproduction_ratios <- 1.1
p <- build_schedule(p)
