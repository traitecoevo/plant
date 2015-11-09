## Build a single schedule:
plant_log_console()
p <- ebt_base_parameters()
p$strategies <- list(FF16_Strategy())
p$seed_rain <- 1.1
p <- build_schedule(p)
