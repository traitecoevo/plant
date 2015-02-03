## Build a single schedule:
p <- ebt_base_parameters()
p$strategies <- list(Strategy())
p$seed_rain <- 1.1
p <- build_schedule(p)
