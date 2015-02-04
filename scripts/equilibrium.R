## Find the equilibrium seed rain a few different ways:
p <- ebt_base_parameters()
p$strategies <- list(Strategy())
p$seed_rain <- 1.1

p$control$equilibrium_solver_name <- "iteration"
ans <- equilibrium_seed_rain(p)
ans$seed_rain

p$control$equilibrium_solver_name <- "nleqslv"
ans <- equilibrium_seed_rain(p)
ans$seed_rain

p$control$equilibrium_solver_name <- "dfsane"
ans <- equilibrium_seed_rain(p)
ans$seed_rain

p$control$equilibrium_solver_name <- "hybrid"
ans <- equilibrium_seed_rain(p)
ans$seed_rain

bounds <- viable_fitness(bounds(lma=c(0.01, 10.0)), p)
traits <- cbind(lma=seq_log_range(bounds, length=50))
w <- fitness_landscape(traits, ans)

plot(traits, w, type="l", log="x", xlab="LMA", ylab="Fitness", las=1)
abline(v=p$strategies[[1]]$lma, h=0, col="grey")
