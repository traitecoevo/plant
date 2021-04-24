## Find the equilibrium seed rain a few different ways:
plant_log_console()
p <- scm_base_parameters()
p$strategies <- list(FF16_Strategy())
p$offspring_arriving <- 1.1

p$control$equilibrium_solver_name <- "iteration"
ans <- equilibrium_offspring_arriving(p)
ans$offspring_arriving

p$control$equilibrium_solver_name <- "nleqslv"
ans <- equilibrium_offspring_arriving(p)
ans$offspring_arriving

p$control$equilibrium_solver_name <- "dfsane"
ans <- equilibrium_offspring_arriving(p)
ans$offspring_arriving

p$control$equilibrium_solver_name <- "hybrid"
ans <- equilibrium_offspring_arriving(p)
ans$offspring_arriving

bounds <- viable_fitness(bounds(lma=c(0.01, 10.0)), p)
traits <- cbind(lma=seq_log_range(bounds, length.out=50))
w <- fitness_landscape(traits, ans)

plot(traits, w, type="l", log="x", xlab="LMA", ylab="Fitness", las=1)
abline(v=p$strategies[[1]]$lma, h=0, col="grey")
