## Start off with a empty landscape:
p <- ebt_base_parameters()
m <- trait_matrix(seq_log(0.01, 10.0, length.out=51), "lma")
w <- fitness_landscape(m, p)
plot(m, w, type="l", log="x", xlab="LMA", ylab="Fitness", las=1)

## Compute the maximum fitness:
bounds_lma <- bounds(lma=c(0.01, 10.0))
w_max <- max_fitness(bounds_lma, p)
plot(m, w, type="l", log="x", xlab="LMA", ylab="Fitness", las=1)
points(w_max, attr(w_max, "fitness"), pch=19)

## Compute viable bounds on fitness:
bounds_lma_w <- viable_fitness(bounds_lma, p)
plot(m, w, type="l", log="x", xlab="LMA", ylab="Fitness", las=1)
abline(v=bounds_lma_w, h=0, col="grey")

## In 2d (not yet implemented)
if (FALSE) {
  bounds_2d <- bounds(lma=c(0.01, 10.0),
                      rho=c(300.0, 2000.0))
  bounds_2d_w <- viable_fitness(bounds_2d, p)
  bounds_2d_w
}

## Infinite bounds:
p <- ebt_base_parameters()
viable_fitness(bounds_infinite("lma"), p)
