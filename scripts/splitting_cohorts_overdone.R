## Cases of overdone wood density.
library(tree)

make_parameters <- function(time_disturbance, B5=0.5) {
  strategy <- new(Strategy,
                  list(B5=B5, c_d1=0, B6=0, c_r1=0.5, c_r2=0))
  p <- tree::ebt_base_parameters()
  p$strategy_default <- strategy
  p$disturbance <- new(Disturbance, time_disturbance)

  ## Tweak seed rain equilibrium things here:
  set_equilibrium_solver("iteration", p)
  p$set_control_parameters(list(equilibrium_progress=FALSE))

  p
}

p <- make_parameters(16)
p$set_control_parameters(list(schedule_progress=TRUE))
s <- p$strategy_default$copy()
s$set_parameters(list(rho=50))
p$add_strategy(s)
p$seed_rain <- 474.9818
res <- build_schedule(p)

progress <- attr(res, "progress")

r <- sapply(progress, function(x) x$seed_rain[[2]])
n <- sapply(progress, function(x) x$schedule$size)

f <- function(m) {
  suppressWarnings(apply(m, 2, max, na.rm=TRUE))
}
total <- lapply(seq_len(p$size), function(idx)
                f(rbind(err_lai[[idx]], err_seed_rain[[idx]])))

eps <- p$control$parameters$schedule_eps

## Problems start about here:
x <- progress[[5]]

t <- x$schedule$times(1)
image(t, t, x$err$lai[[1]])

tmp <- run_ebt_collect(p, x$schedule$copy())

plot(tmp$light_env[[200]], ylim=c(0, 1), type="o", pch=19, cex=.4)
