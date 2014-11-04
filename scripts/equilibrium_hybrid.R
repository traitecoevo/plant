library(tree)
library(tree.assembly)

make_pars <- function(pars, time_disturbance) {
  p <- ebt_base_parameters()
  p$set_control_parameters(list(equilibrium_nsteps=50,
                                equilibrium_eps=1e-3,
                                equilibrium_progress=TRUE))
  set_equilibrium_solver("hybrid", p)
  p$strategy_default <- new(Strategy, pars)
  p$disturbance <- new(Disturbance, time_disturbance)
  p
}

add_strategy <- function(p, pars) {
  s <- p$strategy_default$copy()
  s$set_parameters(pars)
  p$add_strategy(s)
}

time_disturbance <- 2.5
pars <- list(c_r1=0.5, c_r2=0, B5=0, c_d1=1, B6=1)

p <- make_pars(pars, time_disturbance)
add_strategy(p, list(rho=5.1))

res <- tree:::equilibrium_seed_rain2(p)

## Now, try on a hard case from our simulations.  This does take quite
## a while to run, but gets there in the end.
path <- "../../successional_diversity/analysis/experiments/output/lma_grid/simulations"
d <- readRDS(file.path(path, "32.rds"))
sys <- restore_community(d[[6]])
p <- sys$to_parameters()
set_equilibrium_solver("hybrid", p)
p$set_control_parameters(list(equilibrium_eps=1e-3,
                              equilibrium_nsteps=300))

res <- check_inviable(p)

## Then, here manually filter:
p2 <- p$copy()
p2$clear()
keep <- which(!attr(res, "drop"))
for (i in keep) {
  p2$add_strategy(p[[i]])
}
p2$seed_rain <- res[keep, "out"]

options(error=recover)
p2$set_control_parameters(list(equilibrium_nsteps=10)) # too low, but ok
set_equilibrium_solver("hybrid", p2)
res <- tree:::equilibrium_seed_rain2(p2)
