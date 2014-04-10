library(tree)

find_equilibirum <- function(time.disturbance){
  ## Try to establish what the equilubrium seed rain is.
  p <- new(Parameters)
  p$disturbance <- new(Disturbance, time.disturbance)
  traits=cbind(lma=0.1978791, hmat=16.5958691)
  p$add_strategy(new(Strategy, as.list(traits[1,])))
  p$seed_rain <- 1.1                       # Starting rain.
  p$set_parameters(list(patch_area=1.0))   # See issue #13
  p$set_control_parameters(fast.control()) # A bit faster

  t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
  times0 <- cohort.introduction.times(t.max)
  schedule0 <- schedule.from.times(times0)

  ## # 1: Approach to equilibrium:
  res <- equilibrium.seed.rain(p, schedule0, 10, progress=TRUE,
                               build.args=list(verbose=TRUE))
  list(time.disturbance=time.disturbance,
    approach=t(sapply(attr(res, "progress"), "[[", "seed_rain"))
    )
}
