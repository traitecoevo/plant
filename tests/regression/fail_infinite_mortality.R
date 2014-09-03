library(tree)
disturbance <- 5
slope <- 2.5

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$disturbance <- new(Disturbance, disturbance)
p0$set_control_parameters(list(schedule_verbose=TRUE))

s <- new(Strategy, list(c_r1=0.5,
                        c_r2=0,
                        B4=slope))
p0$strategy_default <- s

sys <- list()
sys$traits <- cbind(lma=c(0.1978791, 1.32914554799473,
                      0.0258195785945316, 0.000225482844596086))
sys$seed_rain <- c(629.57180223657, 0.001, 0.001, 0.001)

p <- expand_parameters(colnames(sys$traits), sys$traits[,1], p0,
                       mutant=FALSE)
p$seed_rain <- sys$seed_rain

schedule <- default_cohort_schedule(p)

## From ba6c81b93f1c93a4b12d0b4396985afcb5b108ce and before, this will
## fail:
res <- build_schedule(p)

## The reason why is this:
ebt <- new(EBT, p$copy())
ebt$cohort_schedule <- schedule$copy()
ebt$run_next()
ebt$run_next() # will fail

ebt <- new(EBT, p$copy())
ebt$cohort_schedule <- schedule$copy()
patch <- ebt$patch

for (i in 4:1) {
  patch$add_seedling(i)
}
m <- matrix(patch$ode_values, 4)
rownames(m) <- c("height", "mortality", "fecundity", "log_density")

r <- matrix(patch$ode_rates, 4)
dimnames(r) <- dimnames(m)

## The infinite mortality rate here is the problem.  The new version
## stamps that down.
