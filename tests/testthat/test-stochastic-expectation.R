## The deterministic solver is the expectation of the stochastic one, and where
## nothing competes that expectation can be had from a single plant. Every
## plant then sees the same environment, so a plant born at time s is, at time
## T, the plant grown alone to age T - s. With a constant birth rate b per unit
## area, the expected stand at T per unit area is
##
##   plants alive     b p_est \int_0^T S(a) da
##   total height     b p_est \int_0^T S(a) h(a) da
##   seed produced    S_D b p_est \int_0^T \int_0^a f(u) S(u) du da
##
## with p_est the establishment probability, S(a) = exp(-H(a)) survival on the
## plant's own hazard, h its height and f its fecundity rate. The seed counts
## the dead, with what they made before dying. All three are read off one plant
## grown on a fine grid of ages, independently of both solvers.
##
## Both solvers are held to that reference. The deterministic one, on a dense
## schedule so that reading the stand off its cohorts is accurate enough to
## test, to a relative 1e-3 in either density coordinate, each quantity on its
## own. (One expect_equal() over the three would scale its tolerance by their
## mean, which the seed output dominates, and pass a density several-fold out.)
## The stochastic one to within four standard errors, which follow from the
## process rather than from repeated runs: arrivals are a Poisson process and,
## with nothing competing, each plant establishes and survives independently of
## the rest, so each total is a compound Poisson sum, whose variance the sum of
## the squared terms estimates without bias. With a few hundred plants those
## errors are several per cent, so the stochastic comparison catches gross
## errors only; the deterministic one is the precise check.
##
## Surviving every death check is exactly exp(-H), however often the checks
## come, so the stand at a check time carries no discretisation error. The run's
## end is a check time, which the survivors' hazards, reset to zero there, show
## directly. A run with checks 2.5 yr apart is compared too, though at these
## errors it would catch only a gross dependence on the interval. Seed is not
## compared there, because a plant that dies goes on accruing seed until the
## check that finds it.

## A strategy whose plants do not affect one another, at birth rate `b`.
## FF16 casts no shade with k_I = 0. K93 cannot take that, since it recovers
## the basal area above a plant from the light it sees by dividing by k_I, so
## its competition coefficients are zeroed instead; its hazard is then
## max(-c_0, 0), and a negative c_0 gives it one to test.
uncompeting_strategy <- function(x, b) {
  s <- switch(x, FF16 = FF16_Model(), K93 = K93_Model())
  pars <- s$pars
  if (x == "FF16") {
    pars$k_I <- 0
  } else {
    pars$b_2 <- 0
    pars$c_1 <- 0
    pars$d_1 <- 0
    pars$c_0 <- -0.02
  }
  s$pars <- pars
  s$birth_rate_y <- b
  s
}

trapezium_rule <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

## The expected stand per unit area at `t_max`, from one plant.
expected_stand <- function(x, s, t_max, b, da = 0.02) {
  ind <- Individual(x, environment_type(x))(s)
  p_est <- ind$establishment_probability(Environment(x))
  a <- seq(0, t_max, by = da)
  state <- run(ind, times = a)$state
  surv <- exp(-state[, "mortality"])
  ## Seed made by an established plant up to each age, its survival included.
  seed_to_age <- c(0, cumsum(diff(state[, "fecundity"]) *
                               (head(surv, -1) + tail(surv, -1)) / 2))
  c(density = b * p_est * trapezium_rule(a, surv),
    height = b * p_est * trapezium_rule(a, surv * state[, "height"]),
    seed = s$pars$S_D * b * p_est * trapezium_rule(a, seed_to_age))
}

uncompeting_parameters <- function(x, s, t_max) {
  p <- base_parameters(x)
  p$strategies <- list(s)
  p$max_patch_lifetime <- t_max
  ## Lifetime seed output with no weighting by patch age, which is the seed
  ## produced per unit area over the run.
  p$fitness_weighting <- "none"
  p
}

## The deterministic stand at the end of the run, from cohorts every 0.1 yr.
## The schedule stops one interval short of the end, since a cohort introduced
## at the end has no step of its own, and its birth date could not be
## recovered from the output.
deterministic_stand <- function(x, s, t_max, coordinate) {
  p <- uncompeting_parameters(x, s, t_max)
  p$node_schedule_times <- list(seq(0, t_max - 0.1, by = 0.1))
  ctrl <- Control()
  ctrl$node_density_coordinate <- coordinate
  res <- run(p, control = ctrl, record = "steps")
  stand <- tail(integrate_over_size_distribution(tidy(res)$species,
                                                 coordinate = coordinate), 1)
  expect_equal(stand$time, t_max)
  c(density = stand$density, height = stand$height,
    seed = res$offspring_production)
}

## The stochastic stand at the end of the run, and its standard errors.
stochastic_stand <- function(x, s, t_max, area, death_interval, seed) {
  p <- uncompeting_parameters(x, s, t_max)
  ctrl <- Control(stochastic_patch_area = area)
  ctrl$stochastic_death_interval <- death_interval
  set.seed(seed)
  res <- run(p, control = ctrl, method = "stochastic")
  expect_equal(res$time, t_max)
  plants <- res$solver$patch$state$species[[1]]
  alive <- attr(plants, "is_alive")
  expect_true(all(plants["mortality", alive] == 0))
  height <- plants["height", alive]
  seed_made <- s$pars$S_D * plants["fecundity", ]
  list(n = length(alive),
       value = c(density = sum(alive), height = sum(height),
                 seed = sum(seed_made)) / area,
       se = c(density = sqrt(sum(alive)), height = sqrt(sum(height^2)),
              seed = sqrt(sum(seed_made^2))) / area)
}

test_that("without competition both solvers are the one-plant expectation", {
  t_max <- 25
  b <- 2
  for (x in c("FF16", "K93")) {
    s <- uncompeting_strategy(x, b)
    expected <- expected_stand(x, s, t_max, b)
    ## Plants are still dying and making seed at the end, so the comparison
    ## covers all three processes.
    expect_gt(expected[["seed"]], 0)
    expect_lt(expected[["density"]], 0.95 * b * t_max)

    for (coordinate in c("height", "birth_date")) {
      det <- deterministic_stand(x, s, t_max, coordinate)
      for (q in names(expected)) {
        expect_equal(det[[q]], expected[[q]], tolerance = 1e-3,
                     info = paste(x, coordinate, q))
      }
    }

    sto <- stochastic_stand(x, s, t_max, area = 10, death_interval = 0.05,
                            seed = 1)
    ## Enough plants for the standard errors to be meaningful.
    expect_gt(sto$n, 300)
    z <- (sto$value - expected) / sto$se
    expect_true(all(abs(z) < 4), info = paste(x, toString(signif(z, 3))))

    ## Checks every 2.5 yr, of which the end is one.
    coarse <- stochastic_stand(x, s, t_max, area = 10, death_interval = 2.5,
                               seed = 2)
    keep <- c("density", "height")
    z <- (coarse$value[keep] - expected[keep]) / coarse$se[keep]
    expect_true(all(abs(z) < 4), info = paste(x, toString(signif(z, 3))))
  }
})
