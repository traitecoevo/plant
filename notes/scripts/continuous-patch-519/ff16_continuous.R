## A truly continuous FF16 patch (#519, the issue's original scope): one patch,
## no disturbance, constant light, a cohort introduced every year, run until the
## stand stops changing. No history is recorded; the run keeps one light
## snapshot every few years and the cost is timed per block.
##
##   PLANT_519_BIRTH    birth rate (1)
##   PLANT_519_TMAX     cap on the run (3000)
##   PLANT_519_BLOCK    years per convergence block (50)
##   PLANT_519_TOL      relative change between block means that counts as settled (0.005)
##   PLANT_519_ABSORB   control$absorb_epsilon, the absorbing top layer (0.02)
##   PLANT_519_OUT      output directory

suppressPackageStartupMessages(library(plant))
num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
birth  <- num_env("PLANT_519_BIRTH", 1)
t_max  <- num_env("PLANT_519_TMAX", 3000)
block  <- num_env("PLANT_519_BLOCK", 50)
tol    <- num_env("PLANT_519_TOL", 0.005)
absorb <- num_env("PLANT_519_ABSORB", 0.02)
out_dir <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
stamp("plant %s; FF16 lma 0.0825, birth rate %g, cap %g yr, block %g yr, tol %g, absorb_epsilon %g",
      as.character(packageVersion("plant")), birth, t_max, block, tol, absorb)

p <- scm_base_parameters("FF16")
p$patch_type <- "fixed"
p$max_patch_lifetime <- t_max
p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = birth)
p$node_schedule_times <- list(node_schedule_times_periodic(0, t_max, 1))
ctrl <- Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb)
scm <- plant:::SCM("FF16", "FF16_Env")(p, Environment("FF16"), plant:::empty_events(), ctrl)
scm$snapshot_interval <- 5

rows <- list()
previous <- NULL
streak <- 0L
converged <- FALSE
k <- 0L
repeat {
  k <- k + 1L
  t_end <- k * block
  if (t_end > t_max) break
  w <- system.time(scm$run_to(t_end))[["elapsed"]]
  sn <- tidy_snapshots(scm)$species
  in_block <- sn$time > t_end - block & sn$time <= t_end
  m <- c(density = mean(sn$density[in_block]), height_max = mean(sn$height_max[in_block]),
         n_live = mean(sn$n_live[in_block]))
  rel <- if (is.null(previous)) rep(NA_real_, 3) else abs(m - previous) / previous
  sp <- scm$patch$species[[1]]
  rows[[k]] <- data.frame(time = scm$time, wall_s = w, s_per_year = w / block,
                          steps = length(scm$ode_times), nodes = sp$size, retired = sum(sp$retired),
                          live = sp$size - sum(sp$retired), t(m),
                          rel_density = rel[1], rel_height_max = rel[2], rel_n_live = rel[3])
  if (k %% 4 == 0 || k <= 2) {
    stamp("%5.0f yr: %.2f s/yr, %d cohorts, %d live, %d retired; block means density %.4g, height_max %.4g, live %.1f; rel change %s",
          scm$time, w / block, sp$size, sp$size - sum(sp$retired), sum(sp$retired), m[1], m[2], m[3],
          paste(signif(rel, 2), collapse = "/"))
  }
  streak <- if (!is.null(previous) && all(rel < tol)) streak + 1L else 0L
  previous <- m
  if (streak >= 2L) { converged <- TRUE; break }
}
d <- do.call(rbind, rows)
stamp("%s at %g yr; total wall %.0f s", if (converged) "CONVERGED" else "not converged", scm$time, sum(d$wall_s))

## Fixed point or cycle? Carry on past convergence and look at the snapshots
## from the converged time on.
t_conv <- scm$time
t_verify <- num_env("PLANT_519_VERIFY", 1000)
scm$run_to(t_verify)
sn <- tidy_snapshots(scm)$species
stat <- sn[sn$time >= t_conv, ]
stamp("from %g to %g yr: density mean %.6g, CV %.2e; tallest live cohort %.4f to %.4f m; live cohorts %d to %d",
      t_conv, scm$time, mean(stat$density), sd(stat$density) / mean(stat$density),
      min(stat$height_max), max(stat$height_max), min(stat$n_live), max(stat$n_live))

## Where the density sits: by cohort age.
sp <- scm$patch$species[[1]]
age <- scm$time - sp$node_times
nu <- exp(sp$log_densities_state)
live <- nu > 0
stamp("at %g yr: %d live cohorts, ages %g to %g; share of density in cohorts older than 50 / 100 yr: %.2e / %.2e; tallest live %.3f m",
      scm$time, sum(live), min(age[live]), max(age[live]), sum(nu[age > 50]) / sum(nu),
      sum(nu[age > 100]) / sum(nu), max(sp$heights[live]))
write.csv(d, file.path(out_dir, "ff16-continuous.csv"), row.names = FALSE)
write.csv(sn, file.path(out_dir, "ff16-continuous-snapshots.csv"), row.names = FALSE)

## What the top layer buys: the same 3000 years with and without it.
for (r in c(absorb, 0)) {
  q <- p
  scm_r <- plant:::SCM("FF16", "FF16_Env")(q, Environment("FF16"), plant:::empty_events(),
                                           Control(node_density_in_birth_date = TRUE, absorb_epsilon = r))
  w <- system.time(scm_r$run())[["elapsed"]]
  stamp("3000 yr with absorb_epsilon %g: %.1f s, %d cohorts, %d still integrated, %d steps",
        r, w, scm_r$patch$species[[1]]$size, scm_r$patch$species[[1]]$size - sum(scm_r$patch$species[[1]]$retired),
        length(scm_r$ode_times))
}

## The equilibrium birth rate: in a stationary deterministic stand every arrival
## meets the same stand, so a short window after the stand has settled, run on
## long enough for its arrivals to finish, gives R exactly.
T1 <- num_env("PLANT_519_T1", 400)
W <- num_env("PLANT_519_WINDOW", 10)
L <- num_env("PLANT_519_LIFESPAN", 200)
R_at <- function(b) {
  q <- scm_base_parameters("FF16")
  q$max_patch_lifetime <- T1 + W + L
  q <- add_strategies(q, trait_matrix(0.0825, "lma"), birth_rate = b)
  q <- window_regime(q, T1, T1 + W)
  q$node_schedule_times <- list(node_schedule_times_periodic(0, T1 + W + L, 1))
  s <- plant:::SCM("FF16", "FF16_Env")(q, Environment("FF16"), plant:::empty_events(), ctrl)
  s$run()
  R <- s$net_reproduction_ratios
  sp <- s$patch$species[[1]]
  out <- sp$net_reproduction_ratio_by_node[sp$node_times >= T1 & sp$node_times <= T1 + W] * q$strategies[[1]]$pars$S_D
  stamp("b = %.6g: window R %.6g, the window's arrivals spread %.2e (max/min - 1)", b, R, max(out) / min(out) - 1)
  R
}
## R at stationarity scales as 1/b here (seed production per patch does not
## depend on the birth rate), so b* is read off two evaluations and checked at
## b*. The check is guarded: a stand that dense can stop on FF16's density
## overflow guard, which is reported rather than fatal.
## Birth rates of 10 and above stop on FF16's density overflow in this stand
## on develop too (birth-date coordinate, annual cohorts, no disturbance), so
## the elasticity is read over [b, 3b].
R1 <- R_at(birth)
R3 <- R_at(3 * birth)
elasticity <- log(R3 / R1) / log(3)
b_star <- birth * R1^(-1 / elasticity)
stamp("d log R / d log b = %.5f over [%g, %g]; b* = %.6g", elasticity, birth, 3 * birth, b_star)
check <- tryCatch(sprintf("R(b*) = %.6g", R_at(b_star)), error = function(e) paste("R(b*) not evaluated:", substr(conditionMessage(e), 1, 120)))
stamp("%s", check)
write.csv(data.frame(elasticity = elasticity, b_star = b_star, check = check),
          file.path(out_dir, "ff16-equilibrium.csv"), row.names = FALSE)
stamp("done")
