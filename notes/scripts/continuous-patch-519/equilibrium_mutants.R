## Mutants at a resident's fixed point (#519): a continuous patch under constant
## drivers, no disturbance, a cohort a year and the absorbing top layer, run to
## its fixed point. There the resident makes one field, and each mutant stands
## in it, held, for its whole life (`fixed_point_fitness()`, which runs
## `SCM$run_mutant_held()`): nothing is recorded and nothing replayed.
##
## The script finds the resident's demographic equilibrium birth rate b* (R = 1
## by its own cohort in its own field; regnans owns this in production), reads
## the landscape there, and checks it three ways: the field held at four points
## inside one introduction interval, held a century either side of the reading
## year, and against the same cohort born into the running stand.
##
## ⚠️ THE INTRODUCTION GRID HAS TO BE FINE ENOUGH FOR A FIXED POINT. Between two
## introductions the stand's field swings, so a grid that is too coarse makes
## the stand a periodic orbit, and a field held at an introduction is the top
## of the swing. For FF16 at b* the swing in held R is 1.8 % at a cohort a year,
## 0.17 % at two and 0.02 % at four, where held and in-stand agree to 3e-5.
##
##   PLANT_519_MODEL     "FF16" (default) or "TF24"
##   PLANT_519_ABSORB    control$absorb_epsilon (0.02 for FF16, 0.05 for TF24)
##   PLANT_519_RAIN      TF24's constant rainfall, m per year (1)
##   PLANT_519_DT        years between introductions (0.25)
##   PLANT_519_T_EQ      year at which the field is read (400)
##   PLANT_519_LIFETIME  years each mutant runs; it must outlive its cohort (3000)
##   PLANT_519_LIFE      years the in-stand reference carries its cohorts (400)
##   PLANT_519_LMA       comma-separated mutant lma values
##   PLANT_519_OUT       output directory

suppressPackageStartupMessages(library(plant))
num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
vec_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(strsplit(v, ",")[[1]])
}
model    <- Sys.getenv("PLANT_519_MODEL", "FF16")
tf24     <- identical(model, "TF24")
absorb   <- num_env("PLANT_519_ABSORB", if (tf24) 0.05 else 0.02)
rain     <- num_env("PLANT_519_RAIN", 1)
dt       <- num_env("PLANT_519_DT", 0.25)
t_eq     <- num_env("PLANT_519_T_EQ", 400)
lifetime <- num_env("PLANT_519_LIFETIME", 3000)
life     <- num_env("PLANT_519_LIFE", 400)
lma_grid <- vec_env("PLANT_519_LMA", c(0.04, 0.05, 0.06, 0.07, 0.0825, 0.1, 0.125, 0.15, 0.2))
out_dir  <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
stamp("plant %s; %s, resident lma 0.0825, absorb_epsilon %g%s, a cohort every %g yr, field read at %g yr, mutants run %g yr",
      as.character(packageVersion("plant")), model, absorb,
      if (tf24) sprintf(", rainfall %g m/yr", rain) else "", dt, t_eq, lifetime)

resident_lma <- 0.0825
env_type <- if (tf24) "TF24_Env" else "FF16_Env"
## TF24's soil ODE is stiff for an explicit stepper without a cap on the step.
ctrl <- if (tf24) {
  Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb, ode_step_size_max = 0.05)
} else {
  Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb)
}
make_env <- function() {
  env <- Environment(model)
  if (tf24) env$extrinsic_drivers_set_constant("rainfall", rain)
  env
}
hyperpar <- if (tf24) TF24_hyperpar else FF16_hyperpar
add_lma <- function(p, lma, b) {
  add_strategies(p, trait_matrix(lma, "lma"), hyperpar = hyperpar, birth_rate = b)
}

## A cohort every `dt` years from time zero.
resident_parameters <- function(b, t_max) {
  p <- scm_base_parameters(model)
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- t_max
  p <- add_lma(p, resident_lma, b)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, t_max, dt))
  p
}

## ---- 1. the resident at its fixed point -------------------------------------
## Run to `t_read` and stopped there, whether or not it is an introduction: the
## run ends at its lifetime, which is how a field is read inside an interval.
resident_at <- function(b, t_read = t_eq) {
  p <- resident_parameters(b, t_read)
  scm <- plant:::SCM(model, env_type)(p, make_env(), plant:::empty_events(), ctrl)
  scm$run()
  list(scm = scm, p = p)
}
fitness <- function(r, lma) {
  fixed_point_fitness(r$scm, r$p, trait_matrix(lma, "lma"), hyperpar = hyperpar,
                      lifetime = lifetime)
}

## ---- 2. the reference: the same cohort born into the running stand ----------
## The resident introduced as above, plus the mutant mid-interval at a birth
## rate 1e-6 of the resident's. A window around it, holding no resident cohort,
## exempts it alone from absorption, and it is run out for `life` years.
simultaneous <- function(b, lma) {
  p <- resident_parameters(b, t_eq + life)
  p <- add_lma(p, lma, b * 1e-6)
  p$node_schedule_times[[2]] <- t_eq + 0.5 * dt
  p <- window_regime(p, t_eq + 0.25 * dt, t_eq + 0.75 * dt)
  w <- system.time(scm <- run_scm(p, env = make_env(), ctrl = ctrl))[["elapsed"]]
  sp <- scm$patch$species
  c(R_mutant = sp[[2]]$net_reproduction_ratio_by_node[[1]] * p$strategies[[2]]$pars$S_D, wall = w)
}

## ---- 3. the resident's demographic equilibrium ------------------------------
R_of_b <- function(b) fitness(resident_at(b), resident_lma)
## A bracket by factors of four from b = 1, then a root-find in log b.
R1 <- R_of_b(1)
step <- if (R1 > 1) log(4) else -log(4)
lb <- 0
repeat {
  lb_next <- lb + step
  if (sign(log(R_of_b(exp(lb_next)))) != sign(log(R1))) break
  lb <- lb_next
  if (abs(lb) > log(4^8)) stop("no equilibrium birth rate within a factor 4^8 of 1")
}
w_root <- system.time(root <- uniroot(function(lb) log(R_of_b(exp(lb))), sort(c(lb, lb_next)),
                                      tol = 1e-10))[["elapsed"]]
b_star <- exp(root$root)
stamp("R(b = 1) = %.6g; b* = %.8g after %d evaluations in %.2f s, each a resident run to %g yr and one held mutant",
      R1, b_star, root$iter, w_root, t_eq)

## ---- 4. the landscape at b*, every mutant in the one field ------------------
w_res <- system.time(r <- resident_at(b_star))[["elapsed"]]
sp <- r$scm$patch$species[[1]]
stamp("resident at b*: run to %g yr in %.3f s; %d nodes, %d integrated, top-layer density %.4g",
      t_eq, w_res, sp$size, sp$size - sum(sp$retired), sp$top_layer$density)
rows <- lapply(lma_grid, function(lma) {
  w <- system.time(R <- fitness(r, lma))[["elapsed"]]
  stamp("mutant lma %.4f: R = %.6g in %.3f s", lma, R, w)
  data.frame(lma = lma, R = R, wall_s = w)
})
d <- do.call(rbind, rows)
d$R_resident <- d$R[d$lma == resident_lma]

## ---- 5. inside one introduction interval, and drift -----------------------------
check_lma <- c(0.06, resident_lma, 0.1)
fractions <- c(0, 0.25, 0.5, 0.75)
ph <- sapply(fractions, function(f) {
  rr <- resident_at(b_star, t_eq + f * dt)
  sapply(check_lma, function(lma) fitness(rr, lma))
})
dimnames(ph) <- list(paste0("lma ", check_lma), paste0("+", fractions, " dt"))
stamp("R in the field held at four points inside one introduction interval:")
print(signif(ph, 8))
stamp("relative spread inside the interval: %s",
      paste(signif(apply(ph, 1, function(x) diff(range(x)) / mean(x)), 3), collapse = ", "))
drift <- sapply(c(-100, 0, 100), function(dt) fitness(resident_at(b_star, t_read = t_eq + dt), resident_lma))
stamp("the resident's R read at %g, %g and %g yr: %s", t_eq - 100, t_eq, t_eq + 100,
      paste(signif(drift, 7), collapse = ", "))

## ---- 6. against the full simultaneous run ------------------------------------
ref <- t(sapply(check_lma, function(lma) simultaneous(b_star, lma)))
ref <- data.frame(lma = check_lma, ref)
ref$R_field <- d$R[match(check_lma, d$lma)]
ref$rel_diff <- ref$R_field / ref$R_mutant - 1
stamp("against the full simultaneous run at b* (cohorts born at %g, run out %g yr):", t_eq, life)
print(format(ref, digits = 6), row.names = FALSE)

prefix <- paste0(tolower(model), "-equilibrium-")
write.csv(d, file.path(out_dir, paste0(prefix, "mutants.csv")), row.names = FALSE)
write.csv(ref, file.path(out_dir, paste0(prefix, "check.csv")), row.names = FALSE)
stamp("done")
