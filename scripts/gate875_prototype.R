## #875 gate: one-lane forward tangent of an invader replay vs the reverse sweep,
## on regnans's FF16 configuration (lifetime 30, lma-hmat attractor).
suppressMessages(library(plant))
cat("plant from", find.package("plant"), "\n")
p <- default_parameters("FF16", quiet = TRUE)
p$strategy_default$pars$a_l1 <- 2.17
p$strategy_default$pars$a_l2 <- 0.5
p$strategy_default$pars$hmat <- 10.9
p$max_patch_lifetime <- 30
p$node_schedule_times_default <- node_schedule_times_default(30)
ctrl <- Control(node_density_coordinate = "birth_date")
p <- add_strategies(p, trait_matrix(c(0.0473, 10.9), c("lma", "hmat")), birth_rate = 20)
p$node_schedule_times <- list(p$node_schedule_times_default)
res <- run(p, control = ctrl)
solver <- res$solver
mut <- add_mutant(p, trait_matrix(c(0.0473, 10.9), c("lma", "hmat")))
el <- function(f, n = 5) { t <- numeric(n); for (i in seq_len(n)) t[i] <- system.time(f())[["elapsed"]]; min(t) }
t_run <- el(function() run(p, control = ctrl))
t_mut <- el(function() solver$run_mutant(mut))
solver$run_mutant(mut)
metric <- "1.offspring_production"
g <- stand_gradient(solver, metrics = metric, hyperpar = NULL)
cols <- colnames(g$gradient)
cat("columns:", length(cols), "\n")
traits_all <- c("lma", "hmat", "rho", "omega", "a_l2", "a_f1", "theta", "k_l")
t_sweep <- sapply(c(1, 2, 4, 8), function(k)
  el(function() stand_gradient(solver, metrics = metric, traits = paste0("1.", traits_all[seq_len(k)]), hyperpar = NULL)))
dir <- as.numeric(cols == "1.lma")
t_tan <- el(function() plant:::census_trait_tangent_prototype(solver, "FF16", dir))
tan <- plant:::census_trait_tangent_prototype(solver, "FF16", dir)
mi <- which(names(solver$patch$offspring_production) == "" | TRUE)
cat(sprintf("resident run %.3f s; invader replay %.3f s\n", t_run, t_mut))
cat("reverse sweep k=1,2,4,8:", sprintf("%.3f", t_sweep), "s\n")
cat(sprintf("forward tangent, one lane: %.3f s (= %.2f replays; half sweep = %.3f s)\n",
            t_tan, t_tan / t_mut, min(t_sweep) / 2))
cat("tangent per metric:", format(tan$tangent, digits = 8), "\n")
cat("sweep lma column:  ", format(g$gradient[, "1.lma"], digits = 8), "\n")
cat("\n-- multi-lane forward --\n")
gall <- stand_gradient(solver, hyperpar = NULL)$gradient
for (k in c(1, 2, 4, 8)) {
  tr <- paste0("1.", traits_all[seq_len(k)])
  D <- sapply(tr, function(t) as.numeric(cols == t))
  if (k == 1) D <- matrix(D, ncol = 1)
  t_l <- el(function() plant:::census_trait_tangent_lanes_ff16(solver, D), n = 7)
  out <- plant:::census_trait_tangent_lanes_ff16(solver, D)
  rownames(out) <- rownames(gall)
  err <- max(abs(out - gall[, tr, drop = FALSE]) / pmax(abs(gall[, tr, drop = FALSE]), 1e-300))
  cat(sprintf("k=%d: forward %.3f s (%.2f replays) vs sweep %.3f s -> ratio %.2f; max rel diff vs sweep %.1e\n",
              k, t_l, t_l / t_mut, t_sweep[match(k, c(1, 2, 4, 8))], t_l / t_sweep[match(k, c(1, 2, 4, 8))], err))
}
cat("\n-- central difference of the invader's fitness (2 replays per trait) --\n")
base_tr <- c(lma = 0.0473, hmat = 10.9)
fd_k <- function(k) {
  tr <- traits_all[seq_len(k)]
  vals <- sapply(tr, function(t) p$strategies[[1]]$pars[[t]])
  for (i in seq_along(tr)) for (s in c(-1, 1)) {
    m <- add_mutant(p, trait_matrix(base_tr, names(base_tr)))
    st <- m$strategies[[length(m$strategies)]]; st$pars[[tr[i]]] <- vals[[i]] * (1 + s * 1e-4)
    m$strategies[[length(m$strategies)]] <- st
    solver$run_mutant(m); stand_fitness(solver)
  }
}
for (k in c(1, 2, 4, 8)) cat(sprintf("k=%d: FD %.3f s\n", k, el(function() fd_k(k))))
