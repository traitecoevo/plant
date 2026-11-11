source("helpers.R")
src <- readLines("diag.R"); eval(parse(text = src[2:19]))  # diag_run only
setp <- function(f) function(p) { for (i in seq_along(p$strategies)) p$strategies[[i]]$pars <- f(p$strategies[[i]]$pars); p }
th0 <- scm_base_parameters("TF24")$strategy_default$pars$theta
r <- diag_run("theta = plant default 2.14e-4", p_mod = setp(function(x) { x$theta <- th0; x }))
diag_run("theta default, hmat 8", p_mod = setp(function(x) { x$theta <- th0; x$hmat <- 8; x }))
diag_run("AusTraits theta, hmat 8", p_mod = setp(function(x) { x$hmat <- 8; x }))
# soil water in the baseline sparse stand
p <- make_p(list(pulse_br(0.002, 0, 1, 40)), 40); env <- make_env(MET0, 40)
rr <- run_scm(p, env, control(), collect = TRUE)
e <- rr$env %>% filter(!is.na(soil_moist)) 
print(e %>% group_by(soil_depth) %>% summarise(theta_mean = mean(soil_moist), theta_min = min(soil_moist)) , n = 20)
cat("soil_moist_sat", env$soil_moist_sat, "\n")
