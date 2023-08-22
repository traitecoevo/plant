devtools::load_all(".")

p0 <- scm_base_parameters("FF16w")
p0$max_patch_lifetime <- 30
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,
                        birth_rate_list = list(10))


# init rainfall spline for env
x <- seq(0, 30, length.out = 1000)
rain = list(
  x = x,
  y = 2 + sin(x)
)

env <- make_environment("FF16w",
                        soil_number_of_depths = 1,
                        soil_initial_state = rep(0.4),
                        rainfall = rain)

ctrl <- scm_base_control()
ctrl$ci_niter = 1000
ctrl$GSS_tol_abs <- 1e-1
ctrl$ci_abs_tol <- 1e-8
ctrl$vulnerability_curve_ncontrol <- 100

out <- run_scm_collect(p1, env, ctrl,
                       collect_auxiliary_variables=TRUE)

