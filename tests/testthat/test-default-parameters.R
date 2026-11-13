# Trade-offs travel with Parameters: `p$hyperpar` carries the function the
# strategies were built with, default_parameters() attaches plant's frozen
# testing trade-offs, and the bare base_parameters() attaches none.

tf24_derived <- c("k_l", "d_I", "TF24_cost_scale", "k_s", "r_s", "r_b", "a_f3",
                  "r_l", "stem_c", "stem_P50", "stem_b", "psi_crit")

test_that("base_parameters attaches no trade-offs: traits go straight to pars", {
  p <- base_parameters("TF24")
  expect_null(p$hyperpar)
  expect_true(plant:::is_identity_hyperpar(param_hyperpar(p)))
  # K_s is what TF24's trade-off turns into the stem curve; with none attached
  # it is a plain parameter and the curve stays at the default.
  s <- suppressMessages(generate_strategy(p, trait_matrix(0.5, "K_s")))[[1]]
  expect_equal(s$pars$K_s, 0.5)
  expect_identical(s$pars$stem_P50, TF24_Model()$pars$stem_P50)
})

test_that("generate_strategy says, once, that no trade-offs are attached", {
  # Make rlang's once-per-session message fire on every call for the test.
  rlang::local_options(rlib_message_verbosity = "verbose")
  p <- base_parameters("FF16")
  expect_message(generate_strategy(p, trait_matrix(0.1, "lma")),
                 "No trade-off function is attached")
  # The explicit identity is a choice, and is not remarked on.
  expect_silent(generate_strategy(p, trait_matrix(0.1, "lma"), hyperpar = NULL))
  expect_silent(generate_strategy(default_parameters("FF16", quiet = TRUE),
                                  trait_matrix(0.1, "lma")))
})

test_that("default_parameters warns that it is a testing set and attaches the frozen trade-offs", {
  expect_warning(p <- default_parameters("TF24"), "testing parameter set")
  expect_silent(p <- default_parameters("TF24", quiet = TRUE))
  expect_true(is.function(p$hyperpar))
  expect_identical(param_hyperpar(p), p$hyperpar)
  for (model in c("FF16", "TF24", "TF24f", "TF24t", "TF24w", "K93")) {
    p <- default_parameters(model, quiet = TRUE)
    expect_inherits(p, "Parameters")
    expect_true(is.function(p$hyperpar))
  }
  # K93 has no trade-offs: its testing set carries the identity.
  expect_true(plant:::is_identity_hyperpar(default_parameters("K93", quiet = TRUE)$hyperpar))
  expect_error(default_parameters("NOPE", quiet = TRUE), "Unknown type")
})

test_that("the attached trade-offs reproduce the explicit ones bit for bit", {
  m <- cbind(lma = c(0.1, 0.3), rho = c(600, 900), K_s = c(0.2, 0.8))
  a <- add_strategies(default_parameters("TF24", quiet = TRUE), m, birth_rate = c(1, 1))
  b <- add_strategies(base_parameters("TF24"), m, birth_rate = c(1, 1),
                      hyperpar = default_hyperpar("TF24"))
  expect_identical(a$strategies, b$strategies)
  expect_true(all(tf24_derived %in% names(a$strategies[[1]]$pars)))
  expect_false(identical(a$strategies[[1]]$pars$stem_P50, TF24_Model()$pars$stem_P50))

  m <- cbind(lma = c(0.1, 0.3), narea = c(1.5e-3, 2.5e-3))
  a <- add_strategies(default_parameters("FF16", quiet = TRUE), m, birth_rate = c(1, 1))
  b <- add_strategies(base_parameters("FF16"), m, birth_rate = c(1, 1),
                      hyperpar = default_hyperpar("FF16"))
  expect_identical(a$strategies, b$strategies)
  # The TF24 variants share TF24's trade-offs.
  for (model in c("TF24f", "TF24t", "TF24w")) {
    s <- TF24_Model()
    expect_identical(default_hyperpar(model)(m[, "lma", drop = FALSE], s, filter = FALSE),
                     default_hyperpar("TF24")(m[, "lma", drop = FALSE], s, filter = FALSE))
  }
})

test_that("an explicit hyperpar still overrides the one Parameters carry", {
  p <- default_parameters("TF24", quiet = TRUE)
  s <- generate_strategy(p, trait_matrix(0.5, "K_s"), hyperpar = NULL)[[1]]
  expect_identical(s$pars$stem_P50, TF24_Model()$pars$stem_P50)
  s <- generate_strategy(p, trait_matrix(0.5, "K_s"))[[1]]
  expect_false(identical(s$pars$stem_P50, TF24_Model()$pars$stem_P50))
})

test_that("the trade-offs survive validate(), run() and a saved object", {
  p <- default_parameters("FF16", quiet = TRUE)
  p$max_patch_lifetime <- 3
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  f <- p$hyperpar
  expect_identical(plant:::validate(p)$hyperpar, f)

  scm <- run(p)$solver
  expect_identical(scm$parameters$hyperpar, f)
  expect_identical(param_hyperpar(scm$parameters), f)
  res <- tidy(run(p, record = "steps"))
  expect_identical(res$p$hyperpar, f)

  m <- trait_matrix(c(0.1, 0.2), "lma")
  s <- p$strategy_default
  q <- unserialize(serialize(p, NULL))
  expect_identical(q$hyperpar(m, s, filter = FALSE), f(m, s, filter = FALSE))

  # Removing the element is the same as carrying NULL, through the round trip.
  p$hyperpar <- NULL
  expect_false("hyperpar" %in% names(p))
  expect_null(plant:::validate(p)$hyperpar)
  expect_true(plant:::is_identity_hyperpar(param_hyperpar(plant:::validate(p))))

  # The constructor takes it as a field, and refuses a non-function.
  p <- TF24_Parameters(hyperpar = default_hyperpar("TF24"))
  expect_true(is.function(p$hyperpar))
  p$hyperpar <- 1
  expect_error(param_hyperpar(p), "must be a trade-off function")
})

test_that("the old hyperpar names are gone, and default_hyperpar carries their job", {
  gone <- c("FF16_hyperpar", "TF24_hyperpar", "TF24f_hyperpar", "TF24t_hyperpar",
            "TF24w_hyperpar", "K93_hyperpar", "make_FF16_hyperpar",
            "make_TF24_hyperpar", "make_TF24f_hyperpar", "make_TF24t_hyperpar",
            "make_TF24w_hyperpar", "make_K93_hyperpar", "make_hyperpar", "hyperpar",
            "pulse_series_poisson", "pulse_series_regular", "pulse_series_record",
            "pulse_series_fit", "pulse_series_markov", "soil_widths_graded")
  for (nm in gone) {
    expect_false(exists(nm, envir = asNamespace("plant"), inherits = FALSE), info = nm)
  }
  # Coefficient overrides reach the frozen maker, and are recorded on the closure.
  m <- cbind(lma = 0.1, rho = 700, K_s = 0.3)
  s <- TF24_Model()
  f0 <- default_hyperpar("TF24")
  f1 <- default_hyperpar("TF24", B_hks2 = 1)
  expect_false(identical(f0(m, s, filter = FALSE), f1(m, s, filter = FALSE)))
  expect_identical(attr(f1, "coefficients")$B_hks2, 1)
  expect_identical(attr(f0, "coefficients")$B_hks2, 0.0)
  expect_true("B_kl1" %in% names(attr(default_hyperpar("FF16"), "coefficients")))
  expect_identical(attr(default_hyperpar("FF16"), "input_traits"), "narea")
  expect_error(default_hyperpar("K93", b_0 = 1), "no trade-off coefficients")
  expect_error(default_hyperpar("NOPE"), "Unknown type")
})

test_that("parameter_table groups every parameter of a model", {
  tab <- parameter_table("TF24")
  expect_setequal(tab$name[tab$group == "derived"], tf24_derived)
  expect_setequal(tab$name[tab$group == "trait"],
                  c("lma", "rho", "omega", "K_s", "vcmax_25", "jmax_25"))
  expect_setequal(tab$name[tab$group %in% c("trait", "derived", "parameter")],
                  names(TF24_Model()$pars))
  expect_setequal(tab$name[tab$group == "control"], names(Control()))
  expect_true("g_min" %in% tab$name[tab$group == "parameter"])
  expect_identical(tab$default[tab$name == "g_min"], "0")
  expect_identical(tab$description[tab$name == "eta"], "Crown-shape parameter")
  expect_identical(attr(tab, "model"), "TF24")

  # FF16's narea is read by the trade-off but is not a field of pars.
  tab <- parameter_table("FF16")
  expect_true("narea" %in% tab$name[tab$group == "trait"])
  expect_true(is.na(tab$default[tab$name == "narea"]))
  expect_setequal(tab$name[tab$group == "derived"],
                  c("k_l", "d_I", "k_s", "r_s", "r_b", "a_f3", "a_p1", "a_p2", "r_l"))

  # TF24t's thermal parameters sit beside pars.
  tab <- parameter_table("TF24t")
  expect_true(all(c("topt_offset", "k_i", "k_mat", "c_build") %in%
                    tab$name[tab$group == "strategy"]))
  expect_false("collect_all_auxiliary" %in% tab$name)
  tab <- parameter_table("TF24f")
  expect_true("k_acclim" %in% tab$name[tab$group == "strategy"])

  # With no trade-offs every field of pars is a plain parameter.
  tab <- parameter_table("TF24", hyperpar = NULL)
  expect_false(any(tab$group %in% c("trait", "derived")))
  tab <- parameter_table("K93")
  expect_false(any(tab$group %in% c("trait", "derived")))

  # A Parameters object is classified by what it carries.
  p <- default_parameters("TF24", quiet = TRUE)
  expect_identical(parameter_table(p)$group, parameter_table("TF24")$group)
  p <- base_parameters("TF24")
  expect_false(any(parameter_table(p)$group == "derived"))
  expect_error(parameter_table(1), "model name or a Parameters")
  expect_error(parameter_table("TF24", hyperpar = 1), "must be a trade-off function")
})

test_that("generate_strategy assigns TF24t's thermal fields from a trait matrix", {
  p <- default_parameters("TF24t", quiet = TRUE)
  m <- cbind(lma = 0.1, k_i = 300, topt_offset = 1.5)
  s <- generate_strategy(p, m)[[1]]
  expect_equal(s$k_i, 300)
  expect_equal(s$topt_offset, 1.5)
  expect_equal(s$pars$lma, 0.1)
  expect_false("k_i" %in% names(s$pars))
  # The trade-offs still ran on the pars columns.
  expect_false(identical(s$pars$k_l, TF24t_Model()$pars$k_l))
  # Two rows, two strategies.
  m <- cbind(k_mat = c(0.02, 0.05))
  ss <- generate_strategy(p, m)
  expect_equal(vapply(ss, function(s) s$k_mat, numeric(1)), c(0.02, 0.05))

  # Only numeric scalars beyond the infrastructure are reachable: no model has
  # another kind of top-level field today, so the filter is checked directly.
  expect_identical(
    strategy_field_names(list(pars = list(), control = list(), flag = TRUE,
                              knots = c(1, 2), gain = 0.3)),
    "gain")
  p <- default_parameters("TF24f", quiet = TRUE)
  expect_error(generate_strategy(p, trait_matrix(1, "collect_all_auxiliary")), "Unknown trait")
  expect_equal(generate_strategy(p, trait_matrix(0.3, "k_acclim"))[[1]]$k_acclim, 0.3)
  # Base TF24 has no such fields, so the same name is refused.
  expect_error(generate_strategy(default_parameters("TF24", quiet = TRUE),
                                 trait_matrix(300, "k_i")), "Unknown trait")
})

test_that("g_min is a TF24 parameter that reaches the leaf", {
  s0 <- TF24_Model()
  expect_identical(s0$pars$g_min, 0)
  s1 <- TF24_Model()
  s1$pars$g_min <- 5e-3
  theta <- rep(0.15, 3)
  a0 <- tf24_probe(s0, theta)
  a1 <- tf24_probe(s1, theta)
  expect_true(all(is.finite(a1)))
  expect_false(isTRUE(all.equal(a0, a1)))
  # Settable as a trait, and a plain parameter (no trade-off derives it).
  p <- default_parameters("TF24", quiet = TRUE)
  s <- generate_strategy(p, trait_matrix(5e-3, "g_min"))[[1]]
  expect_identical(s$pars$g_min, 5e-3)
})
