strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {

  # initialise birth rate per species
  s <- strategy_types[[x]]()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  
  e <- environment_types[[x]]
  plant <- Individual(x, e)(s)
  node <- Node(x, e)(s)
  
  p <- Parameters(x, e)(strategies=list(s),
                        patch_type = 'meta-population')
  env <- Environment(x)
  
  # The FF16 references below are blessed on the birth-date coordinate, named
  # here rather than inherited: Control() defaults to height, and a fixture that
  # reads its coordinate from the default silently re-points at another function
  # when the default moves.
  #
  # It has to be named in two places. A Patch imposes its own control on its
  # species, but a standalone Node reads the coordinate from its own strategy's
  # control, so the comparison node below is only the same function as the
  # patch's node if the strategy carries the coordinate too.
  ctrl <- Control(node_density_in_birth_date = TRUE)
  p$strategies[[1]]$control <- ctrl
  patch <- Patch(x, e)(p, env, ctrl)
  cmp <- Node(x, e)(p$strategies[[1]])

  test_that(sprintf("Basics %s", x), {
    ## Default birth rate (see #333): a fresh strategy should be a constant
    ## rate of 1 with the variable-rate (spline) path switched off.
    s_default <- strategy_types[[x]]()
    expect_equal(s_default$birth_rate_y, 1)
    expect_false(s_default$is_variable_birth_rate)

    expect_equal(patch$size, 1)
    expect_identical(patch$height_max, cmp$height)
    expect_equal(patch$parameters, p)
    expect_equal(patch$get_area, 1.0)
    expect_inherits(patch$environment, c(paste0(x, "_Environment"), "R6"))
    expect_identical(patch$environment$time, 0.0)

    expect_equal(length(patch$species), 1)
    expect_inherits(patch$species[[1]], sprintf("Species<%s,%s>",x,e))
    
    # with no nodes, we only expect env vars
    env_size <- env$ode_size
    #fails here
    env_state <- patch$ode_state
    env_rates <- patch$ode_rates
    expect_equal(patch$ode_size, env_size)
    
    # either 0 or numeric(0)
    if(x %in% c("FF16", "K93")) {
      expect_equal(patch$ode_state, numeric(0))
      expect_equal(patch$ode_rates, numeric(0))
    }
    if(x %in% c("TF24")) {
      length_odes <- env$get_soil_number_of_depths()
      soil_moist_inits <- c(rep(env$soil_moist_sat, length_odes)/2, rep(0,5))
      expect_equal(patch$ode_state, soil_moist_inits)
      expect_equal(patch$ode_rates, c(3.312786717, rep(0,4), 1.000000000, 0.996093750, 0.002257735, 0.000000000, 0.000000000), tolerance = 1e-4)
    }
    
    expect_identical(patch$ode_state, env_state)
    expect_identical(patch$ode_rates, env_rates)
    
    ## Empty environment:
    patch$compute_environment()
    expect_identical(patch$compute_competition(0), 0)
    
    expect_error(patch$introduce_new_node(0, 0), "Invalid value")
    expect_error(patch$introduce_new_node(2, 0), "out of bounds")
    
    # introduce a node and expect different results
    node_size <- Node(x, e)(s)$ode_size
    ode_size = node_size + env_size
    patch$introduce_new_node(1, 0)
    expect_equal(patch$node_ode_size, node_size)
    expect_equal(patch$ode_size, ode_size)
    if (x == "FF16") {
      expect_equal(patch$node_ode_size, 7)
      expect_equal(patch$ode_size, 7)
    }
    ## Then pull this out:
    cmp$compute_initial_conditions(patch$environment, patch$pr_survival(0.0), 
                                   patch$species[[1]]$extrinsic_drivers$evaluate("birth_rate", 0))
     
    ode_state <- c(cmp$ode_state, env_state)
    ode_rates <- c(cmp$ode_rates, env_rates)
    expect_identical(patch$ode_state, ode_state)
    # The patch's boundary node is evaluated by the field build, in a field that
    # excludes the boundary interval; this comparison node is seeded afterwards in
    # the completed field. So the two agree except through
    # pr_patch_survival_at_birth, which divides the fecundity rate and is fixed at
    # seeding -- about 1% of a quantity of order 1e-21.
    expect_equal(patch$ode_rates, ode_rates, tolerance = 1e-2)
    if (x == "FF16") {
      # Blessed on the birth-date coordinate, which ctrl names above. The two
      # coordinates are different functions rather than two discretisations of
      # one, and the last entry of each vector is where they part company -- so
      # the two literals below are asserted alongside the relations that produce
      # them, which is what says WHICH function this is a reference for.
      expect_true(ctrl$node_density_in_birth_date)
      expect_equal(ode_state, c(0.3441947, 0.009159, 0, 0, 0, 0, -0.009159),
                   tolerance = 1e-4)
      expect_equal(ode_rates,
                   c(0.3341652, 0.01000000, 0, 5.1781e-09, 9.60270e-07, 0, -0.01),
                   tolerance = 1e-4)

      # The density rate is minus the mortality rate and nothing else. On the
      # height coordinate it carries a compression term as well, which is what
      # made this entry -0.78726 and cost a displaced-height solve to compute.
      expect_equal(ode_rates[[7]], -ode_rates[[2]])

      # And the boundary condition is the seed arrival times the establishment
      # probability, with no division by the growth rate. At a birth rate of one
      # that makes the log density exactly minus the cumulative mortality, since
      # the mortality a node is seeded with is -log(pr_estab). On the height
      # coordinate the division put this entry -log(g) higher, at 1.08695.
      expect_equal(patch$species[[1]]$extrinsic_drivers$evaluate("birth_rate", 0), 1)
      expect_equal(ode_state[[7]], -ode_state[[2]])
    }
    y <- patch$ode_state
    patch$set_ode_state(y, 0)
    expect_identical(patch$ode_state, y)
    
    ## NOTE: These should be identical, but are merely equal...
    expect_equal(patch$derivs(y, 0), ode_rates)

    patch$reset()
    expect_equal(patch$ode_size, env_size)
    expect_identical(patch$environment$time, 0.0)
  })

  test_that("change patch size", {
  
    ctrl <- Control()
    e <- environment_types[[x]]
    env <- Environment(x)

    # TF24 is much costlier per step, so use a shorter horizon -- but with a low
    # height-at-maturity (hmat = 5) so plants actually reproduce within it. With
    # the default hmat over so short a patch the TF24 reproductive ratios
    # underflow to ~1e-15 and the checks below would merely compare zeros;
    # hmat = 5 lifts them to O(10) so this is a genuine demographic check.
    # FF16/K93 reproduce fine at their default hmat over the longer horizon.
    if (x == "TF24") {
      max_patch_lifetime <- 4
      mk_strategy <- function()
        add_strategies(scm_base_parameters("TF24"),
                       trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = 1)$strategies[[1]]
    } else {
      max_patch_lifetime <- 30
      mk_strategy <- function() strategy_types[[x]]()
    }

    mk_pars <- function(area) Parameters(x, e)(strategies = list(mk_strategy()),
                                               patch_area = area,
                                               max_patch_lifetime = max_patch_lifetime)

    p2  <- mk_pars(2)
    p10 <- mk_pars(10)
    expect_equal(p2$patch_area, 2)
    expect_equal(Patch(x, e)(p2, env, ctrl)$get_area, 2)
    expect_equal(p10$patch_area, 10)
    expect_equal(Patch(x, e)(p10, env, ctrl)$get_area, 10)

    # Same birth rate: the larger patch has less competition, so a higher
    # net reproductive ratio.
    p2$strategies[[1]]$birth_rate_y <- 1
    p10$strategies[[1]]$birth_rate_y <- 1
    expect_warning(scm2   <- run_scm(p2))
    expect_warning(scm10a <- run_scm(p10))
    expect_gt(scm10a$net_reproduction_ratios, scm2$net_reproduction_ratios)

    # Birth rate scaled x5 with the x5 patch area gives the same seed input per
    # unit area, hence the same reproductive ratio as the small patch. (p2 is
    # unchanged from above, so it is not re-run.)
    p10b <- p10
    p10b$strategies[[1]]$birth_rate_y <- 5
    expect_warning(scm10b <- run_scm(p10b))
    expect_equal(scm2$net_reproduction_ratios,
                 scm10b$net_reproduction_ratios, tolerance = 1e-4)
    expect_gt(scm10a$net_reproduction_ratios, scm10b$net_reproduction_ratios)
  })
  
  test_that("Weibull Disturbance as default", {
    expect_equal(patch$time, 0.0)
    expect_equal(patch$pr_survival(patch$time), 1.0)
    
    expect_equal(patch$disturbance_mean_interval(), 30)
    disturbance <- Weibull_Disturbance_Regime(105.32)
    
    patch$set_time(10)
    expect_equal(patch$time, 10)
    
    expect_identical(patch$pr_survival(10), 
                     disturbance$pr_survival(10))

    # This is how we'd like it but Rcpp wouldn't handle a disturbance pointer
    #expect_inherits(patch$disturbance_regime, "Disturbance")
  })
  
  test_that("No Disturbance for fixed-time patches", {
    p$patch_type <- "fixed"
    env <- Environment(x)
    ctrl <- Control()

    patch <- Patch(x, e)(p, env, ctrl)
    
    expect_identical(patch$time, 0.0)
    expect_identical(patch$pr_survival(patch$time), 1.0)
    
    expect_true(is.na(patch$disturbance_mean_interval()))
    disturbance <- No_Disturbance()
    
    patch$set_time(10)
    expect_equal(patch$time, 10)
    
    expect_identical(patch$pr_survival(10), 1)
    
    expect_identical(patch$pr_survival(10), 
                     disturbance$pr_survival(10))
    
    # This is how we'd like it but Rcpp wouldn't handle a disturbance pointer
    #expect_inherits(patch$disturbance_regime, "Disturbance")
  })

  test_that(sprintf("Patch aux carries one slot per resource %s", x), {
    # The environment's slots follow the whole species range, so a model with no
    # resources keeps the width it had.
    n_resources <- if (x == "TF24") Environment(x)$get_soil_number_of_depths() else 0
    per_node <- Individual(x, e)(s)$aux_size

    pa <- Patch(x, e)(p, Environment(x), Control())
    expect_equal(length(pa$ode_aux), n_resources)
    for (n in 1:2) {
      pa$introduce_new_node(1, n * 1.0)
      expect_equal(length(pa$ode_aux), n * per_node + n_resources)
    }
  })
}

test_that("TF24 patch aux reports the per-layer uptake", {
  s <- get_list_of_strategy_types()$TF24()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- Parameters("TF24", "TF24_Env")(strategies = list(s),
                                      patch_type = "meta-population")
  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  patch$introduce_new_node(1, 0)

  # A node at its birth height takes up nothing measurable, so grow it and give
  # it a density before reading the uptake.
  names_ode <- Node("TF24", "TF24_Env")(s)$ode_names
  y <- patch$ode_state
  y[[which(names_ode == "height")]] <- 10
  y[[which(names_ode == "log_density")]] <- log(5)

  n_layers <- Environment("TF24")$get_soil_number_of_depths()
  # A dry layer under a wet one, because that is the arrangement signed uptake
  # needs and it is not the one an arbitrary state supplies. Water moves on the
  # difference between a layer's potential and the collar's, so where a deep
  # layer is drier than the collar the flux runs the other way and the plant
  # releases into it. On a uniform column every layer is drawn from and the
  # signedness below would be asserting a phenomenon the state does not contain.
  #
  # The profile is measured rather than chosen: at these contents the deepest two
  # layers come back at -0.158 and -0.169 against +1.7 at the top. It is written
  # for the five-layer column and refuses a wider one rather than padding with
  # NAs, because a silently truncated profile would be a uniform column again.
  moisture <- c(0.35, 0.32, 0.25, 0.12, 0.06)
  expect_lte(n_layers, length(moisture))
  at <- length(y) - Environment("TF24")$ode_size
  y[at + seq_len(n_layers)] <- moisture[seq_len(n_layers)]
  patch$set_ode_state(y, 1.0)
  # The auxiliaries are written by a rate evaluation and not by loading a state,
  # so a read taken before one returns the previous evaluation's values -- zeros,
  # on a patch whose rates have never been computed at this state. This was
  # passing by accident: the height coordinate's density rate solves the
  # physiology again at a displaced height, which populated them on the way.
  invisible(patch$ode_rates)
  uptake <- tail(patch$ode_aux, n_layers)

  # The soil's cumulative total-uptake rate is the last environment rate and is
  # accumulated on its own, so it checks the published values and not the width.
  expect_equal(sum(uptake), tail(patch$ode_rates, 1))
  expect_true(all(is.finite(uptake)))
  expect_gt(uptake[[1]], 0)
  # Uptake is signed, and the state above is what makes that observable rather
  # than incidental: the deepest layer is drier than the collar, so the plant
  # releases into it instead of drawing from it.
  expect_lt(uptake[[n_layers]], 0)
})

# A rate evaluation is a function of the state and the time it is given, and of
# nothing the patch is carrying from the last one.
#
# ⚠️ THIS IS THE ONE PROPERTY NOTHING ELSE CHECKS, and it has already been broken
# three times. `Internals` holds what an evaluation is given beside what it
# produces, so a produced value read before it is written this pass reads the
# PREVIOUS evaluation's -- and on the tape it is worse than stale, because an
# active read before it is written is an unregistered input to the recorded
# function whose row is then dropped in silence. The rates and consumption rates
# are filled with NA at construction and so say so on the first pass; the
# auxiliaries are filled with zero, which is a plausible number, and they have the
# most readers.
#
# Revisiting a state is what makes it observable: an evaluation that reads
# something it did not write reads what the intervening state left, so the two
# visits disagree. Both published channels, because they are written in different
# places -- the rates by compute_rates, the auxiliaries by compute_rates AND by
# the state load.
for (x in names(strategy_types)) {
  test_that(sprintf("a rate evaluation carries nothing between calls (%s)", x), {
    s <- strategy_types[[x]]()
    s$birth_rate_y <- 1
    s$is_variable_birth_rate <- FALSE
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies = list(s), patch_type = "meta-population")
    patch <- Patch(x, e)(p, Environment(x), Control())
    patch$introduce_new_node(1, 0)

    y <- patch$ode_state
    # Displaced by a relative step, floored so a zero state still moves.
    away <- y + 0.05 * pmax(abs(y), 1e-3)
    time <- 0.1

    first <- patch$derivs(y, time)
    first_aux <- patch$ode_aux
    invisible(patch$derivs(away, time))
    again <- patch$derivs(y, time)
    again_aux <- patch$ode_aux

    # Bit-exact, not close: the two visits run the same arithmetic on the same
    # numbers, so anything else is state that survived the first one.
    expect_identical(again, first)
    expect_identical(again_aux, first_aux)
    # The displaced state must actually reach the rates, or the two visits agree
    # for the reason that nothing happened.
    expect_false(identical(patch$derivs(away, time), first))
  })
}

test_that("TF24 patch aux goes back the way it came", {
  s <- get_list_of_strategy_types()$TF24()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- Parameters("TF24", "TF24_Env")(strategies = list(s),
                                      patch_type = "meta-population")
  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  patch$introduce_new_node(1, 0)
  patch$set_ode_state(patch$ode_state, 1.0)

  aux <- patch$ode_aux
  expect_gt(length(aux), 0)

  # Restoring what a rate evaluation published is what the reverse pass does with a
  # stage's aux, so the round trip has to be exact rather than close.
  patch$set_ode_aux(rev(aux))
  expect_identical(patch$ode_aux, rev(aux))
  patch$set_ode_aux(aux)
  expect_identical(patch$ode_aux, aux)

  # The width is what the solver checks, so a wrong one stops here.
  expect_error(patch$set_ode_aux(aux[-1]), "Incorrect length input")
})

