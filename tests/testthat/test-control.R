
test_that("Defaults", {
  expected <- list(
    node_gradient_eps = 1e-6,
    node_gradient_direction = -1L,
    node_gradient_richardson = FALSE,
    node_gradient_richardson_depth = 4, # size_t, so not int
    node_density_in_birth_date = FALSE,

    ode_a_dydt = 0.0,
    ode_a_y = 1.0,
    ode_step_size_initial = 1e-6,
    ode_step_size_max = 5,
    ode_step_size_min = 1e-6,
    ode_tol_abs = 1e-4,
    ode_tol_rel = 1e-4,
    fixed_time_step = 0.0,

    function_integration_rule = 21, # size_t so not int
    shading_model = "", # empty = each strategy's own default
    ppa_layer_optical_depth = 0.5,
    ppa_layer_smoothing = 0.3,

    offspring_production_iterations = 1000, # size_t
    offspring_production_tol = 1e-8, # 1e-6, Had to change this...

    schedule_nsteps   = 20, # size_t
    schedule_eps      = 2e-2,
    schedule_verbose  = FALSE,
    
    # The smallest profit curvature the reverse-mode collar response will divide
    # by. Measured rather than chosen: the smallest magnitude over 1351 solved
    # interior leaf states is 0.0623, so this sits sixty times below the range the
    # model occupies.
    gradient_curvature_floor = 1e-3,

    GSS_tol_abs = 1e-1,
    # Read rather than restated: this is the leaf's choice, and it lived as four
    # separate numbers -- here, plant's Control, phylloptim's Leaf and
    # phylloptim's R control -- of which plant's was the one every stand ran on.
    vulnerability_curve_ncontrol =
      phylloptim::leaf_control()$vulnerability_curve_ncontrol,
    ci_abs_tol = 1e-3,
    ci_niter = 1000
  )

  keys <- sort(names(expected))

  ctrl <- Control()
  expect_inherits(ctrl, "Control")

  expect_identical(sort(names(ctrl)), keys)
  expect_identical(unclass(ctrl)[keys], expected[keys])
})
