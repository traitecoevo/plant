

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

#for (x in names(strategy_types)) {
for (x in c()){ #"FF16", "FF16r")) {

  context(sprintf("Fitness-%s",x))

  test_that("Viable fitness", {
    params <- scm_base_parameters(x)
    
    # setup with infinite bounds
    expect_silent(bounds1 <- bounds_infinite("lma"))
    expect_type(bounds1, "double")
    expect_is(bounds1, "matrix")
    expect_equal(bounds1[1,], c(lower=-Inf, upper=Inf))

    # Now calaculate bounds
    expect_silent(bounds <- viable_fitness(bounds1, params))
    expect_type(bounds, "double")
    expect_is(bounds, "matrix")
    expect_equal(bounds[1, ], c(lower = 0.02533822, upper = 4.989169))
  })

  test_that("Fitness landscape", {

    expect_silent(
      lma <- trait_matrix(seq_log_range(bounds, 10), "lma")
    )

    expect_length(lma[, 1], 10)
    expect_is(lma, "matrix")
    expect_equal(lma[c(1, 10), 1], as.numeric(bounds))


    expect_silent(
      fitness <- fitness_landscape(lma, params) 
    )
    expect_length(fitness, 10)

    vals <- c(-0.000511, 7.804816, 10.096565, 10.860015, 11.085425, 10.997738, 10.557026, 9.484135, 6.973809, -0.000147)
    expect_equal(fitness, vals, tolerance = 1e-4)
  })
}
