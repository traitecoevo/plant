for (x in c("FF16", "K93")) {
  
  context(sprintf("Tidy-patch-%s", x))

  p0 <- scm_base_parameters(x)

  if(x == "FF16")
    p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, mutant = FALSE)
  if (x == "K93")
    p1 <- expand_parameters(trait_matrix(0.059, "b_0"), p0, mutant = FALSE)

  env <- make_environment(x)
  ctrl <- scm_base_control()

  results <- run_scm_collect(p1, env, ctrl)

  test_that(sprintf("Tidy outputs %s", x), {

    expect_silent(
      results_tidy <- results %>% tidy_patch()
    )

    expect_silent(
      results_tidy_integrated <-
        results_tidy$species %>% integrate_over_size_distribution()
    )

    expect_equal(names(results), c("time", "species", "env", "net_reproduction_ratios", "patch_density", "p"))

    expect_equal(names(results_tidy), c("time", "species", "env", "net_reproduction_ratios", "p", "n_spp"))

    expect_is(results_tidy, "list")
    expect_is(results_tidy$time, "numeric")
    expect_is(results_tidy$species, "data.frame")
    expect_is(results_tidy$species, "tbl")
    expect_is(results_tidy$net_reproduction_ratios, "numeric")

    expect_is(results_tidy$env, "data.frame")
    expect_is(results_tidy$env, "tbl")

    # check values transferred correctly
    expect_equal(results_tidy$time, results$time, info = "time")
    expect_equal(results_tidy$net_reproduction_ratio, results$net_reproduction_ratios, info = "net_reproduction_ratios")
    expect_equal(results_tidy$p, results$p, info= "p")

    expect_equal(results_tidy$time, results$time, info = "time")
    expect_equal(results_tidy$n_spp, length(results$species))

    core_vars <- c("step", "time", "patch_density", "species", "cohort", "height", "mortality", "fecundity", "offspring_produced_survival_weighted", "log_density", "density")
    expect_true(all(core_vars %in% names(results_tidy$species)))

    n_spp <- length(results$species)
    for(n in seq_len(n_spp))
      for(v in c("height", "mortality", "fecundity"))
        for(coh in c(1, 5, 10))
          expect_equal(
            results_tidy$species %>% dplyr::filter(species==n, cohort==coh) %>% dplyr::pull(v),
            results$species[[n]][v, ,coh] 
            ) 
    
    core_vars_integrated <- setdiff(c("cohort"), names(results_tidy$species))

    expect_true(all(core_vars_integrated %in% names(results_tidy_integrated)))

    # TODO: test integration gives correct number?
    

    # TODO: test tidy individual

    times <- c(1, 5, 10)
    expect_silent(
      tidy_species_new <- interpolate_to_times(results_tidy$species, times)
    )
    expect_true(all(names(tidy_species_new) %in% setdiff(names(results_tidy$species), c("step"))))
    # TODO: test interpolation actually works, also failures

    if(x == "FF16")
      heights <- c(1, 5, 10)
    if(x == "K93")
      heights <- c(1, 2, 4)
    
    expect_silent(
      tidy_species_new <- interpolate_to_heights(results_tidy$species, heights)
    )
    expect_true(all(names(tidy_species_new) %in% setdiff(names(results_tidy$species), c("cohort"))))
    # TODO: test interpolation actually works, also failures
  })
}
