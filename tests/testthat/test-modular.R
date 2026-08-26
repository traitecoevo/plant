
test_that("Construction", {

  ## This is a *minimal* set of tests that checks that it is possible
  ## to create the full set of different object types (Individual, Species,
  ## Patch etc) with different underlying Strategy types.  It doesn't
  ## actually try to run them though, so do that elsewhere.

  strategy_types <- get_list_of_strategy_types()
  environment_types <- get_list_of_environment_types()

  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    expect_inherits(s, paste0(x, "_Strategy"))

    p <- Individual(x, e)(s)

    expect_inherits(p, "Individual")
    expect_inherits(p, sprintf("Individual<%s,%s>", x, e))
    expect_equal(class(p$strategy), class(s))
    expect_equal(p$strategy, s)

    node <- Node(x, e)(s)

    expect_inherits(node, "Node")
    expect_inherits(node, sprintf("Node<%s,%s>", x, e))
    expect_equal(class(node$individual), class(p))

    sp <- Species(x, e)(s)

    expect_inherits(sp, "Species")
    expect_inherits(sp, sprintf("Species<%s,%s>", x, e))
    expect_equal(class(sp$new_node), class(node))

    par <- Parameters(x, e)(strategies=list(s))
    expect_inherits(par, "Parameters")
    expect_inherits(par, sprintf("Parameters<%s,%s>", x, e))
    expect_equal(par$strategies[[1]], s)
    
    env <- Environment(x)

    ctrl <- Control()
    
    pat <- Patch(x, e)(par, env, ctrl)
    expect_inherits(pat, "Patch")
    expect_inherits(pat, sprintf("Patch<%s,%s>", x, e))
    expect_equal(class(pat$species[[1]]), class(sp))

    scm <- SCM(x, e)(par, env, empty_events(), ctrl)
    expect_inherits(scm, "SCM")
    expect_inherits(scm, sprintf("SCM<%s,%s>", x, e))
    expect_equal(class(scm$patch), class(pat))

    ## Stochastic model:
    s_sp <- StochasticSpecies(x, e)(s)
    expect_inherits(s_sp, "StochasticSpecies")
    expect_inherits(s_sp, sprintf("StochasticSpecies<%s,%s>", x, e))
    expect_equal(class(s_sp$new_node), class(p))

    s_pat <- StochasticPatch(x, e)(par, env, ctrl)
    expect_inherits(s_pat, "StochasticPatch")
    expect_inherits(s_pat, sprintf("StochasticPatch<%s,%s>", x, e))
    expect_equal(class(s_pat$species[[1]]), class(s_sp))

    s_pr <- StochasticPatchRunner(x, e)(par, env, ctrl)
    expect_inherits(s_pr, "StochasticPatchRunner")
    expect_inherits(s_pr, sprintf("StochasticPatchRunner<%s,%s>", x, e))
    expect_equal(class(s_pr$patch), class(s_pat))
  }
})
