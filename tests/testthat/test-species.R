strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

## Each block below puts a species into a state through a patch of one species
## (one_species_patch(), set_node_states()) and reads it back from
## patch$species[[1]], comparing it against cmp_compute_competition, a
## height-coordinate oracle. They are therefore tests of the height trapezium
## and name that coordinate: TF24 resolves Control's "auto" to birth date. The
## nodes are all introduced at time zero, which the height coordinate tolerates
## because it does not integrate over birth dates.
height_patch <- function(x) {
  one_species_patch(x, strategy_types[[x]](), coordinate = "height")
}

for (x in names(strategy_types)) {
  e <- environment_types[[x]]


  test_that("Basics", {
    patch <- height_patch(x)
    sp <- patch$species[[1]]
    plant <- Individual(x, e)(strategy_types[[x]]())
    boundary <- sp$new_node
    h0 <- boundary$height

    expect_equal(sp$size, 0)
    expect_identical(sp$height_max, h0)
    expect_identical(sp$nodes, list())
    expect_identical(sp$heights, numeric(0))
    expect_identical(sp$log_densities, numeric(0))
    expect_identical(sp$compute_competition_effect_by_nodes, numeric(0))
    expect_identical(patch$compute_competition_effect_error_by_node_for_species_i(1),
                     numeric(0))
    expect_equal(sp$ode_size, 0)
    expect_identical(sp$ode_state, numeric(0))
    expect_identical(sp$ode_rates, numeric(0))

    ## Introducing a node moves the boundary node into the node list.
    patch$introduce_new_node(1, 0)
    sp <- patch$species[[1]]
    expect_equal(sp$size, 1)

    nodes <- sp$nodes
    expect_inherits(nodes, "list")
    expect_equal(length(nodes), 1)
    expect_inherits(nodes[[1]], sprintf("Node<%s,%s>", x, e))
    expect_identical(nodes[[1]]$ode_state, boundary$ode_state)
    expect_equal(sp$heights, boundary$height)
    expect_equal(sp$log_densities, boundary$log_density)
    expect_equal(sp$compute_competition_effect_by_nodes,
                 nodes[[1]]$compute_competition(0.0))

    env <- patch$environment
    expect_identical(sp$new_node$individual$establishment_probability(env),
                     plant$establishment_probability(env))

    ## The node and the boundary node stand at the same height, so the
    ## trapezium between them has no width.
    expect_equal(sp$compute_competition(0), 0)

    set_node_states(patch, "height", 1)
    sp <- patch$species[[1]]

    h <- 0
    xx <- c(sp$new_node$height, sp$heights)
    y <- c(sp$new_node$compute_competition(h),
           sp$nodes[[1]]$compute_competition(h))

    expect_identical(sp$compute_competition(h), trapezium(xx, y))

    ## A reset empties it again.
    patch$reset()
    expect_equal(patch$species[[1]]$size, 0)
  })

  ## 1: empty species (no nodes) has no leaf area above any height:
  test_that("Empty species has no leaf area", {
    sp <- height_patch(x)$species[[1]]
    expect_equal(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(10), 0)
    expect_equal(sp$compute_competition(Inf), 0)
  })

  ## 2: Node up against boundary has no leaf area:
  test_that("species with only boundary node no leaf area", {
    patch <- height_patch(x)
    patch$introduce_new_node(1, 0)
    sp <- patch$species[[1]]
    expect_equal(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(10), 0)
    expect_equal(sp$compute_competition(Inf), 0)
  })

  cmp_compute_competition <- function(h, sp) {
    x <- c(sp$heights, sp$new_node$height)
    y <- c(sapply(sp$nodes, function(p) p$compute_competition(h)),
           sp$new_node$compute_competition(h))
    trapezium(rev(x), rev(y))
  }

  ## A species of nodes at `fractions` of four times the boundary height.
  species_at <- function(patch, fractions) {
    for (i in seq_along(fractions)) patch$introduce_new_node(1, 0)
    h_top <- patch$species[[1]]$height_max * 4
    set_node_states(patch, "height", h_top * fractions)
    patch$species[[1]]
  }

  ## 3: Single node; one round of trapezium:
  test_that("Leaf area sensible with one node", {
    patch <- height_patch(x)
    sp <- species_at(patch, 1)
    h_top <- sp$heights[[1]]

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))

    expect_identical(sp$compute_competition(h_top), 0.0)

    ## Part way up (and above bottom offspring boundary condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))

    ode_size <- sp$new_node$ode_size
    ode_state <- sp$ode_state
    p <- sp$nodes[[1]]
    expect_equal(sp$ode_size, ode_size)
    expect_equal(length(ode_state), ode_size)
    expect_identical(ode_state, p$ode_state)
  })

  test_that("Leaf area sensible with two nodes", {
    patch <- height_patch(x)
    sp <- species_at(patch, c(1, .6))
    h_top <- sp$heights[[1]]

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))
    expect_equal(sp$compute_competition(h_top), 0)
    ## Part way up (below bottom node, above boundarty condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))
    ## Within the top pair (excluding the offspring)
    expect_equal(sp$compute_competition(h_top * .8), cmp_compute_competition(h_top * .8, sp))

    ode_size <- sp$new_node$ode_size
    ode_state <- sp$ode_state
    nodes <- sp$nodes
    expect_equal(sp$ode_size, ode_size * sp$size)
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(nodes, function(p) p$ode_state)))
  })

  test_that("Leaf area sensible with three nodes", {
    patch <- height_patch(x)
    sp <- species_at(patch, c(1, .75, .6))
    h_top <- sp$heights[[1]]

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))
    expect_equal(sp$compute_competition(h_top), 0)
    ## Part way up (below bottom node, above boundarty condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))
    ## Within the top pair (excluding the offspring)
    expect_equal(sp$compute_competition(h_top * .8), cmp_compute_competition(h_top * .8, sp))

    cmp_competition_effect <- vapply(sp$nodes,
                                     function(nd) nd$compute_competition(0.0),
                                     numeric(1))
    expect_identical(sp$compute_competition_effect_by_nodes, cmp_competition_effect)

    ## The refinement signal: the local error of the competition integral,
    ## scaled by the patch's total competition.
    cmp <- local_error_integration(sp$heights, cmp_competition_effect,
                                   patch$compute_competition(0))
    expect_identical(patch$compute_competition_effect_error_by_node_for_species_i(1),
                     cmp)

    ode_size <- sp$new_node$ode_size
    ode_state <- sp$ode_state
    nodes <- sp$nodes
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(nodes, function(p) p$ode_state)))
  })
}
