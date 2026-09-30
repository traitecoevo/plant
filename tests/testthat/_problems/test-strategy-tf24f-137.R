# Extracted from test-strategy-tf24f.R:137

# test -------------------------------------------------------------------------
s <- TF24f_Strategy()
p <- TF24f_Individual(s)
expect_equal(p$aux_size, 13)
expect_equal(length(p$internals$auxs), 13)
expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "root_mass",
    "opt_psi_stem",
    "opt_root_psi",
    "transpiration",
    "E_up_",
    "profit",
    "shadow_cost",
    "stom_cond_CO2",
    "assimilation",
    "leaf_marginal_return"
  ))
s <- TF24f_Strategy(collect_all_auxiliary=TRUE)
expect_true(s$collect_all_auxiliary)
p <- TF24f_Individual(s)
expect_equal(p$aux_size, 14)
