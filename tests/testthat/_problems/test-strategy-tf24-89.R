# Extracted from test-strategy-tf24.R:89

# test -------------------------------------------------------------------------
expected_pars <- list(
    a_l2     = 0.306,
    S_D   = 0.25,
    a_y      = 0.7,
    a_l1     = 5.44,
    a_r1     = 0.07,
    a_b1      = 0.17,
    r_b   = 8024 / 608,
    r_l   = 39.27 / 0.1978791,
    r_r   = 217,
    r_s   = 4012/608,
    a_f3  = 3.0*3.8e-5,
    a_bio  = 0.0245,
    d_I   = 0.01,
    a_dG1   = 5.5,
    a_dG2   = 20,
    a_st1 = 0.10,
    a_st2 = 0.10,
    a_st3 = 0.8,
    a_pl0 = 0.0,
    a_pl1 = 0.0,
    a_pl2 = 0.2,
    a_pl3 = 0.05,
    a_pl4 = 3.0,
    a_p1   = 151.177775377968,
    a_p2   = 0.204716166503633,
    a_f1   = 1,
    a_f2   = 50,
    a_d0   = 0.1,
    eta    = 12,
    hmat   = 16.5958691,
    k_b    = 0.2,
    k_l   = 0.4565855,
    k_r    = 1,
    k_s   = 0.2,
    lma    = 0.1978791,
    rho    = 608,
    omega  = 3.8e-5,
    theta  = 1.0/4669,
    k_I = 0.5,
    vcmax_25 = 96,
    stem_P50 = 1.85,
    K_s = 1,
    stem_c = log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16)),
    stem_b = 1.85 /((-log(1 - 50.0 / 100.0))^(1 / (log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16))))),
    psi_crit = (1.85 /((-log(1 - 50.0 / 100.0))^(1 / (log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16))))))*log(1/0.05)^(1/(log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16)))),
    beta1 = 20000,
    TF24_beta2 = 1.5,
    TF24_cost_scale = 7.5,
    TF24_floor_lambda_o = 0,
    jmax_25 = 157.44,
    a = 0.3,
    curv_fact_elec_trans = 0.7,
    curv_fact_colim = 0.99,
    var_sapwood_volume_cost = 1,
    nmass_l = 0.013,
    nmass_s = 0.00198,
    nmass_b = 0.0034,
    nmass_r = 0.00335,
    dmass_dN = 0,
    root_depth_shape_eta = 0.2,
    root_c = 2.680147,
    root_b = 3.898245,
    root_psi_crit = 3.898245 * log(1 / 0.05)^(1 / 2.680147),
    rooting_depth_max = 1.5,
    recruitment_decay = 0,
    use_energy_balance = 0,
    d = 0.05)
expected_top <- c("pars", "control", "collect_all_auxiliary",
                    "birth_rate_x", "birth_rate_y", "is_variable_birth_rate")
s <- TF24_Strategy()
expect_inherits(s, "TF24_Strategy")
expect_identical(sort(names(s)), sort(expected_top))
expect_identical(s$control, Control())
expect_identical(s$collect_all_auxiliary, FALSE)
expect_identical(s$birth_rate_x, numeric(0))
expect_identical(s$birth_rate_y, c(1.0))
expect_identical(s$is_variable_birth_rate, FALSE)
pars_keys <- sort(names(expected_pars))
expect_identical(sort(names(s$pars)), pars_keys)
