## Override the RcppR6-generated Leaf() constructor to catch misspelled argument
## names.  R's partial matching silently accepts any unambiguous abbreviation
## (e.g. `vcma` → `vcmax_25`), which turns typos into silent data errors.
## This wrapper inspects the raw call via sys.call() before partial matching has
## been applied and errors if any supplied name is not an exact match.

Leaf <- function(vcmax_25, c, b, psi_crit, root_c, root_b, root_psi_crit, beta2, jmax_25,
                 hk_s, a, curv_fact_elec_trans, curv_fact_colim, GSS_tol_abs,
                 vulnerability_curve_ncontrol, ci_abs_tol, ci_niter, g1_TF24,
                 beta_R_H, beta_R_V) {
  call_names <- names(as.list(sys.call())[-1])
  call_names <- call_names[nzchar(call_names)]
  valid_args <- c(
    "vcmax_25", "c", "b", "psi_crit", "root_c", "root_b", "root_psi_crit",
    "beta2", "jmax_25", "hk_s", "a", "curv_fact_elec_trans", "curv_fact_colim",
    "GSS_tol_abs", "vulnerability_curve_ncontrol", "ci_abs_tol", "ci_niter",
    "g1_TF24", "beta_R_H", "beta_R_V"
  )
  bad <- setdiff(call_names, valid_args)
  if (length(bad) > 0) {
    stop(sprintf(
      "Unknown argument(s) to Leaf(): %s\nCheck for misspelled parameter names. Valid names are:\n  %s",
      paste(bad, collapse = ", "),
      paste(valid_args, collapse = ", ")
    ), call. = FALSE)
  }
  Leaf__ctor(
    vcmax_25, c, b, psi_crit, root_c, root_b, root_psi_crit, beta2, jmax_25,
    hk_s, a, curv_fact_elec_trans, curv_fact_colim, GSS_tol_abs,
    vulnerability_curve_ncontrol, ci_abs_tol, ci_niter, g1_TF24,
    beta_R_H, beta_R_V
  )
}
