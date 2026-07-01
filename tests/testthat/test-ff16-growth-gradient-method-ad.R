# Milestone C increment 6 (#472 scope B / #537 A1): the exact growth-rate
# gradient as a real FF16_Strategy method. growth_rate_gradient_height_ad uses
# forward-mode AD (header-only, like Leaf::dprofit_droot_collar_psi) over the
# growth kernel, compiled into plant.so -- the direct replacement for the finite
# difference in Node::growth_rate_gradient. Validated vs a fine FD of the live
# crown-top model.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_growth_method <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/models/ff16_strategy.h"))
  testthat::skip_if(!any(has_hdr), "FF16 headers not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  plant_so <- system.file("libs", "plant.so", package = "plant")
  odelia_so <- system.file("libs", "odelia.so", package = "odelia")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found.")
  testthat::skip_if(!nzchar(system.file("include", package = "BH")),
                    "BH include dir not resolvable (e.g. R CMD check sandbox).")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc)),
                         paste0("-I", shQuote(system.file("include", package = "BH")))),
    PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                     shQuote(normalizePath(odelia_so))))

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      #include <plant/models/ff16_strategy.h>
      static plant::FF16_Strategy mk() {
        plant::FF16_Strategy s;
        s.control.shading_model = "crown-centre";  // crown-top assimilation
        s.prepare_strategy();
        return s;
      }
      // [[Rcpp::export]]
      double ad_growth_grad(double height, double light_E) {
        plant::FF16_Strategy s = mk();
        plant::FF16_Environment env; env.set_fixed_environment(light_E, 1e4);
        return s.growth_rate_gradient_height_ad(height, env);
      }
      // [[Rcpp::export]]
      double live_growth_rate(double height, double light_E) {
        plant::FF16_Strategy s = mk();
        plant::FF16_Environment env; env.set_fixed_environment(light_E, 1e4);
        double area_leaf = s.area_leaf(height);
        double net = s.net_mass_production_dt(env, height, area_leaf);
        if (net <= 0.0) return 0.0;
        return s.dheight_darea_leaf(area_leaf) *
               (net * s.fraction_allocation_growth(height) *
                s.darea_leaf_dmass_live(area_leaf));
      }

      // A VARYING light profile (rising with height) -- the real patch case,
      // where d(light)/d(height) at the moving crown point is non-zero.
      static plant::FF16_Environment varying_env() {
        plant::FF16_Environment env;
        std::vector<double> st = {0,5,10,15,22, 0.30,0.50,0.68,0.85,1.0};
        env.r_init_interpolators(st);
        return env;
      }
      // [[Rcpp::export]]
      double ad_growth_grad_varying(double height) {
        plant::FF16_Strategy s = mk();
        plant::FF16_Environment env = varying_env();
        return s.growth_rate_gradient_height_ad(height, env);
      }
      // [[Rcpp::export]]
      double live_growth_rate_varying(double height) {
        plant::FF16_Strategy s = mk();
        plant::FF16_Environment env = varying_env();
        double area_leaf = s.area_leaf(height);
        double net = s.net_mass_production_dt(env, height, area_leaf);
        if (net <= 0.0) return 0.0;
        return s.dheight_darea_leaf(area_leaf) *
               (net * s.fraction_allocation_growth(height) *
                s.darea_leaf_dmass_live(area_leaf));
      }

      // The DEFAULT model: deep-crown assimilation (the GK crown integral).
      static plant::FF16_Strategy mk_deep() {
        plant::FF16_Strategy s; s.prepare_strategy(); return s;  // default shading
      }
      // [[Rcpp::export]]
      double ad_growth_grad_deep(double height) {
        plant::FF16_Strategy s = mk_deep();
        plant::FF16_Environment env = varying_env();
        return s.growth_rate_gradient_height_ad(height, env);
      }
      // [[Rcpp::export]]
      double live_growth_rate_deep(double height) {
        plant::FF16_Strategy s = mk_deep();
        plant::FF16_Environment env = varying_env();
        double area_leaf = s.area_leaf(height);
        double net = s.net_mass_production_dt(env, height, area_leaf);
        if (net <= 0.0) return 0.0;
        return s.dheight_darea_leaf(area_leaf) *
               (net * s.fraction_allocation_growth(height) *
                s.darea_leaf_dmass_live(area_leaf));
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("FF16_Strategy::growth_rate_gradient_height_ad matches a fine FD (A1)", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 growth-gradient method in pkgload load_all sessions.")
  compile_ff16_growth_method()

  for (height in c(2, 5, 9)) {
    light_E <- 0.85
    g <- ad_growth_grad(height, light_E)
    h <- 1e-6
    fd <- (live_growth_rate(height + h, light_E) -
           live_growth_rate(height - h, light_E)) / (2 * h)
    expect_equal(g, fd, tolerance = 1e-6)
  }
})

testthat::test_that("growth-rate gradient is exact in a VARYING light profile (A1, patch case)", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 growth-gradient method in pkgload load_all sessions.")
  compile_ff16_growth_method()

  # Here the crown sampling point moves through a non-flat light profile as
  # height changes, so the d(light)/d(height) term is load-bearing (a fixed-light
  # gradient would be wrong). The exact AD gradient still matches the live FD.
  for (height in c(4, 8, 12)) {
    g <- ad_growth_grad_varying(height)
    h <- 1e-6
    fd <- (live_growth_rate_varying(height + h) -
           live_growth_rate_varying(height - h)) / (2 * h)
    expect_equal(g, fd, tolerance = 1e-6)
  }
})

testthat::test_that("growth-rate gradient is exact for DEEP-CROWN (the default model)", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 growth-gradient method in pkgload load_all sessions.")
  compile_ff16_growth_method()

  # FF16's default assimilation: the Gauss-Kronrod crown integral. The AD
  # gradient differentiates through the moving GK nodes (bounds scale with
  # height), the canopy density q(z/h, z), and the light's d/dz at each node.
  for (height in c(4, 8, 12)) {
    g <- ad_growth_grad_deep(height)
    h <- 1e-6
    fd <- (live_growth_rate_deep(height + h) -
           live_growth_rate_deep(height - h)) / (2 * h)
    expect_equal(g, fd, tolerance = 1e-6)
  }
})
