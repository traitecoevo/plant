# Milestone C increment 4 (#472 scope B / #537): exact AD growth-rate gradient.
# The templated kernel ff16_height_dt_crown_top reproduces a LIVE crown-top
# FF16_Strategy's dheight/dt (faithfulness), and reverse-mode AD gives the exact
# d(growth rate)/d(height) -- the quantity Node::growth_rate_gradient currently
# obtains by finite difference (#537 item A1) -- plus exact trait gradients.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_growth_ad <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/models/ff16_strategy.h"))
  testthat::skip_if(!any(has_hdr), "FF16 headers not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  odelia_so <- system.file("libs", "odelia.so", package = "odelia")
  plant_so <- system.file("libs", "plant.so", package = "plant")
  testthat::skip_if(!nzchar(odelia_so) || !file.exists(odelia_so),
                    "odelia shared library not found for tape linking.")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found (live FF16_Strategy needs it).")
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
      #include <XAD/XAD.hpp>
      #include <plant/models/ff16_strategy.h>
      using adt = xad::adj<double>::active_type;

      static plant::FF16_Strategy make_crown_top() {
        plant::FF16_Strategy s;
        s.control.shading_model = "crown-centre";  // -> crown-top assimilation
        s.prepare_strategy();
        return s;
      }

      // Live crown-top dheight/dt, assembled from FF16_Strategy public methods.
      // [[Rcpp::export]]
      double live_height_dt(double height, double light_E) {
        plant::FF16_Strategy s = make_crown_top();
        plant::FF16_Environment env;
        env.set_fixed_environment(light_E, 1e4);
        double area_leaf = s.area_leaf(height);
        double net = s.net_mass_production_dt(env, height, area_leaf);
        if (net <= 0.0) return 0.0;
        return s.dheight_darea_leaf(area_leaf) *
               (net * s.fraction_allocation_growth(height) *
                s.darea_leaf_dmass_live(area_leaf));
      }

      // Kernel dheight/dt (double) from the live strategy params.
      // [[Rcpp::export]]
      double kernel_height_dt(double height, double light_E) {
        plant::FF16_Strategy s = make_crown_top();
        return plant::ff16_height_dt_crown_top<double>(s.prod_pars(), height, light_E);
      }

      // Reverse-mode d(height_dt)/d(height) (A1) and d/d(lma), from live params.
      // [[Rcpp::export]]
      Rcpp::NumericVector kernel_height_dt_grad(double height, double light_E) {
        plant::FF16_Strategy s = make_crown_top();
        plant::FF16ProdPars<double> p0 = s.prod_pars();
        xad::adj<double>::tape_type tape;
        adt h = height, lma = p0.lma;
        tape.registerInput(h); tape.registerInput(lma);
        tape.newRecording();
        plant::FF16ProdPars<adt> p;
        p.lma=lma; p.rho=p0.rho; p.theta=p0.theta; p.a_b1=p0.a_b1; p.a_r1=p0.a_r1;
        p.eta_c=p0.eta_c; p.a_p1=p0.a_p1; p.a_p2=p0.a_p2;
        p.r_l=p0.r_l; p.r_s=p0.r_s; p.r_b=p0.r_b; p.r_r=p0.r_r;
        p.k_l=p0.k_l; p.k_b=p0.k_b; p.k_s=p0.k_s; p.k_r=p0.k_r;
        p.a_bio=p0.a_bio; p.a_y=p0.a_y; p.a_l1=p0.a_l1; p.a_l2=p0.a_l2;
        p.a_f1=p0.a_f1; p.a_f2=p0.a_f2; p.hmat=p0.hmat;
        adt y = plant::ff16_height_dt_crown_top<adt>(p, h, adt(light_E));
        tape.registerOutput(y); xad::derivative(y) = 1.0; tape.computeAdjoints();
        return Rcpp::NumericVector::create(xad::derivative(h), xad::derivative(lma));
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD tape symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("AD growth-rate gradient matches the live model and finite differences", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 growth-rate AD in pkgload load_all sessions.")
  compile_ff16_growth_ad()

  height <- 5; light_E <- 0.85

  # Faithfulness: kernel reproduces the live crown-top dheight/dt.
  expect_equal(kernel_height_dt(height, light_E),
               live_height_dt(height, light_E), tolerance = 1e-10)

  # A1: exact d(growth rate)/d(height) vs a fine central FD of the live model.
  g <- kernel_height_dt_grad(height, light_E)
  h <- 1e-6
  g_height_fd <- (live_height_dt(height + h, light_E) -
                  live_height_dt(height - h, light_E)) / (2 * h)
  expect_equal(g[1], g_height_fd, tolerance = 1e-6)
  expect_true(is.finite(g[2]))  # d/d(lma) available in the same sweep
})
