# Milestone C increment 3 (#472 scope B / #537): drive the AD net-production
# kernel from a LIVE, prepared FF16_Strategy. FF16_Strategy::prod_pars() gathers
# the kernel's parameter set from the strategy's actual pars + derived eta_c, so
# reverse-mode trait gradients come from the real model configuration rather than
# hand-supplied numbers. Validated vs finite differences of the same live model.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_live_ad <- function() {
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
  # FF16_Strategy pulls in the full plant/odelia/BH/Rcpp include surface, and
  # constructing one needs its compiled methods/vtable from plant.so (plus the
  # XAD tape from odelia.so).
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

      // Net production gradient w.r.t. lma and a_l1, from a LIVE prepared
      // strategy, plus the value; compared to FD of the same model in R.
      // [[Rcpp::export]]
      Rcpp::NumericVector live_netprod(double dlma, double da_l1,
                                       double height, double light_E) {
        plant::FF16_Strategy s;
        s.pars.lma += dlma; s.pars.a_l1 += da_l1;
        s.prepare_strategy();
        plant::FF16ProdPars<double> p0 = s.prod_pars();

        xad::adj<double>::tape_type tape;
        adt lma = p0.lma, a_l1 = s.pars.a_l1;
        tape.registerInput(lma); tape.registerInput(a_l1);
        tape.newRecording();
        plant::FF16ProdPars<adt> p;
        p.lma=lma; p.rho=p0.rho; p.theta=p0.theta; p.a_b1=p0.a_b1; p.a_r1=p0.a_r1;
        p.eta_c=p0.eta_c; p.a_p1=p0.a_p1; p.a_p2=p0.a_p2;
        p.r_l=p0.r_l; p.r_s=p0.r_s; p.r_b=p0.r_b; p.r_r=p0.r_r;
        p.k_l=p0.k_l; p.k_b=p0.k_b; p.k_s=p0.k_s; p.k_r=p0.k_r;
        p.a_bio=p0.a_bio; p.a_y=p0.a_y;
        adt area_leaf = plant::ff16_area_leaf<adt>(a_l1, adt(s.pars.a_l2), adt(height));
        adt y = plant::ff16_net_mass_production_crown_top<adt>(p, adt(height), area_leaf, adt(light_E));
        tape.registerOutput(y); xad::derivative(y) = 1.0; tape.computeAdjoints();
        return Rcpp::NumericVector::create(xad::value(y),
                                           xad::derivative(lma), xad::derivative(a_l1));
      }
      // [[Rcpp::export]]
      double live_netprod_value(double dlma, double da_l1, double height, double light_E) {
        plant::FF16_Strategy s;
        s.pars.lma += dlma; s.pars.a_l1 += da_l1;
        s.prepare_strategy();
        plant::FF16ProdPars<double> p = s.prod_pars();
        double area_leaf = plant::ff16_area_leaf<double>(s.pars.a_l1, s.pars.a_l2, height);
        return plant::ff16_net_mass_production_crown_top<double>(p, height, area_leaf, light_E);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD tape symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("live FF16_Strategy drives AD net-production trait gradient", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping live FF16 AD in pkgload load_all sessions.")
  compile_ff16_live_ad()

  height <- 5; light_E <- 0.8
  out <- live_netprod(0, 0, height, light_E)
  expect_equal(out[1], live_netprod_value(0, 0, height, light_E), tolerance = 1e-10)

  h <- 1e-6
  g_lma_fd <- (live_netprod_value(h, 0, height, light_E) -
               live_netprod_value(-h, 0, height, light_E)) / (2 * h)
  g_al1_fd <- (live_netprod_value(0, h, height, light_E) -
               live_netprod_value(0, -h, height, light_E)) / (2 * h)
  expect_equal(out[2], g_lma_fd, tolerance = 1e-6)
  expect_equal(out[3], g_al1_fd, tolerance = 1e-6)
})
