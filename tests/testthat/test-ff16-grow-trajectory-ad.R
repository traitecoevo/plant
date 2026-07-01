# Milestone C increment 9 (#472 scope B / #537): a time-integrated trait
# gradient. Integrate a single plant's height trajectory (dheight/dt over fixed
# RK4 steps, in a fixed light) and reverse-mode-differentiate the final height
# w.r.t. traits -- a calibration gradient through the growth ODE, the bridge from
# instantaneous-rate gradients to emergent (time-integrated) outputs. Driven by a
# live FF16_Strategy's params; validated vs finite differences of the same
# trajectory.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_traj_ad <- function() {
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
  testthat::skip_if(!nzchar(odelia_so) || !file.exists(odelia_so),
                    "odelia shared library not found for tape linking.")
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

      // d(final height)/d(lma, a_p1) over a growth trajectory, one reverse sweep.
      // [[Rcpp::export]]
      Rcpp::NumericVector traj_grad(double h0, double light_E, double t_end, int n) {
        plant::FF16_Strategy s; s.prepare_strategy();
        plant::FF16ProdPars<double> p0 = s.prod_pars();
        xad::adj<double>::tape_type tape;
        adt lma = p0.lma, a_p1 = p0.a_p1;
        tape.registerInput(lma); tape.registerInput(a_p1);
        tape.newRecording();
        plant::FF16ProdPars<adt> p;
        p.lma=lma; p.rho=p0.rho; p.theta=p0.theta; p.a_b1=p0.a_b1; p.a_r1=p0.a_r1;
        p.eta_c=p0.eta_c; p.a_p1=a_p1; p.a_p2=p0.a_p2;
        p.r_l=p0.r_l; p.r_s=p0.r_s; p.r_b=p0.r_b; p.r_r=p0.r_r;
        p.k_l=p0.k_l; p.k_b=p0.k_b; p.k_s=p0.k_s; p.k_r=p0.k_r;
        p.a_bio=p0.a_bio; p.a_y=p0.a_y; p.a_l1=p0.a_l1; p.a_l2=p0.a_l2;
        p.a_f1=p0.a_f1; p.a_f2=p0.a_f2; p.hmat=p0.hmat;
        adt hT = plant::ff16_grow_height<adt>(p, adt(h0), adt(light_E), t_end, n);
        tape.registerOutput(hT); xad::derivative(hT) = 1.0; tape.computeAdjoints();
        return Rcpp::NumericVector::create(xad::value(hT),
                                           xad::derivative(lma), xad::derivative(a_p1));
      }
      // [[Rcpp::export]]
      double traj_value(double dlma, double da_p1, double h0, double light_E,
                        double t_end, int n) {
        plant::FF16_Strategy s; s.prepare_strategy();
        plant::FF16ProdPars<double> p = s.prod_pars();
        p.lma += dlma; p.a_p1 += da_p1;
        return plant::ff16_grow_height<double>(p, h0, light_E, t_end, n);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD tape symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("time-integrated growth trajectory differentiates w.r.t. traits", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 trajectory AD in pkgload load_all sessions.")
  compile_ff16_traj_ad()

  h0 <- 0.4; light_E <- 0.9; t_end <- 5; n <- 60
  out <- traj_grad(h0, light_E, t_end, n)
  expect_equal(out[1], traj_value(0, 0, h0, light_E, t_end, n), tolerance = 1e-10)
  expect_gt(out[1], h0)  # the plant grew

  e <- 1e-6
  g_lma_fd <- (traj_value(e, 0, h0, light_E, t_end, n) -
               traj_value(-e, 0, h0, light_E, t_end, n)) / (2 * e)
  g_ap1_fd <- (traj_value(0, e, h0, light_E, t_end, n) -
               traj_value(0, -e, h0, light_E, t_end, n)) / (2 * e)
  expect_equal(out[2], g_lma_fd, tolerance = 1e-6)
  expect_equal(out[3], g_ap1_fd, tolerance = 1e-6)
})
