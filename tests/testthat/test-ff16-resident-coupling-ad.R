# Milestone B increment 2 (#472 scope B / #537): the full resident self-shading
# coupling. A trait flows through the templated FF16 area_leaf into BOTH (a) the
# plant's own mass cascade + assimilation and (b) the competition it exerts,
# which sets the resident light spline's knot values y_i = exp(-density * k_I *
# area_leaf * Qfrac_i) (Beer's law; Qfrac_i = fraction of leaf above knot i, a
# frozen geometric weight). Net production through the real deep-crown integral
# over that AD light spline differentiates correctly w.r.t. a_l1 (everywhere) and
# k_I (through the light only), in one reverse sweep -- composing the FF16 kernel,
# odelia's differentiable interpolator, and a frozen quadrature.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_resident_ad <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/models/ff16_production_kernel.h"))
  testthat::skip_if(!any(has_hdr),
                    "FF16 AD kernel header not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  odelia_so <- system.file("libs", "odelia.so", package = "odelia")
  testthat::skip_if(!nzchar(odelia_so) || !file.exists(odelia_so),
                    "odelia shared library not found for tape linking.")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc))),
    PKG_LIBS = shQuote(normalizePath(odelia_so)))

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      #include <vector>
      #include <cmath>
      #include <XAD/XAD.hpp>
      #include <odelia/interpolator.hpp>
      #include <plant/models/ff16_production_kernel.h>
      using adt = xad::adj<double>::active_type;

      template <typename S>
      S resident_net(S a_l1, S k_I, const plant::FF16ProdPars<S>& p,
                     double a_l2, double height, double density,
                     std::vector<double> xk, std::vector<double> qfrac_knot,
                     std::vector<double> z, std::vector<double> wq) {
        S area_leaf = plant::ff16_area_leaf<S>(a_l1, S(a_l2), S(height));
        // Resident light spline: Beer\'s law over its own competition.
        std::vector<S> yk(xk.size());
        for (size_t i = 0; i < xk.size(); ++i)
          yk[i] = exp(-density * k_I * area_leaf * qfrac_knot[i]);
        odelia::interpolator::basic_interpolator<S> light;
        light.init(xk, yk);
        S assim = plant::ff16_assimilation_deep_crown_replay<S>(
            p.a_p1, p.a_p2, area_leaf, z, wq,
            [&](double zz){ return light.eval(zz); });
        return plant::ff16_net_from_components<S>(p, S(height), area_leaf, assim);
      }

      template <typename S>
      plant::FF16ProdPars<S> mkpod(std::vector<double> v) {
        plant::FF16ProdPars<S> p;
        p.lma=v[0];p.rho=v[1];p.theta=v[2];p.a_b1=v[3];p.a_r1=v[4];p.eta_c=v[5];
        p.a_p1=v[6];p.a_p2=v[7];p.r_l=v[8];p.r_s=v[9];p.r_b=v[10];p.r_r=v[11];
        p.k_l=v[12];p.k_b=v[13];p.k_s=v[14];p.k_r=v[15];p.a_bio=v[16];p.a_y=v[17];
        return p;
      }

      // [[Rcpp::export]]
      Rcpp::NumericVector resident_grad(double a_l1, double k_I, std::vector<double> v,
          double a_l2, double height, double density, std::vector<double> xk,
          std::vector<double> qfrac_knot, std::vector<double> z, std::vector<double> wq) {
        xad::adj<double>::tape_type tape;
        adt al1 = a_l1, ki = k_I;
        tape.registerInput(al1); tape.registerInput(ki);
        tape.newRecording();
        adt y = resident_net<adt>(al1, ki, mkpod<adt>(v), a_l2, height, density,
                                  xk, qfrac_knot, z, wq);
        tape.registerOutput(y); xad::derivative(y) = 1.0; tape.computeAdjoints();
        return Rcpp::NumericVector::create(xad::derivative(al1), xad::derivative(ki));
      }
      // [[Rcpp::export]]
      double resident_value(double a_l1, double k_I, std::vector<double> v,
          double a_l2, double height, double density, std::vector<double> xk,
          std::vector<double> qfrac_knot, std::vector<double> z, std::vector<double> wq) {
        return resident_net<double>(a_l1, k_I, mkpod<double>(v), a_l2, height, density,
                                    xk, qfrac_knot, z, wq);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD tape symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("FF16 resident self-shading coupling differentiates end-to-end", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 resident-coupling AD in pkgload load_all sessions.")
  compile_ff16_resident_ad()

  H <- 8; density <- 1.5
  xk <- seq(0, H, length.out = 11)
  qfrac_knot <- pmax(0, 1 - xk / H)        # leaf-above fraction at each spline knot
  nq <- 41; z <- seq(0, H, length.out = nq); dz <- H / (nq - 1)
  simp <- ifelse(seq_len(nq) %in% c(1, nq), 1, ifelse(seq_len(nq) %% 2 == 0, 4, 2))
  q <- pmax(0, 1 - z / H)
  wq <- simp * dz / 3 * q

  a_l1 <- 0.306; a_l2 <- 0.75; k_I <- 0.5
  v <- c(0.1978791, 608, 0.0002141786, 0.17, 0.07, 0.5805, 151.177, 0.204,
         0.01979, 0.0859, 0.04, 0.2086, 0.4565, 0.2, 0.0, 1.0, 0.0245, 0.7)

  g <- resident_grad(a_l1, k_I, v, a_l2, H, density, xk, qfrac_knot, z, wq)
  f <- function(al1, ki) resident_value(al1, ki, v, a_l2, H, density, xk,
                                        qfrac_knot, z, wq)
  h <- 1e-6
  g_al1_fd <- (f(a_l1 + h, k_I) - f(a_l1 - h, k_I)) / (2 * h)
  g_kI_fd  <- (f(a_l1, k_I + h) - f(a_l1, k_I - h)) / (2 * h)

  expect_equal(g[1], g_al1_fd, tolerance = 1e-6)  # a_l1: area_leaf -> everything
  expect_equal(g[2], g_kI_fd, tolerance = 1e-6)   # k_I: through the light spline only
})
