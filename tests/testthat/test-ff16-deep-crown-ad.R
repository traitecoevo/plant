# Milestone B increment 1 (#472 scope B / #537): the deep-crown assimilation
# integral as a frozen-replay weighted sum over an AD-capable resident light
# spline. Validates that net production through the REAL default assimilation
# variant differentiates correctly w.r.t. (a) a physiology trait (a_p1) and
# (b) a stand trait `theta` that enters through the light spline's knot values
# (the resident self-shading coupling) -- composing the FF16 kernel, odelia's
# differentiable interpolator, and a frozen quadrature, all in one reverse sweep.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_deep_crown_ad <- function() {
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

      // theta enters the resident light spline knot values: y_i = exp(-theta*s_i).
      template <typename S>
      S net_with_deep_crown(S theta, S a_p1, double a_p2, double area_leaf,
                            std::vector<double> xk, std::vector<double> shade,
                            std::vector<double> z, std::vector<double> wq,
                            double resp, double turn, double a_bio, double a_y) {
        std::vector<S> yk(xk.size());
        for (size_t i = 0; i < xk.size(); ++i) yk[i] = exp(-theta * shade[i]);
        odelia::interpolator::basic_interpolator<S> light;
        light.init(xk, yk);
        S assim = plant::ff16_assimilation_deep_crown_replay<S>(
            a_p1, S(a_p2), S(area_leaf), z, wq,
            [&](double zz){ return light.eval(zz); });
        return plant::ff16_net_production_A<S>(S(a_bio), S(a_y), assim, S(resp), S(turn));
      }

      // [[Rcpp::export]]
      Rcpp::NumericVector deep_crown_grad(double theta, double a_p1, double a_p2,
          double area_leaf, std::vector<double> xk, std::vector<double> shade,
          std::vector<double> z, std::vector<double> wq, double resp, double turn,
          double a_bio, double a_y) {
        xad::adj<double>::tape_type tape;
        adt th = theta, ap1 = a_p1;
        tape.registerInput(th); tape.registerInput(ap1);
        tape.newRecording();
        adt y = net_with_deep_crown<adt>(th, ap1, a_p2, area_leaf, xk, shade, z, wq,
                                         resp, turn, a_bio, a_y);
        tape.registerOutput(y); xad::derivative(y) = 1.0; tape.computeAdjoints();
        return Rcpp::NumericVector::create(xad::derivative(th), xad::derivative(ap1));
      }
      // [[Rcpp::export]]
      double deep_crown_value(double theta, double a_p1, double a_p2,
          double area_leaf, std::vector<double> xk, std::vector<double> shade,
          std::vector<double> z, std::vector<double> wq, double resp, double turn,
          double a_bio, double a_y) {
        return net_with_deep_crown<double>(theta, a_p1, a_p2, area_leaf, xk, shade,
                                           z, wq, resp, turn, a_bio, a_y);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD tape symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("deep-crown frozen-replay assimilation differentiates (trait + light coupling)", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 deep-crown AD in pkgload load_all sessions.")
  compile_ff16_deep_crown_ad()

  H <- 6
  xk <- seq(0, H, length.out = 9)
  shade <- 0.15 * (H - xk)               # taller canopy above -> more shade low down
  # Frozen quadrature schedule (stand-in for the discovered crown integral):
  # composite Simpson nodes/weights on [0,H], folding a leaf-density factor q(z).
  nq <- 41
  z <- seq(0, H, length.out = nq); dz <- H / (nq - 1)
  simp <- ifelse(seq_len(nq) %in% c(1, nq), 1, ifelse(seq_len(nq) %% 2 == 0, 4, 2))
  q <- pmax(0, 1 - z / H)                # leaf-area density weight, frozen
  wq <- simp * dz / 3 * q

  theta <- 0.8; a_p1 <- 151.177; a_p2 <- 0.204; area_leaf <- 0.3
  resp <- 0.05; turn <- 0.08; a_bio <- 0.0245; a_y <- 0.7

  g <- deep_crown_grad(theta, a_p1, a_p2, area_leaf, xk, shade, z, wq,
                       resp, turn, a_bio, a_y)
  f <- function(th, ap1) deep_crown_value(th, ap1, a_p2, area_leaf, xk, shade, z, wq,
                                          resp, turn, a_bio, a_y)
  h <- 1e-6
  g_theta_fd <- (f(theta + h, a_p1) - f(theta - h, a_p1)) / (2 * h)
  g_ap1_fd   <- (f(theta, a_p1 + h) - f(theta, a_p1 - h)) / (2 * h)

  expect_equal(g[1], g_theta_fd, tolerance = 1e-6)  # through the light spline
  expect_equal(g[2], g_ap1_fd, tolerance = 1e-6)    # through assimilation_leaf
})
