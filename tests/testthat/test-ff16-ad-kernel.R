# Milestone A (#472 scope B / #537): reverse-mode AD gradient of FF16 single-plant
# net mass production (crown-top variant) w.r.t. traits, through the templated
# kernel (plant/models/ff16_production_kernel.h) that FF16_Strategy's double
# methods delegate to. Faithfulness of that kernel to the live FF16 numerics is
# covered by the existing FF16 suite (the reference-comparison test passes
# because the double methods now route through the kernel); here we check that
# the AD gradient is correct (vs finite differences).
#
# The XAD adjoint Tape<T,N> is instantiated in odelia's shared library, so this
# sourceCpp build links against it via PKG_LIBS (the odelia AD-test recipe).

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_ad_kernel <- function() {
  # Prefer the source tree (has the header in dev/load_all) then the installed
  # package (CI on the built branch); skip if neither carries the kernel header.
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
      #include <XAD/XAD.hpp>
      #include <plant/models/ff16_production_kernel.h>

      static plant::FF16ProdPars<xad::adj<double>::active_type>
      pod_ad(const std::vector<xad::adj<double>::active_type>& v) {
        plant::FF16ProdPars<xad::adj<double>::active_type> p;
        p.lma=v[0];p.rho=v[1];p.theta=v[2];p.a_b1=v[3];p.a_r1=v[4];p.eta_c=v[5];
        p.a_p1=v[6];p.a_p2=v[7];p.r_l=v[8];p.r_s=v[9];p.r_b=v[10];p.r_r=v[11];
        p.k_l=v[12];p.k_b=v[13];p.k_s=v[14];p.k_r=v[15];p.a_bio=v[16];p.a_y=v[17];
        return p;
      }
      static plant::FF16ProdPars<double> pod_d(const std::vector<double>& v) {
        plant::FF16ProdPars<double> p;
        p.lma=v[0];p.rho=v[1];p.theta=v[2];p.a_b1=v[3];p.a_r1=v[4];p.eta_c=v[5];
        p.a_p1=v[6];p.a_p2=v[7];p.r_l=v[8];p.r_s=v[9];p.r_b=v[10];p.r_r=v[11];
        p.k_l=v[12];p.k_b=v[13];p.k_s=v[14];p.k_r=v[15];p.a_bio=v[16];p.a_y=v[17];
        return p;
      }

      // [[Rcpp::export]]
      Rcpp::NumericVector ff16_netprod_grad(std::vector<double> v, double height,
                                            double area_leaf, double light_E) {
        using adt = xad::adj<double>::active_type;
        xad::adj<double>::tape_type tape;
        std::vector<adt> va(v.begin(), v.end());
        adt h=height, al=area_leaf, le=light_E;
        for (auto& x : va) tape.registerInput(x);
        tape.registerInput(h); tape.registerInput(al); tape.registerInput(le);
        tape.newRecording();
        adt y = plant::ff16_net_mass_production_crown_top(pod_ad(va), h, al, le);
        tape.registerOutput(y); xad::derivative(y)=1.0; tape.computeAdjoints();
        Rcpp::NumericVector g(v.size()+3);
        for (size_t i=0;i<v.size();++i) g[i]=xad::derivative(va[i]);
        g[v.size()]=xad::derivative(h); g[v.size()+1]=xad::derivative(al);
        g[v.size()+2]=xad::derivative(le);
        return g;
      }
      // [[Rcpp::export]]
      double ff16_netprod_value(std::vector<double> v, double height,
                                double area_leaf, double light_E) {
        return plant::ff16_net_mass_production_crown_top(pod_d(v), height, area_leaf, light_E);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD tape symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("FF16 net-production AD gradient matches finite differences", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 AD kernel in pkgload load_all sessions.")
  compile_ff16_ad_kernel()

  v <- c(0.1978791, 608, 0.0002141786, 0.17, 0.07, 0.5805, 151.177, 0.204,
         0.01979, 0.0859, 0.04, 0.2086, 0.4565, 0.2, 0.0, 1.0, 0.0245, 0.7)
  height <- 5; area_leaf <- 0.3; light_E <- 0.78

  g <- ff16_netprod_grad(v, height, area_leaf, light_E)
  inp <- c(v, height, area_leaf, light_E)
  fd <- vapply(seq_along(inp), function(i) {
    h <- 1e-6 * max(1, abs(inp[i]))
    up <- inp; dn <- inp; up[i] <- up[i] + h; dn[i] <- dn[i] - h
    f <- function(z) ff16_netprod_value(z[1:18], z[19], z[20], z[21])
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))

  expect_equal(g, fd, tolerance = 1e-6)
})
