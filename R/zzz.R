## plant links against the odelia package for its ODE solver, interpolator and
## ODE-interface code (see DESCRIPTION LinkingTo: odelia). Some of those symbols
## live in odelia's compiled library (e.g. the XAD autodiff `Tape` runtime), so
## plant's shared object carries undefined references that must resolve to
## odelia's DLL at runtime. odelia registers its DLL privately (local), so we
## ask it to (re)load globally before any plant routine that needs those symbols
## is called.
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("odelia", quietly = TRUE) &&
      exists("odelia_load_dll", asNamespace("odelia"))) {
    odelia::odelia_load_dll(local = FALSE)
  }
}
