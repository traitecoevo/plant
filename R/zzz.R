.onLoad <- function(libname, pkgname) {
  # Wrap Leaf() to error on partial argument name matches.
  # R's partial matching silently accepts any unambiguous prefix (e.g. `vcma`
  # resolves to `vcmax_25`), turning typos into wrong-value bugs.  We derive
  # valid_args from formals() so the check stays correct automatically whenever
  # RcppR6.R is regenerated with new parameters.
  ns <- asNamespace(pkgname)
  original_Leaf <- get("Leaf", envir = ns, inherits = FALSE)
  valid_args <- names(formals(original_Leaf))

  leaf_wrapper <- function() {
    supplied <- names(as.list(sys.call())[-1])
    supplied <- supplied[nzchar(supplied)]
    bad <- setdiff(supplied, valid_args)
    if (length(bad) > 0) {
      stop(sprintf(
        "Unknown argument(s) to Leaf(): %s\nCheck for misspelled parameter names.",
        paste(bad, collapse = ", ")
      ), call. = FALSE)
    }
    env <- environment()
    do.call(original_Leaf, lapply(valid_args, function(nm) get(nm, envir = env)))
  }
  formals(leaf_wrapper) <- formals(original_Leaf)
  assign("Leaf", leaf_wrapper, envir = ns)
}
