## plant uses odelia's ODE solver and interpolator entirely through C++ (see
## DESCRIPTION `LinkingTo: odelia`). The XAD autodiff runtime symbols those
## templates reference live in odelia's compiled library, which odelia registers
## globally in its own `.onLoad` (odelia #28). For that to take effect before
## plant's own shared object is loaded, odelia's namespace must load first --- R
## processes a package's NAMESPACE imports before its `useDynLib`. Importing a
## symbol from odelia here guarantees that ordering, so a fresh `library(plant)`
## resolves the XAD symbols with no consumer-side linker hacks. This mirrors the
## `@importFrom Rcpp evalCpp` idiom (imported to load the namespace, not called).
##' @importFrom odelia odelia_load_dll
NULL
