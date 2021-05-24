##' Thin wrapper around \code{nleqslv} and \code{dfsane}
##' @title Thin wrapper around nleqslv and dfsane
##' @param x Starting point
##' @param fn Function to solve
##' @param tol Tolerance (for \code{nleqslv} this will be both
##' absolute and relative)
##' @param maxit Maximum number of iterations.  The number of function
##' evaluations will likely exceed this.
##' @param solver The solver to use.  Either "nleqslv" or "dfsane" for now.
##' @export
nlsolve <- function(x, fn, tol=1e-6, maxit=100, solver="nleqslv") {
  solver <- match.arg(solver, c("nleqslv", "dfsane"))

  res <- switch(solver,
                nleqslv=nlsolve_nleqslv(x, fn, tol, maxit),
                dfsane=nlsolve_dfsane(x, fn, tol, maxit),
                stop("Unknown solver ", solver))

  if (!attr(res, "converged")) {
    stop(sprintf("Solver has likely failed: code=%d, msg: %s",
                 attr(res, "code"), attr(res, "message")),
         immediate.=TRUE)
  }
  
  res
}

nlsolve_nleqslv <- function(x, fn, tol=1e-6, maxit=100) {
  control <- list(xtol=tol, ftol=tol, maxit=maxit)
  sol <- nleqslv::nleqslv(x, fn, global="none", control=control)
  code <- sol$termcd
  res <- sol$x
  attributes(res) <- nlsolve_nleqslv_attr(sol)
  res
}

nlsolve_nleqslv_attr <- function(sol) {
  list(y=sol$fvec, # different to dfsane
       iter=sol$iter,
       feval=sol$nfcnt, # does not include jacobian evals
       code=sol$termcd,
       message=sol$message,
       converged=!(sol$termcd > 2 || sol$termcd < 0),
       solver="nleqslv")
}

nlsolve_dfsane <- function(x, fn, tol=1e-6, maxit=100) {
  control <- list(tol=tol, maxit=maxit, trace=FALSE)
  ## This works around `is.vector`, which returns FALSE if x has any
  ## attribute, which confuses dfsane.
  fn_vector <- function(x) as.numeric(fn(x))
  sol <- BB::dfsane(x, fn_vector, control=control, quiet=TRUE)
  res <- sol$par
  attributes(res) <- nlsolve_dfsane_attr(sol)
  res
}

nlsolve_dfsane_attr <- function(sol) {
  list(y=sol$residual,
       iter=sol$iter,
       feval=sol$feval,
       code=sol$convergence,
       message=sol$message,
       converged=sol$convergence == 0,
       solver="dfsane")
}

failed <- function(x) {
  inherits(x, "try-error")
}
