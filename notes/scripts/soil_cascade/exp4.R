source("soil_cascade.R"); source("exp2_ctl.R"); options(digits=6)

run <- function(h_start, mode, fix_nan=FALSE, t_end=0.5, h_max=5) {
  y <- c(P$theta_sat, 0.05,0.05,0.05,0.05); t <- 0; h <- h_start
  ctl <- CTL; ctl$h_max <- h_max
  a_psi <- 1.78e-3*1e6
  psi_of <- function(th){ tt <- pmax(th, P$theta_res); a_psi*(tt/P$theta_sat)^(-P$n_psi)/1e6 }
  upt <- function(th){ psi <- psi_of(th); bad <- !is.finite(psi)|psi<=0|th>P$theta_sat
    u <- rep(0.05,length(th))
    if(any(bad)){ if(mode=="throw") stop("find_root_psi failed"); if(mode=="nan") u[bad] <- NaN }
    u }
  ns<-0; nr<-0; tmax<-max(y); status<-"completed"
  while (t < t_end && ns < 2e5) {
    hh <- if (t+h > t_end) t_end-t else h
    st <- tryCatch(stage_with_uptake(y, hh, "baseline", upt),
                   error=function(e) structure(list(msg=conditionMessage(e)), class="thrown"))
    if (inherits(st,"thrown")) { status <- paste("DIED:", st$msg); break }
    a <- if (fix_nan && (any(!is.finite(st$yerr))||any(!is.finite(st$y))))
           list(h=max(hh*0.2, ctl$h_min), shrank=TRUE) else adjust(hh,5,st$y,st$yerr,st$dydt_out,ctl)
    if (a$shrank) { nr<-nr+1; if (a$h>=hh) { status<-"DIED: cannot achieve accuracy"; break }; h<-a$h }
    else { t<-t+hh; y<-st$y; ns<-ns+1; tmax<-max(tmax,max(y)); h<-a$h
           if (any(!is.finite(y))) { status<-"DIED: non-finite state ACCEPTED"; break } }
  }
  list(ns=ns,nr=nr,tmax=tmax,status=status,t=t)
}

cat("=== Q12. Entering a leg with a LARGE inherited step (step_size_last leak) ===\n")
for (h0 in c(1e-6, 1e-3, 0.08, 1.0, 5.0)) {
  for (mode in c("throw","nan")) {
    r <- run(h0, mode)
    cat(sprintf("  h_start=%6.4g mode=%-5s : t=%.4f steps=%5d rej=%4d max theta=%11.5g  %s\n",
                h0, mode, r$t, r$ns, r$nr, r$tmax, r$status))
  }
}
cat("\n=== Q13. Same, with the NaN controller bug fixed ===\n")
for (h0 in c(0.08, 1.0, 5.0)) for (mode in c("throw","nan")) {
  r <- run(h0, mode, fix_nan=TRUE)
  cat(sprintf("  h_start=%6.4g mode=%-5s : t=%.4f steps=%5d rej=%4d max theta=%11.5g  %s\n",
              h0, mode, r$t, r$ns, r$nr, r$tmax, r$status))
}
