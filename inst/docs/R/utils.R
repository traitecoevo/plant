copy_assets_figures <- function(cleanup=FALSE) {
  dest <- "vignettes"
  copy <- c("figure", "assets")
  unlink(file.path(dest, copy), recursive=TRUE)
  if (!cleanup) {
    file.copy(copy, dest, recursive=TRUE, overwrite=TRUE)
  }
}

#TODO - everything from here down can be deleted I think

empty_box <- function() {
  par(oma=c(0,0,0,0), mar=rep(0.1,4))
  plot(1,1, ann=FALSE, axes=FALSE, type='n')
  box()
}

label <- function(text, px=0.03, py=NULL, ..., adj=c(0, 1)) {
  if (is.null(py)) {
    fin <- par("fin")
    r <- fin[[1]] / fin[[2]]
    if (r > 1) { # x is longer.
      py <- 1 - px
      px <- (1 - py) / r
    } else {
      py <- 1 - px * r
    }
  }
  usr <- par("usr")
  x <- usr[1] + px*(usr[2] - usr[1])
  y <- usr[3] + py*(usr[4] - usr[3])

  ## NOTE: base 10 log:
  if (par("xlog")) {
    x <- 10^x
  }
  if (par("ylog")) {
    y <- 10^y
  }

  text(x, y, text, adj=adj, ...)
}

label_panel <- function(i, ...) {
  label(sprintf("%s)", letters[i]), ...)
}


mix_colours <- function(col1, col2, p) {
  m1 <- drop(col2rgb(col1))
  m2 <- drop(col2rgb(col2))
  m3 <- outer(m1, as.numeric(p)) + outer(m2, as.numeric(1 - p))
  ret <- rep(NA_character_, length(p))
  i <- colSums(is.na(m3)) == 0L
  ret[i] <- rgb(m3[1, i], m3[2, i], m3[3, i], maxColorValue=255)
  ret
}

## no-extrapolation version of splinefun_loglog
splinefun_loglog2 <- function(x, y, xx){
  ret <- splinefun_loglog(x, y)(xx)
  ii <- xx <= max(x, na.rm=TRUE) & xx >= min(x, na.rm=TRUE)
  ret[!ii] <- 0
  ret
}

