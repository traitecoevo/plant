## Towards evolutionary assembly.

library(Revolve)

## Plot the community over time; binned into
cols <- grey((32:0)/32)
col <- "#00000055"
output.dir <- "output"
time.disturbance <- c(7.5, 10.6, 20, 30, 60, 120)
xlim <- c(-2,0)
nplots <- length(time.disturbance)

op <- par(mfrow=c(nplots,2), mar=c(.5,.5,.5,.5), oma=c(5, 5, 2.5, 5))

for(i in seq_len(nplots)){

  # LOAD DATA
  res <- readRDS(sprintf("%s/assembly-%s-%s.rds", output.dir, time.disturbance[i], 1))

  # RECAST AS X,Y,T, standard form used in revolve
  res.std <- lapply(res[2:length(res)], function(x)list(x=log10(x$traits[, "lma"]), y=x$seed_rain))
  for(ii in 1:length(res.std))
      res.std[[ii]]$t <- ii

  # discretise into bins
  img <- discretise(res.std)

  # plot trait evolution
  ylim <- c(0,200)
  plot(NA, xlim=xlim, ylim=ylim, ann=FALSE, axes=FALSE)
  image(img$x, img$t, img$y, col=cols, ylab="", xlim=xr,
      xaxs="r", las=1, xlab="", add=TRUE); box()

  # X-AXIS - only print number on last plot
  labels <- FALSE
  if(i==nplots)
  labels <- 10^(-4:0)
  axis(1, at=-4:0, labels=labels)
  axis(2, at=pretty(ylim), las=1)
  text(-2, 180, labels=paste0("T=",time.disturbance[i],"yrs"), pos=4)
  if(i==1)
    mtext("Assembly", 3, line=1)

  # plot mixture at end
  x <- res[[length(res)]]
  ylim=c(-4,3)
  plot(NA, xlim=xlim, ylim=ylim, ann=FALSE, axes=FALSE)
  points(log10(x$traits[,"lma"]), log10(x$seed_rain))

# X-AXIS - only print number on last plot
  labels <- FALSE
  if(i==nplots)
    labels <- 10^(-4:0)
  axis(1, at=-4:0, labels=labels)
  axis(4, at=c(-4,-2,0,2), labels=10^c(-4,-2,0,2), las=1)
  box()
  if(i==1)
    mtext("End community", 3, line=1)
}

mtext("Time", 2, 3, outer=TRUE)
mtext("Seed rain", 4, 3, outer=TRUE)
mtext("LMA", 1, 3, outer=TRUE)
par(op)
