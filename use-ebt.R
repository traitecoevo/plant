library(tree)

p <- new(Parameters)
p$add_strategy(new(Strategy))

patch.c <- new(PatchCohortTop, p)
r <- pi/2
patch.c$seed_rain <- seed_rain(r)

## Add a cohort every '1'.
dt <- 1
t.next <- 1
patch.c$reset()
patch.c$add_seedling(1)
t <- patch.c$time
h <- list(patch.c$height[[1]])

while (patch.c$time < 30) {
  patch.c$step()
  if (patch.c$time > t.next) {
    patch.c$add_seedling(1)
    t.next <- t.next + dt
  }
  t <- c(t, patch.c$time)
  h <- c(h, list(patch.c$height[[1]]))
}

n <- length(h[[length(h)]])
h <- t(sapply(h, function(x) c(x, rep(NA, n-length(x)))))
matplot(t, h, type="l", col="black", lty=1)
