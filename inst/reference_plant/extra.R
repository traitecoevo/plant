## Extra functions added by RGF.

leaf.pdf <- function(z, h) {
  tmp <- (z / h)^p.eta
  2 * p.eta * (1 - tmp) * tmp / z
}

leaf.icdf <- function(x, h) {
  ((1 - sqrt(x))^(1/p.eta)) * h
}

## Assume env is a list of x and y values:
assimilation.plant <- function(h, light.env) {
  a <- LeafArea(h)
  f <- function(x)
    leaf.pdf(x, h) * Assim(a, light.env(x))
  integrate(f, 0, h)$value
}

respiration.given.height <- function(traits, h) {
  a <- LeafArea(h)
  mass.sapwood <- SapwoodMass(traits$rho, a, h)
  mass.bark <- BarkMass(traits$rho, a, h)
  mass.root <- RootMass(a)
  Respiration(a, mass.sapwood/traits$rho, mass.bark/traits$rho,
              mass.root)
}

turnover.given.height <- function(traits, h) {
  a <- LeafArea(h)
  mass.leaf <- LeafMass(traits$lma, a)
  mass.sapwood <- SapwoodMass(traits$rho, a, h)
  mass.bark <- BarkMass(traits$rho, a, h)
  mass.root <- RootMass(a)
  Turnover(traits, mass.leaf, mass.sapwood, mass.bark, mass.root)
}

net.production <- function(traits, h, env) {
  assimilation.plant(h, env) -
    respiration.given.height(traits, h) -
      turnover.given.height(traits, h)
}

fecundity.rate <- function(traits, h, env) {
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  f <- r * p / (p.c_acc + traits$s)
  f[p < 0] <- 0
  f
}

## See detail.md for the derivation of this.
##
## This verison is more explicit than the version in the C++ code for
## didactic purposes and because
leaf.fraction <- function(traits, h) {
  1 / (1 + sapwood.per.leaf.mass(traits, h)
    + bark.per.leaf.mass(traits, h)
    + root.per.leaf.mass(traits, h))
}

sapwood.per.leaf.mass <- function(traits, h){
  rho <- traits$rho
  phi <- traits$lma
  etac <- etac(p.eta)
  ml <- LeafMass(traits$lma, LeafArea(h))

  rho * etac * p.a1 * (p.B1 + 1) / (p.theta * phi) *
    (ml / phi)^p.B1
}

bark.per.leaf.mass <- function(traits, h){
  p.b * sapwood.per.leaf.mass(traits, h)
}

root.per.leaf.mass <- function(traits, h){
  phi <- traits$lma
  p.a3 / phi
}


mass.leaf.growth.rate <- function(traits, h, env) {
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  l <- leaf.fraction(traits, h)
  g <- (1 - r) * p * l
  g[p < 0] <- 0
  g
}

## Based on Daniel's code for computing dh/da.
height.growth.rate <- function(traits, h, env) {
  a <- LeafArea(h)
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  g <- dHdA(a) * dAdMt(traits, a) * p * (1-r)
  g[p < 0] <- 0
  g
}

dleaf_area_dt <- function(traits, h, env){
  a <- LeafArea(h)
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  g <- dAdMt(traits, a) * p * (1-r)
  g[p < 0] <- 0
  g
}

## sapwood area growth rate
dsapwood_area_dt <- function(traits, h, env){
  dleaf_area_dt(traits, h, env)/p.theta
}

## bark area growth rate
dbark_area_dt <- function(traits, h, env){
 p.b*dleaf_area_dt(traits, h, env)/p.theta
}

## heartwood area growth rate
dheartwood_area_dt <- function(traits, h, env){
  p.k_s*LeafArea(h)/p.theta
}

## basal area growth rate
dbasal_area_dt <- function(traits, h, env){
  dheartwood_area_dt(traits, h, env) + dsapwood_area_dt(traits, h, env) + dbark_area_dt(traits, h, env)
}

## change in basal diameter per basal area
dbasal_diam_dbasal_area <- function(basal_area){
  sqrt(pi/basal_area)
}

## basal diameter growth rate
dbasal_diam_dt <- function(traits, h, env){
 dbasal_diam_dbasal_area(basal_area(h)) * dbasal_area_dt(traits, h, env)
}

## basal area
basal_area <- function(h){
  sapwood_area(h) + bark_area(h) + heartwood_area(h)
}

## heartwood area
heartwood_area <- function(h){
  0
}

## sapwood area
sapwood_area <- function(h){
  LeafArea(h) / p.theta
}

## bark area
bark_area <- function(h){
  p.b * LeafArea(h) / p.theta
}

## Based on the above function, same algorithm as used in C++ version.
height.growth.rate.via.mass.leaf <- function(traits, h, env) {
  dmdt <- mass.leaf.growth.rate(traits, h, env)
  a <- LeafArea(h)
  p.a1 * p.B1 * (a)^(p.B1 - 1) * dmdt / traits$lma
}

p.c_d0 <- 0.01
p.c_d2 <- 5.5
p.c_d3 <- 20.0
p.c_s0 <- 0.1

mortality.rate <- function(traits, h, env) {
  p <- net.production(traits, h, env)
  a <- LeafArea(h)
  p.c_d0 +
    p.c_d2 * exp(-p.c_d3 * p / a)
}

height.at.birth <- function(traits) {
  ## TODO: These are hard-coded for now, but the upper limit could be
  ## got at directly if we had a couple of extra functions to compute
  ## height from total area.
  hmin <- 1e-16
  hmax <- 1
  f <- function(h)
    LiveMass(traits, LeafArea(h)) - traits$s
  uniroot(f, c(hmin, hmax))$root
}

leaf.mass.at.birth <- function(traits) {
  LeafMass(traits$lma, LeafArea(height.at.birth(traits)))
}

germination.probability <- function(traits, env) {
  m <- leaf.mass.at.birth(traits)
  h <- height.at.birth(traits)
  P <- net.production(traits, h, env)
  if (P > 0)
    ((P / LeafArea(h))^(-2) * p.c_s0^2 + 1)^(-1)
  else
    0
}
