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
  mass.leaf <- a*traits$lma
  mass.sapwood <- SapwoodMass(traits$rho, a, h)
  mass.bark <- BarkMass(traits$rho, a, h)
  mass.root <- RootMass(a)
  Respiration(mass.leaf, mass.sapwood, mass.bark,
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

fecundity_dt <- function(traits, h, env) {
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  f <- r * p / (p.a_f3 + traits$omega)
  f[p < 0] <- 0
  f
}

## See detail.md for the derivation of this.
##
## This verison is more explicit than the version in the C++ code for
## didactic purposes and because
leaf.area.deployment <- function(traits, h) {
  1 / (traits$lma + sapwood.per.leaf.area(traits, h)
    + bark.per.leaf.area(traits, h)
    + root.per.leaf.area(traits, h))
}

sapwood.per.leaf.area <- function(traits, h){
  rho <- traits$rho
  phi <- traits$lma
  etac <- etac(p.eta)
  al <- LeafArea(h)

  rho * etac * p.a_l1 * (p.a_l2 + 1) * (p.theta) * al^p.a_l2
}

bark.per.leaf.area <- function(traits, h){
  p.a_b1 * sapwood.per.leaf.area(traits, h)
}

root.per.leaf.area <- function(traits, h){
  p.a_r1
}


area.leaf.growth.dt <- function(traits, h, env) {
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  l <- leaf.area.deployment(traits, h)
  g <- (1 - r) * p * l
  g[p < 0] <- 0
  g
}

## Based on Daniel's code for computing dh/da.
height.growth.dt <- function(traits, h, env) {
  a <- LeafArea(h)
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  g <- dHdA(a) * dAdMt(traits, a) * p * (1-r)
  g[p < 0] <- 0
  g
}

competition_effect_dt <- function(traits, h, env){
  a <- LeafArea(h)
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  g <- dAdMt(traits, a) * p * (1-r)
  g[p < 0] <- 0
  g
}

## sapwood area growth rate
area_sapwood_dt <- function(traits, h, env){
  competition_effect_dt(traits, h, env)* p.theta
}

## bark area growth rate
area_bark_dt <- function(traits, h, env){
 p.a_b1 * competition_effect_dt(traits, h, env)* p.theta
}

## heartwood area growth rate
area_heartwood_dt <- function(traits, h, env){
  p.k_s*LeafArea(h)* p.theta
}

## basal area growth rate
area_stem_dt <- function(traits, h, env){
  area_heartwood_dt(traits, h, env) + area_sapwood_dt(traits, h, env) + area_bark_dt(traits, h, env)
}

## change in basal diameter_stem per basal area
ddiameter_stem_darea_stem <- function(area_stem){
  sqrt(1.0/(pi*area_stem))
}

## stem diameter growth rate
diameter_stem_dt <- function(traits, h, env){
 ddiameter_stem_darea_stem(area_stem(h)) * area_stem_dt(traits, h, env)
}

## basal area
area_stem <- function(h){
  area_sapwood(h) + area_bark(h) + area_heartwood(h)
}

## heartwood area
area_heartwood <- function(h){
  0
}

## sapwood area
area_sapwood <- function(h){
  LeafArea(h) * p.theta
}

## bark area
area_bark <- function(h){
  p.a_b1 * LeafArea(h) * p.theta
}

## Based on the above function, same algorithm as used in C++ version.
height.growth.dt.via.area.leaf <- function(traits, h, env) {
  daldt <- area.leaf.growth.dt(traits, h, env)
  a <- LeafArea(h)
  p.a_l1 * p.a_l2 * (a)^(p.a_l2 - 1) * daldt
}

p.d_I <- 0.01
p.a_dG1 <- 5.5
p.a_dG2 <- 20.0
p.a_d0 <- 0.1

mortality.dt <- function(traits, h, env) {
  p <- net.production(traits, h, env)
  a <- LeafArea(h)
  p.d_I +
    p.a_dG1 * exp(-p.a_dG2 * p / a)
}

height.at.birth <- function(traits) {
  ## TODO: These are hard-coded for now, but the upper limit could be
  ## got at directly if we had a couple of extra functions to compute
  ## height from total area.
  hmin <- 1e-16
  hmax <- 1
  f <- function(h)
    LiveMass(traits, LeafArea(h)) - traits$omega
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
    ((P / LeafArea(h))^(-2) * p.a_d0^2 + 1)^(-1)
  else
    0
}
