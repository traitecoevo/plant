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
  f <- r * p / (p.c_acc * traits$s)
  f[p < 0] <- 0
  f
}

leaf.fraction <- function(traits, h) {
  a <- LeafArea(h)
  1/(1.0 + p.a3/traits$lma +
     (traits$rho / p.theta * p.a1 * etac(p.eta) * (1.0 +p.b) *
      (1.0+p.B1) * a^p.B1 / traits$lma +
      traits$rho * p.a2 * etac(p.eta) * p.B2 * a^(p.B2-1) / traits$lma))
}

growth.rate <- function(traits, h, env) {
  r <- ReproductiveAllocation(traits$hmat, h)
  p <- net.production(traits, h, env)
  l <- leaf.fraction(traits, h)
  g <- (1 - r) * p * l
  g[p < 0] <- 0
  g
}
