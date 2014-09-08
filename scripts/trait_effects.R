library(tree)

base <- new(Strategy)
p0 <- base$parameters

d0 <- p0$c_d0
d1 <- p0$c_d1
rho_0 <- p0$rho
B6 <- p0$B6

mortality <- Vectorize(function(d0,d1, rho, rho_0,height, B6) tree::mortality_growth_independent(d0,d1, rho, rho_0,height, B6) , c("d0","d1","rho", "rho_0","height","B6"))

# Model 1: Effect of wood density but not of height

rho=seq(100,1000)
plot(rho, mortality(rho, height=1, d0 =d0 , d1 = 0, rho_0 =rho_0, B6=B6),
	ylim = c(1E-4, 10), type='l', log='xy', xlab= "Wood density (kg /m3",
		ylab= "Mortality")
d1 <- seq(0, 2, length.out=10)
tmp <- lapply(d1, function(x) points(rho, mortality(rho, height=1, d0 =d0 , d1 = x, rho_0 =rho_0, B6=B6), type='l', col="red"))


# Model 2: Effect of wood density interacts with height

B6 = 1
h=seq(0.1,20, length.out=400)
plot(h, mortality(rho=608, height=h, d0 =d0 , d1 = 1, rho_0 =rho_0, B6=B6),
	ylim = c(1E-4, 10), type='l',  xlab= "Plant height (m)",
		ylab= "Mortality")
tmp <- lapply(seq(0, 2, length.out=10), function(x) points(h, mortality(rho=608, height=h, d0 =d0 , d1 = 1, rho_0 =rho_0, B6=x), type='l', col="blue"))


B6 = 1
rho=seq(200,1200)
plot(rho, mortality(rho, height=1, d0 =d0 , d1 = 0, rho_0 =rho_0, B6=B6),
	ylim = c(1E-4, 10), type='l', log='xy', xlab= "Wood density (kg /m3",
		ylab= "Mortality")
d1 <- seq(0, 2, length.out=10)
tmp <- lapply(d1, function(x) points(rho, mortality(rho, height=1, d0 =d0 , d1 = x, rho_0 =rho_0, B6=B6), type='l', col="red"))
tmp <- lapply(d1, function(x) points(rho, mortality(rho, height=10, d0 =d0 , d1 = x, rho_0 =rho_0, B6=B6), type='l', col="blue"))
tmp <- lapply(d1, function(x) points(rho, mortality(rho, height=0.1, d0 =d0 , d1 = x, rho_0 =rho_0, B6=B6), type='l', col="green"))
