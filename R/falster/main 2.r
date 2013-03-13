#load data & libaries
setwd("~/Dropbox/Documents/_research/Falster - traits growth/R")
#library(smatr)
source('growthModel.r')
source('plots.r')
source('params.r')

#baseline traits
traits<-list()
traits$lma =1.11E-01
traits$rho = 608;
traits$hmat = 100; p.c_r1=0.75; p.c_r2=10;



#Mass fraction
plotMassFraction(traits, seqLog(0.2, 40, 100))

#plot production per leaf area
plotProductionComponents(traits, h=seqLog(0.2, 40, 100), env=c(1.0))
plotProductionComponents(traits, h=seqLog(0.2, 40, 100), env=c(0.2, 0.5, 1.0))

#plot marginal costs
plotMarginalCostsLfArea(traits, seqLog(0.01, 10, 100))
plotMarginalCostsHeight(traits, seqLog(0.01, 10, 100))
  
#reproduction
traits$hmat=25; p.c_r1=1.0;
plotReproductiveAllocation(traits, seqLog(1, 50, 100))
traits$hmat=10; p.c_r1=1.0;
plotReproductiveAllocation(traits, seqLog(1, 40, 100), add=TRUE)
traits$hmat=4; p.c_r1=1.0;
plotReproductiveAllocation(traits, seqLog(1, 40, 100), add=TRUE)
points(c(10,4, 25), c(0.5, 0.5, 0.5), type='p', pty='19', col="red")
arrows(13.3, 0.5, 7.5, 0.5, col="red", code=3)

arrows(10, 0.5, 13.3, 0.5, col="red")

#GROWTH RATE -----------------------------------------------
#LMA
X<-seqLog(0.005, 0.5, 100)
plotTraitvHeightGrowth(X,traits, h=0.25, env=1, option=1)
plotTraitvHeightGrowth(X,traits, h=0.25, env=seq(0.1,1,0.1), option=1)
plotTraitvHeightGrowth(X,traits, h=c(0.25,0.5, 1, 2, 4,8,16,32), env=1, option=1)

#WOOD DENSITY GROWTH RATE
X<-seq(200, 1000, 50)
plotTraitvHeightGrowth(X,traits, h=c(0.25, 0.5, 1, 5, 10, 20), env=1, option=2)

#SEED MASS GROWTH RATE
X<-seqLog(1E-6, 1E-1, 10)
plotTraitvHeightGrowth(X,traits, h=0.25, env=seq(1), option=0)

#theta
X<-seqLog(5E3, 1E4, 10)
plotTraitvHeightGrowth(X,traits, h=c(0.25, 0.5, 1, 5, 10, 20), env=1, option=3)


#plot height growth vs size
traits$hmat=50; p.c_r1=0.0;
RGR=FALSE
traits$lma = 0.01
plotSizevHeightGrowth(traits, h=seqLog(0.2, 25, 100), env=1.0, add=FALSE, RGR=RGR)
traits$lma = 0.05
plotSizevHeightGrowth(traits, h=seqLog(0.2, 25, 100), env=1.0, add=TRUE, RGR=RGR)
traits$lma = 0.1
plotSizevHeightGrowth(traits, h=seqLog(0.2, 25, 100), env=1.0, add=TRUE, RGR=RGR)


#Change in optimum with height
X<-seqLog(0.005, 0.5, 100)
plotTraitvHeightGrowth(X,traits, h=0.25, env=1, option=1)
OPT<-optForGrowth(traits, h=0.25, env=1, Range=c(min(X), max(X)))
points(OPT$maximum, OPT$objective)

h=c(0.2,0.5, 1, 2, 5, 10, 15, 20)
OPT<-NULL
for(i in 1:length(h))
  OPT[i]<-optForGrowth(traits, h=h[i], env=1, Range=c(min(X), max(X)))$maximum
plot(h, OPT, xlab="Height", ylab="optimum LMA", log='xy', type='l')


#Light compensation point-------------------------------------------------------
#LMA
plotTraitvWPLCP(seqLog(0.005, 0.5, 100), traits, h=c(0.25, 0.5, 1, 5, 10, 20), option = 1)

#also plot surival against light env

#WD
plotTraitvWPLCP(seq(200, 1000, 50), traits, h=c(0.25, 0.5, 1, 5, 10, 15), option = 2)

#WPLCP V size
traits$hmat=50
traits$lma = 0.1
plotSizevWPLCP(traits, h=seqLog(0.2, 25, 100), add=FALSE)
traits$lma = 0.25
plotSizevWPLCP(traits, h=seqLog(0.2, 25, 100), add=TRUE)
traits$lma = 0.5
plotSizevWPLCP(traits, h=seqLog(0.2, 25, 100), add=TRUE)


#figure for talk
par(mfrow=c(2,1))
X<-seqLog(0.02, 1, 100)
plotTraitvHeightGrowth(X,traits, h=0.25, env=seq(0.1,1,0.1), option=1, lty="solid")
plotTraitvWPLCP(X, traits, h=c(0.25, 0.5, 1, 5, 10, 20), option = 1, lty="solid")



