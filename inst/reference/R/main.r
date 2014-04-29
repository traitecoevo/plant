#load data & libaries
## setwd("~/Dropbox/Documents/_projects/Falster - traits growth review/R/")
#library(smatr)
source('growthModel.r')
source('plots.r')
source('params.r')

#baseline traits
traits<-list()
traits$lma =1.11E-01
traits$rho = 608;
traits$hmat = 20; p.c_r1=0.75; p.c_r2=10;



#compare to matlab
h = 10;
dHdt(traits, h, env=1.0)
X<-seqLog(0.05, 1, 100)
plotTraitvHeightGrowth(X,traits, h=c(0.2,0.5, 1, 2, 5, 10, 15, 20), env=1, option=1)
WPLCP(traits, h)

#
plotMassFraction(traits, seqLog(0.2, 40, 100))
plotProductionComponents(traits, h=seqLog(0.2, 40, 100), env=c(0.2, 0.5, 1.0))

plotMarginalCostsLfArea(traits, seqLog(0.01, 10, 100))
plotMarginalCostsHeight(traits, seqLog(0.01, 10, 100))

plotReproductiveAllocation(traits, seqLog(1, 40, 100))

plotDMtdt(traits, seqLog(0.01, 20, 100), 0.8)

#LMA GROWTH RATE 
p.a4=2.86E-2*0.5;
X<-seqLog(0.005, 0.2, 100)
X<-seqLog(0.05, 1, 100)
plotTraitvHeightGrowth(X,traits, h=0.25, env=seq(0.1,1,0.1), option=1)
plotTraitvHeightGrowth(X,traits, h=c(0.2,0.5, 1, 2, 5, 10, 15, 20), env=1, option=1)
points(traits$lma, WPLCP(traits, h=0.25)$f.root)

plotTraitvHeightGrowth(X,traits, h=0.25, env=1, option=1)
OPT<-optForGrowth(traits, h=0.25, env=1, Range=c(min(X), max(X)))
points(OPT$maximum, OPT$objective)

#WPLCP
plotTraitvWPLCP(X, traits, h=c(0.25, 0.5, 1, 5, 10, 20), option = 1)

#WOOD DENSITY GROWTH RATE
X<-seq(200, 1000, 50)
plotTraitvHeightGrowth(X,traits, h=c(0.25, 0.5, 1, 5, 10, 20), env=1, option=2)
OPT<-optForGrowth(traits, h=0.25, env=1, Range=c(min(X), max(X)), option=2)
points(OPT$maximum, OPT$objective)

#WPLCP
plotTraitvWPLCP(X, traits, h=c(0.25, 0.5, 1, 5, 10, 15), option = 2)
traits2<-traits
traits2$rho<-X
plot(X, Production(traits2, h=20, env=1, print=1))

#SEED MASS
X<-seqLog(1E-6, 1E-1, 10)
plotTraitvHeightGrowth(X,traits, h=0.25, env=seq(0.1,1,0.1), option=0)

#WPLCP
plotTraitvWPLCP(X, traits, option = 0)


#plot height growth rate at smallest size
traits$lma = 0.1
plotSizevHeightGrowth(traits, h=seqLog(0.2, 25, 100), env=1.0, add=FALSE)
traits$lma = 0.01
plotSizevHeightGrowth(traits, h=seqLog(0.2, 25, 100), env=1.0, add=TRUE)
traits$lma = 0.5
plotSizevHeightGrowth(traits, h=seqLog(0.2, 25, 100), env=1.0, add=TRUE)




#Test derivative functions----------------------------------------
dHdA(10)
(Height(10.01)- Height(10.0))/0.01

dMldA(traits$lma, 10)
(LeafMass(traits$lma, 10.01)-LeafMass(traits$lma, 10))/0.01

dMsdA(traits$rho, 10)
(SapwoodMass(traits$rho, 10.01, Height(10.01))-SapwoodMass(traits$rho, 10, Height(10)))/0.01

dMhdA(traits$rho, 10)
(HeartwoodMass(traits$rho, 10.001)-HeartwoodMass(traits$rho, 10))/0.001

dMbdA(traits$rho, 10)
(BarkMass(traits$rho, 10.01,Height(10.01))-BarkMass(traits$rho, 10, Height(10.0)))/0.01

dMrdA(10)
(RootMass(10.01)-RootMass(10))/0.01

dMtdA(traits, 10)
(LiveMass(traits, 10.001)-LiveMass(traits, 10))/0.001


Production(traits, 10, 1,1)

#plot change in marginal cost against height
seqLog(0.2, 40, 100)
plot(dAdMt(traits, 10)
     
#Next steps
- compare numbers for production to matlab [error check - incorporates alls functions used thus far]
- convert m to h / ml
- growth rate
- decide on naming conventions ofr pars, functions, variables (once have seen Remko's code)

