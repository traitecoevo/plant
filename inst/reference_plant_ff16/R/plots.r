#plot mass growth rate
plotDMtdt<-function(traits, h, env){
  plot(h, dMtdt(traits, h, env), type='l', log="x",col='black', xaxs="i", yaxs="i",  ylab="mass growth rate (kg /yr)")
  }


#plot plotMarginalCosts
plotMarginalCostsHeight<-function(traits, h){
dml<-dMldH(traits$lma, h)
 dms<-dMsdA(traits$rho, h)
 dmb<-dMbdA(traits$rho, h)
 dmh<-dMhdA(traits$rho, h)
 dmr<-dMrdA(h)
 dmt<-dml+dmr+dmb+dms+dmh
  plot(h,dml, type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, max(dmt)*1.25), ylab="marginal cost of growth (kg/m)")
  points(h,dml, type='l', col='green')
  points(h,(dml+dmr), type='l', col='yellow')
  points(h,(dml+dmr+dmb), type='l', col='blue')
  points(h,(dml+dmr+dmb+dms), type='l', col='brown')
  }


#plot plotMarginalCosts
plotMarginalCostsLfArea<-function(traits, A){
dml<-dMldA(traits$lma, A)
dms<-dMsdA(traits$rho, A)
dmb<-dMbdA(traits$rho, A)
dmh<-dMhdA(traits$rho, A)
dmr<-dMrdA(A)
dmt<-dml+dmr+dmb+dms+dmh
  plot(A,dmt, type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, max(dmt)*1.25), ylab="marginal cost of growth (kg/m2)")
  points(A,dml, type='l', col='green')
  points(A,(dml+dmr), type='l', col='yellow')
  points(A,(dml+dmr+dmb), type='l', col='blue')
  points(A,(dml+dmr+dmb+dms), type='l', col='brown')
  }


#plot mass fractions
plotMassFraction<-function(traits, h){
  A= LeafArea(h)
  ml =A*traits$lma
  ms =SapwoodMass(traits$rho,A,h)
  mb =BarkMass(traits$rho,A,h) 
  mr =RootMass(A)
  mh= HeartwoodMass(traits$rho,A)
  mt = ml+ms+mb+mr+mh

  plot(h,ml/mt, type='l', ylim=c(0,1), col='green', xaxs="i", yaxs="i", ylab="fraction of an individual’s mass")
  points(h,(ml+mr)/mt, type='l', col='yellow')
  points(h,(ml+mr+mb)/mt, type='l', col='blue')
  points(h,(ml+mr+mb+ms)/mt, type='l', col='brown')
}


#plot mass fractions
plotProductionComponents<-function(traits, h, env){
  A= LeafArea(h)
  ml =A*traits$lma
  ms =SapwoodMass(traits$rho,A,h)
  mb =BarkMass(traits$rho,A,h)
  mr =RootMass(A)
  
  
  plot(h, h*0, type='n', ylim=c(0,2.5), col='black', xaxs="i", yaxs="i", ylab= "Income or cost (kg m^-2 yr^-1)")

  for(i in 1:length(env))
    points(h, Assim(A,env[i])/A, type='l',col='black')
  
  Cost =   Respiration.leaf(A)
  points(h,Cost/A, type='l', lty="dashed", col='green')
  Cost =  Cost +  Turnover.leaf(traits$lma, ml)
  points(h,Cost/A, type='l', lty="solid", col='green')
    
  Cost =  Cost+ Respiration.root(mr)
  points(h,Cost/A, type='l', lty="dashed", col='yellow')
  Cost =  Cost +  Turnover.root(mr)
  points(h,Cost/A, type='l', lty="solid", col='yellow')

  Cost =  Cost+ Respiration.bark(mb/traits$rho)
  points(h,Cost/A, type='l', lty="dashed", col='blue')
  Cost =  Cost +  Turnover.bark(mb)
  points(h,Cost/A, type='l', lty="solid", col='blue')

  Cost =  Cost+ Respiration.sapwood(ms/traits$rho)
  points(h,Cost/A, type='l', lty="dashed", col='brown')
  Cost =  Cost +  Turnover.sapwood(ms)
  points(h,Cost/A, type='l', lty="solid", col='brown')
}
#plot    ReproductiveAllocation
plotReproductiveAllocation<-function(traits, h, add=FALSE){
  if(add==FALSE)
    plot(h,   ReproductiveAllocation(traits$hmat,h), type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, 1), ylab="fraction mass allocated to reproduction", xlab="heigth (m)")
  else
    points(h,   ReproductiveAllocation(traits$hmat,h), type='l', col='black') 
  }

#plot Height Growth v size
plotSizevHeightGrowth<-function(traits, h, env=1.0, option = 1, add=FALSE, RGR=FALSE){
  nLines<-length(env)
  #solve height growth rate for env and height combination
  y<-matrix(0,length(h),nLines)
  for(i in 1:nLines){y[,i]<-dHdt(traits, h, env[i])}  
  
  YLAB="height growth rate (m/yr)"
  if(RGR==TRUE){YLAB="relative height growth rate (m/m/yr)"; y<-y/h;}
  
  #plot all curves
   if(add==FALSE)
    matplot(cbind(h),y, type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, max(y)*1.25),ylab=YLAB, xlab="height (m)")
 else  
    matpoints(cbind(h),y, type='l',col='black')
  }

#plot WPLCP v size
plotSizevWPLCP<-function(traits, h, option = 1, add=FALSE){
  nLines<-length(h)
  #solve height growth rate for env and height combination
  y<-matrix(0, length(h),1)
  for(i in 1:nLines){y[i,1]<-WPLCP(traits, h[i])}  
  #plot all curves
   if(add==FALSE)
    matplot(cbind(h),y, type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, max(y)*1.25),ylab="WPLCP (0-1)", xlab="height (m)")
 else  
    matpoints(cbind(h),y, type='l',col='black')
  }


#plot plot Trait v Height Growth: 0= seed mass, 1 = LMA, 2= WOOD DENSITY
plotTraitvHeightGrowth<-function(x, traits, h, env=1.0, option = 1, add=FALSE, ...){

  if(option==0){
  XLAB = "SEED MASS (kg)";
  h<-Height.mt(traits, x)
  nLines<-length(env)
  y<-matrix(0,length(x),nLines)

  for(i in 1:nLines){
    for(j in 1:length(x)){y[j,i]<-dHdt(traits, h[j], env[i])}}
  } 

if(option>0){ #not seed mass
  #make new object with variation in focal trait
  if(option==1){ XLAB = "LMA (kg/m2)";}
  if(option==2){XLAB = "RHO (kg/m3)";}
  if(option==3){XLAB = "LA:SA (m2/m2)"; theta = p.theta}

  
  nLines<-max(length(h), length(env))
  if(length(h)<nLines){h<-rep(h[1], nLines)}
  if(length(env)<nLines){env<-rep(env[1], nLines)}
  
  #solve height growth rate for env and height combination
  y<-matrix(0,length(x),nLines)
  for(i in 1:nLines){
    for(j in 1:length(x)){
      if(option==1){traits$lma<-x[j]}
      if(option==2){traits$rho<-x[j]}
      if(option==3){p.theta <-x[j]; cat(p.theta)}
      y[j,i]<-dHdt(traits, h[i], env[i])}}
  }
  
  #plot all curves
  if(add==FALSE)
    matplot(cbind(x),y, type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, max(y)*1.25),ylab="height growth rate (m/yr)", xlab=XLAB, ...)
 else  
    matpoints(cbind(x),y, type='l',col='black', ...)
  
  if(option==3){p.theta=theta}
  }

#plot plot Trait v Height Growth: 1 = LMA, 2= WOOD DENSITY, 3=theta
plotTraitvWPLCP<-function(x, traits, h=1, option = 1, ...){

if(option==0){
  XLAB = "SEED MASS (kg)";
  h<-Height.mt(traits, x)
  y<-matrix(0,length(x),1)
  for(j in 1:length(x)){y[j,1]<-WPLCP(traits, h[j])}
  } 

if(option>0){ #not seed mass
  if(option==1){XLAB = "LMA (kg/m2)";}
  if(option==2){XLAB = "RHO (kg/m3)";}
  if(option==3){XLAB = "SA:LA (m2/m2)"; theta = p.theta}
  
  nLines<-max(length(h))
  if(length(h)<nLines){h<-rep(h[1], nLines)}
  
  #solve height growth rate for env and height combination
  y<-matrix(0,length(x),nLines)
  for(i in 1:nLines){
    for(j in 1:length(x)){
      if(option==1){traits$lma<-x[j]}
      if(option==2){traits$rho<-x[j]}
      if(option==3){p.theta <-x[j]}
 
      y[j,i]<-WPLCP(traits, h[i])
    }
  }
  }
  #plot all curves
  matplot(cbind(x),y, type='l', log="x",col='black', xaxs="i", yaxs="i", ylim=c(0, 0.5),ylab="Light compensation point", xlab=XLAB, ...)

 if(option==3){p.theta=theta}
}


