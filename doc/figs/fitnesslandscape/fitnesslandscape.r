#!/usr/bin/Rscript


loadFitnessDataFromFile <- function( filename=filename, zOffset = 0, zlim=c( -Inf, Inf ) ){ 
  
  Data<-read.table( filename, skip=1, header=FALSE )
  
  names( Data )<-switch(as.character(dim( Data )[2]),
                        "2" = c( "x", "z" ),
                        "3" = c( "x", "y", "z" ))
  
  #apply offset 
  Data$z= Data$z-zOffset;
  #trim within specified limits
  Data$z[Data$z > zlim[2]] =zlim[2]
  Data$z[Data$z < zlim[1]] =zlim[1]
  
  
  Data
}


plotFitnessLandscapeFromFile <- function(filename, dim=2, title="", ...){
  
  if( !file.exists( filename ) ){
    plot.new()   #file not available, plot blank
  } else {  
    
    zlim=c(-6,6)
    
    Data<-try(loadFitnessDataFromFile( filename=filename, zlim=zlim ), silent=TRUE)
    # use try, because sometimes load fails, if errors in file  
    if(inherits(Data, "try-error"))
      plot.new()   #file not available, plot blank
    else{
      Residents <- loadResidnetTraitFromFitnessFile( filename = filename)
      plotFitnessLandscape(Data, Residents, dim=dim, title=title,...)
    }
  }
} 


plotFitnessLandscape <- function(Data, Residents, title="", axisInfoFn=axisInfoFn, dim=2, 
                                 cex=1, outer=FALSE, line=3, zlim=c(-6,6)){
  
  if( dim==2 )
    plotFitnessSurf2D(Data,Residents, zlim=zlim)  
  else
    plotFitnessSurf1D(Data,Residents, dim, zlim=zlim)
  
  #add title
  mtext( title, side =3, line =1, cex=cex )
  
  #X axis - add ticks, tick labels, title
  axis( 1, at =-2:2, labels= powerTenExpression(-2:2) )
  mtext("Leaf mass per area ( kg/m2 )", side = 1, line=line, outer =outer, cex=cex )    
  
  #Yaxis - add ticks, tick labels, title
  hspTick = log10( c( 1,5,10,20,40,80 ) )
  axis( 2, at = hspTick, labels= 10^hspTick, las=1)
  mtext("Height at maturation (m)", side = 2, line=line, outer =outer, cex=cex )    
  
}

to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

# returns 10^x as expression. useful for labelling axes
powerTenExpression<-function(x){
  do.call(expression, lapply(x, function(i) bquote(10^.(i))))  
}



plotFitnessSurf1D <- function(Data, Residents, dim, zlim=c( -6,6 ),xlab="", ylab="", main="" ){
  # Plots trait vs fitness landscape for data stored in filename
  
  #Make plot
  plot( Data$x, Data$z, type='l',
        ylim=zlim,
        xlab=xlab, ylab=ylab, main=main,
        xaxs = "i", yaxs = "i", las = 1, ann=FALSE, xaxt="n", yaxt="n" )  
  
  xlim <- par( "usr" )[1:2]  #upper and lower limits of x-axis
  points( xlim, xlim*0, type='l', lty='dashed' )   
  
  #Add resident location
  #interpolate to get fitness at resident point
  points(Residents[,dim+1] ,  approx(Data$x, Data$z, Residents[,dim+1])$y , type='p', pch=21, col="black", bg="white", cex=1.5 )        
  
}


plotFitnessSurf2D <- function( Data, Residents, zlim=c( -6,6 ), add = FALSE ){
  
  Data<-collapse.grid( Data )
  
  #load colorMap
  myColorMap<-loadColorMap( "colormap.csv", zlim=zlim )
  
  #Make plot
  if(!add){
    plot.new( )
    plot.window( xlim= range( Data$x, finite = TRUE ), ylim= range( Data$y, finite = TRUE ), "", xaxs = "i", yaxs = "i", las = 1 )    
    box( )
  }

  .filled.contour( Data$x, Data$y, Data$z, levels= myColorMap$levels, col = myColorMap$col )
  
  #Add resident location
  points( Residents[,1],Residents[,2], type='p', pch=21, col="black", bg="white", cex=1.5 )
  
  return( myColorMap )
  #add residents  
}


loadColorMap <- function( filename="R/colormap.csv", zlim=c( -1,1 ) ){
  #Load custom colormap
  #read colormap
  tmp<-t( read.csv( filename, h=TRUE ) )
  
  #structure to store output
  out<-list( )
  out$col<- rgb( tmp[1,], tmp[2,], tmp[3,], alpha=1 )
  
  #Set number of levels to match colormap - should be one more
  nlevels <- length( out$col ) +1
  out$levels <- seq( from =zlim[1], to=zlim[2], length.out=nlevels ) 
  return( out )
}



loadResidnetTraitFromFitnessFile <- function( filename, ncol=3){
  con <- file( filename, "r", blocking = FALSE )
  Res.str<-strsplit( readLines( con,1 ), "%" )[[1]][2]
  close( con )
  Res<-as.numeric( unlist( strsplit( Res.str, "\t" ) ) )
  
  matrix( Res, ncol=ncol, nrow = length( Res )/3, byrow = TRUE )
}


is.wholenumber <-  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


collapse.grid<-function(Data){
  # Does reverse of expand.grid function
  # Takes a dataframe with up to 2 or 3 cols and converts into x,y,z vectors
  # such that Data[,c("x","y)] = expand.grid(x,y)
  
  #Find dimensions of X and Y vectors
  ncol<-match(FALSE, Data$y[2:length(Data$y)] > Data$y[1:(length(Data$y)-1)])
  nrow<-length(Data$x)/ncol
  #Extract x and y vectors
  out<-NULL
  out$y <-Data$y[1:ncol]
  out$x <-Data$x[ncol*seq(0,(nrow-1))+1]
  #Reshape Z data into matrix
  out$z<-matrix(Data$z, ncol =ncol, nrow=nrow, byrow=TRUE)
  out
}

setwd("figs/fitnesslandscape/")
to.pdf(plotFitnessLandscapeFromFile("T-7400-2D.txt", title = ""),
       filename = "plot.pdf",
       height = 6, width = 8)
  