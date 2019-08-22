rm(list = ls())
getwd()

library( MeanShift )

sample <- read.csv("/Users/jsgr/Downloads/sample/USER1_1.txt", sep = " ", skip = 1, header = FALSE)

my.files <- list.files(path = "/Users/jsgr/Downloads/sample/", pattern = "USER*", full.names = TRUE)

my.data <- lapply(my.files, read.csv, sep = " ", header = FALSE, skip = 1)

x.matrix <- lapply((lapply(my.data, '[[', 1)), as.double)
y.matrix <- lapply((lapply(my.data, '[[', 2)), as.double)
t.matrix <- lapply((lapply(my.data, '[[', 3)), as.double)
button.list <- lapply((lapply(my.data, '[[', 4)), as.double)
azimuth.list <- lapply((lapply(my.data, '[[', 5)), as.double)
altitude.list <- lapply((lapply(my.data, '[[', 6)), as.double)
pressure.list <- lapply((lapply(my.data, '[[', 7)), as.double)

plot( x.matrix[[1]], y.matrix[[1]], type="o", pch=16, main="Type 1 signature", xlab="x", ylab="y" )

signatures.labels <- rep( 1:5, rep( 40, 5 ) )

max.grid.length <- max( sapply( t.matrix, length ) )

grid.length = 2^(ceiling(log(max.grid.length)/log(2)))

## standardize x and y coordinates
standardize <- function( x )
{
  range <- range( x )
  output <- ( x - range[1] ) / diff( range )
  return( output )
}

x.list <- lapply( x.matrix, standardize )
y.list <- lapply( y.matrix, standardize )
t.list <- lapply( t.matrix, standardize )
  
wave.x <- mapply( projectCurveWavelets, x=t.list, y=x.list, MoreArgs=list( grid.length=grid.length, levels=8 ), SIMPLIFY=FALSE )
wave.y <- mapply( projectCurveWavelets, x=t.list, y=y.list, MoreArgs=list( grid.length=grid.length, levels=8 ), SIMPLIFY=FALSE )
wave.button <- mapply( projectCurveWavelets, x=t.list, y=button.list, MoreArgs=list( grid.length=grid.length, levels=8 ), SIMPLIFY=FALSE )
wave.azimuth <- mapply( projectCurveWavelets, x=t.list, y=azimuth.list, MoreArgs=list( grid.length=grid.length, levels=8 ), SIMPLIFY=FALSE )
wave.altitude <- mapply( projectCurveWavelets, x=t.list, y=altitude.list, MoreArgs=list( grid.length=grid.length, levels=8 ), SIMPLIFY=FALSE )
wave.pressure <- mapply( projectCurveWavelets, x=t.list, y=pressure.list, MoreArgs=list( grid.length=grid.length, filter.number=4, levels=8 ), SIMPLIFY=FALSE )

## wavelet coefficients
extractCoefficients <- function( x )
{
  output <- sapply( x , "[[", "coefficients" )
  return( output )
}

## combine wavelet objects into a unique list
wave.list <- list( wave.x, wave.y, wave.button, wave.azimuth, wave.altitude, wave.pressure )

## get coefficients list
wave.coefficients <- lapply( wave.list, extractCoefficients )

## combine coefficients into a unique matrix
wave.coefficients <- do.call( rbind, wave.coefficients )

## 3066 wavelet coefficients: ( 512 - 1 ) * 6 for each one of the 100 signatures
dim( wave.coefficients )

## note that the matrix is sparse
## print the proportion of non-zero wavelet coefficient for each curve
round( apply( wave.coefficients, 2, function( x ){ mean( x != 0 ) } ), 2 )

## bandwidth candidates
h.cand <- quantile( dist( t( wave.coefficients ) ), seq( 0.01, 0.17, by=0.01 ) )

## clustering using the blurring mean shift algorithm
system.time( clustering <- lapply( h.cand, function( h ){ bmsClustering( wave.coefficients, h=h ) } ) )

lapply( lapply( lapply( clustering, "[[", "labels" ), table ), sort, decreasing=TRUE )

index <- 2

## compare cluster labels and true labels
cluster1 <- which( clustering[[index]]$labels == 1 ) # USER2 21:40
print( cluster1 )
print( signatures.labels[cluster1] )

cluster2 <- which( clustering[[index]]$labels == 4 ) # USER1 1:20
print( cluster2 )
print( signatures.labels[cluster2] )

cluster3 <- which( clustering[[index]]$labels == 15 ) # USER3 40:60
print( cluster3 )
print( signatures.labels[cluster3] )

cluster4 <- which( clustering[[index]]$labels == 37 ) # USER5 80:100
print( cluster4 )
print( signatures.labels[cluster4] )

cluster5 <- which( clustering[[index]]$labels == 67 ) # USER4 60:80
print( cluster5 )
print( signatures.labels[cluster5] )

## modal coefficients for interesting clusters
## (each column corresponds to a feature)
mode1.coeffs <- matrix( clustering[[index]]$components[,1], nrow=1023 )
mode2.coeffs <- matrix( clustering[[index]]$components[,4], nrow=1023 )
mode3.coeffs <- matrix( clustering[[index]]$components[,15], nrow=1023 )
mode4.coeffs <- matrix( clustering[[index]]$components[,37], nrow=1023 )
mode5.coeffs <- matrix( clustering[[index]]$components[,67], nrow=1023 )

## put in a list
modal.coeffs <- list( mode1.coeffs, mode2.coeffs, mode3.coeffs, mode4.coeffs, mode5.coeffs )

## we need "wd" objects on which to apply the inverse DWT:
## we can extract the wd object from the wave.xxx lists!
wd.x.object <- wave.x[[1]]$y.wdT
wd.y.object <- wave.y[[1]]$y.wdT

## we used all of the 6 features for clustering, but
## for visualization we only care about the projection
## of the functional modes on the x-y plane!
##
## wd.button.object <- wave.button[[1]]$y.wdT
## wd.azimuth.object <- wave.azimuth[[1]]$y.wdT
## wd.altitude.object <- wave.altitude[[1]]$y.wdT
## wd.pressure.object <- wave.pressure[[1]]$y.wdT

## function for inverse DWT
invDWT <- function( wd.object, modal.coefficients )
{
  wd.object$D <- modal.coefficients
  output <- wr( wd.object )
  return( output )
}

modes.x <- vector( mode="list", length=5 )
modes.y <- vector( mode="list", length=5 )
for( i in 1:5 ){
  
  modes.x[[i]] <- invDWT( wd.x.object, modal.coeffs[[i]][,1] )
  modes.y[[i]] <- invDWT( wd.y.object, modal.coeffs[[i]][,2] )
  
}

## cluster 1
plot( NULL, main="Cluster 1", xlab="x", ylab="y", xlim=c( 0, 1 ), ylim=c( -0.1, 1.1 ) )
for( i in cluster1 )
{
  points( x.list[[i]], y.list[[i]], type="o", pch=16, cex=0.5, lwd=0.5 )
}
lines( modes.x[[1]], modes.y[[1]], col=2, lwd=4 )

## cluster 2
plot( NULL, main="Cluster 2", xlab="x", ylab="y", xlim=c( 0, 1 ), ylim=c( 0, 1 ) )
for( i in cluster2 )
{
  points( x.list[[i]], y.list[[i]], type="o", pch=16, cex=0.5, lwd=0.5 )
}
lines( modes.x[[2]], modes.y[[2]], col=2, lwd=4 )

## cluster 3
plot( NULL, main="Cluster 3", xlab="x", ylab="y", xlim=c( -0.1, 1 ), ylim=c( -0.4, 1 ) )
for( i in cluster3 )
{
  points( x.list[[i]], y.list[[i]], type="o", pch=16, cex=0.5, lwd=0.5 )
}
lines( modes.x[[3]], modes.y[[3]], col=2, lwd=4 )

## cluster 4
plot( NULL, main="Cluster 4", xlab="x", ylab="y", xlim=c( -0.1, 1 ), ylim=c( -0.2, 1 ) )
for( i in cluster4 )
{
  points( x.list[[i]], y.list[[i]], type="o", pch=16, cex=0.5, lwd=0.5 )
}
lines( modes.x[[4]], modes.y[[4]], col=2, lwd=4 )

## cluster 5
plot( NULL, main="Cluster 5", xlab="x", ylab="y", xlim=c( -0.1, 1.1 ), ylim=c( -0.1, 1 ) )
for( i in cluster5 )
{
  points( x.list[[i]], y.list[[i]], type="o", pch=16, cex=0.5, lwd=0.5 )
}
lines( modes.x[[5]], modes.y[[5]], col=2, lwd=4 )
