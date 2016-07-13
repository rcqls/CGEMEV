library(CGEMEV)
library(fields)
library(pdist)

#setting for simulation of a Gaussian Matern field on a grid n1*n1

ex1.missing.domains <- list(
	list(half.center=c(0.67841,0.67841),radius=0.17841),
	list(half.center=c(0.278412, 0.228412),radius=0.071365)
)

ex2.missing.domains <- c(ex1.missing.domains,list(half.center=c(0.8, 0.25),radius=0.03568))


generate.missing.pixels <- function(missing.domains,grid.size=c(256,256),oversampling.size=2.5) {
	n1 <- grid.size[1] #TODO: extend to other domain
	xFull <- (1./(n1-1)) * matrix( c(rep(1: n1, n1),rep(1: n1,each= n1)), ncol=2)
	n1Xn1  <- n1*n1
	n.disks <- length(missing.domains$listOfHalfCenters)


	listOfbelong.to.OneOfTheDisks <- rep(FALSE,  n1Xn1)
	for(indexPixel in 1:n1Xn1) for(i.disk in 1:n.disks) listOfbelong.to.OneOfTheDisks[indexPixel] <- listOfbelong.to.OneOfTheDisks[indexPixel] | any(belong.to.disk(xFull[indexPixel,] ,  missing.domains[[i.disk]]$half.center,missing.domains[[i.disk]]$radius^2 ))
	listOfbelong.to.OneOfTheDisks
}

generate.distant.pixels <- function(missing.domains,grid.size=c(256,256),oversampling.size=2.5,roundHalfWidth) {
	n1 <- grid.size[1] #TODO: extend to other domain
	xFull <- (1./(n1-1)) * matrix( c(rep(1: n1, n1),rep(1: n1,each= n1)), ncol=2)
	n1Xn1  <- n1*n1
	n.disks <- length(missing.domains)

	## fix roundHalfWidth=0 when calling this function to choose no border
	if(missing(roundHalfWidth)) roundHalfWidth  <-  round(oversampling.size/n1)
	halfWidthNeighborood  <-  oversampling.size/n1
	augmentedfRadius  <-  sapply(missing.domains,function(md) md$half.center + halfWidthNeighborood)

	listOfDistantEnoughPixels.from.BoundaryAndTheDisks <- rep(FALSE,  n1Xn1)
	for(i in (roundHalfWidth+1):(n1-roundHalfWidth)) for(j in (roundHalfWidth+1):(n1-roundHalfWidth))  {
		indexPixel <- j + (i-1)*n1
		for(i.disk in 1:n.disks) listOfDistantEnoughPixels.from.BoundaryAndTheDisks[indexPixel] <- listOfbelong.to.OneOfTheDisks[indexPixel] | !any(belong.to.disk(xFull[indexPixel,] ,  missing.domains[[i.disk]]$half.center,augmentedRadius[i.disk]^2 )
	}
	listOfDistantEnoughPixels.from.BoundaryAndTheDisks
}

missing.pixels <- generate.missing.pixels(ex1.missing.domains)
distant.pixels <- generate.distant.pixels(ex1.missing.domains)


### the following function will be used in the "PCG via FFT" and in
### precomputation of the preconditioning
###################### expand to the full grid (0 at missing) ##########
expand.to.fullGrid  <- function(grid.size, z, listOfbelong.to.OneOfTheDisks) {
	if(length(grid.size)===1) grid.size <- c(grid.size,grid.size)
	## it should be checked that length(z)+ (nb of true booleans) = n1*n1
	zWithMissingEqNA  <- matrix(0,nr=grid.size[1],nc=grid.size[2])
	zWithMissingEqNA[!(listOfbelong.to.OneOfTheDisks)]  <- z
	matrix(zWithMissingEqNA,nr=grid.size[1],nc=grid.size[2])
}
#######

### 1 example of parameters setting:
bTrue<- 1.0    #variance (or "power") of the field
# noise level =0
nu <-1./2    #this “smoothness index” will be assumed known
#
trueTheta=10.
trueRange<-1./trueTheta
########
#simulate.missing.gaussian.matern <- function() {
#one simulates 1 realization using fctions of fields -package
factor <-1    #factor =2  or 4 or even 8 may be required  for some other settings
sizeGridsimul<-factor*n1
#Simulate a Gaussian random field with an exponential covariance function, 
#range parameter =1./10 and the domain is [0, factor]X [0, factor] evaluating the 
#field at a sizeGridsimul x sizeGridsimul grid. 
grid<- list( x= seq( 0, factor,, sizeGridsimul), y= seq(0, factor,, sizeGridsimul)) 
# system.time(obj<-stationary.image.cov( setup=TRUE, 
             # grid=grid, 
             # range= trueRange,smoothness=nu))
#
#XXXXXX attention en 2016: il faut utiliser theta au lieu de range et matern ds le nom de la fct
system.time(obj<-matern.image.cov( setup=TRUE, 
             grid=grid, 
             theta= trueRange,smoothness=nu))
# utilisateur     système      écoulé
#      1.369       0.072       1.531
set.seed(321)  # so that it is reproducible #
system.time(look <- sim.rf( obj) )
########### plot of this realization
x1 <- x2 <- seq(1, n1, 1)
trueZ <-look[1:n1,1:n1]
image(x1,x2,trueZ)
# only observed outside the disks:
prov <- c(trueZ)
z  <-   prov[!listOfbelong.to.OneOfTheDisks]
length(z)
sum(z**2)/n
### precomputation of the preconditioning-matrix sparseG :
### source("precondSettingWithMissingDisks.R")

listOfPixels.in.OneOfTheDisks  <- listOfbelong.to.OneOfTheDisks
#
prov <-  listOfPixels.in.OneOfTheDisks
arrayOfIndicesOfObservationInVectorZ  <- expand.to.fullGrid(n1,1:length(prov[!prov]), prov)
#
array.listOfPixels.in.OneOfTheDisks   <- matrix(listOfPixels.in.OneOfTheDisks,ncol=n1)
#
array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks   <- matrix(listOfDistantEnoughPixels.from.BoundaryAndTheDisks,ncol=n1)
listOfIndicesOfObservationInVectorZ <- c(arrayOfIndicesOfObservationInVectorZ)
#

ut<-system.time(
{
	preconditioning(arrayOfIndicesOfObservationInVectorZ,array.listOfPixels.in.OneOfTheDisks,c(64,64))

}
)

# Message d'avis :
# In pdist(observationSitesWithMissingDisk[indexOfPixelInZ, ], observationSitesWithMissingDisk[max(1,  :
  # Y is the same as X, did you mean to use dist instead?
# > ut
   # user  system elapsed
# 352.688  35.782 388.394
G  <-  newMoreMoreSpeeded.G.efficientComput
transpOfG <-  t(G)

candidateThetas.Grid <- (10.) * 10**seq(-1.1,1.1,,15)
#
out <- fsai11Precond.GEevalOnThetaGrid(n1,z, candidateThetas.Grid,nu, listOfbelong.to.OneOfTheDisks , G, transpOfG,tolPGC=1e-03)
out
