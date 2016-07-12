cd Documents/testFSAIinR
R
library(fields)
#setting for simulation of a Gaussian Matern field on a grid n1*n1
n1<-256
#    (x_1,x_2)(i), i=1,..,n1*n1     equispaced on the unit square:
xFull<- (1./(n1-1))*matrix( c(rep(1: n1, n1),rep(1: n1,each= n1)), ncol=2)
n1Xn1  <- n1*n1
#each realization of the field will be observed on a given uncomplete grid
#### initialization of the array of booleans
####  listOfbelong.to.OneOfTheDisks
#### by the following "source" :
##source("design-setting2missingDisks.R")
########## design setting in cases of a few missing disks
### initialize the booleans-array listOfbelong.to.OneOfTheDisks     ###
### and the booleans-array  listOfDistantEnoughPixels.from.BoundaryAndTheDisks
### eaxmaple :
listOfHalfDiameters  <-  c(0.17841  ,  0.071365)
listOfHalfCenters  <-  c(    c(0.67841,0.67841),    c(0.278412, 0.228412) )
listOfAugmentedfHalfDiameters  <-  listOfHalfDiameters + rep(2.5/n1, 2)
##########
belong.to.disk  <- function(x, center, halfDiameter.power2) {
	dif <- x -center
	t(dif) %*% dif  <halfDiameter.power2 }
#######
listOfbelong.to.OneOfTheDisks <- rep(FALSE,  n1Xn1)
for(indexPixel in 1:n1Xn1) {
	listOfbelong.to.OneOfTheDisks[indexPixel] <-any(
	c(
	belong.to.disk(xFull[indexPixel,] ,  c(0.67841,0.67841),  (0.17841)^2) ,
	belong.to.disk(xFull[indexPixel,] ,  c(0.278412, 0.228412),  (0.071365)^2)
	))
}
###############################################################
nbOfDisks <- 3
listOfHalfDiameters  <-  c(0.17841  ,  0.071365,  0.03568)
listOfHalfCenters  <- t( matrix(    c(c(0.67841,0.67841),    c(0.278412, 0.228412),    c(0.8, 0.25) ), ncol= nbOfDisks) )
halfWidthNeighborood  <-  2.5/n1
#halfWidthNeighborood  <-  0.2
listOfAugmentedfHalfDiameters  <-  listOfHalfDiameters + rep(halfWidthNeighborood, nbOfDisks)
###################
####### precomputation of booleans array  for the distant-enough sites ##################
roundHalfWidth <- round(halfWidthNeighborood*n1)
listOfDistantEnoughPixels.from.BoundaryAndTheDisks <- rep(FALSE,  n1Xn1)
for(i in (roundHalfWidth+1):(n1-roundHalfWidth)){
for(j in (roundHalfWidth+1):(n1-roundHalfWidth))  {
	indexPixel <-j +(i-1)*n1
	listOfDistantEnoughPixels.from.BoundaryAndTheDisks[indexPixel] <-!any(
	c(
	belong.to.disk(xFull[indexPixel,] ,  listOfHalfCenters[1,],  (listOfAugmentedfHalfDiameters[1])^2) ,
	belong.to.disk(xFull[indexPixel,] , listOfHalfCenters[2,],  (listOfAugmentedfHalfDiameters[2])^2)
	))
}}
########################
###################
observationSitesWithMissingDisk <- xFull[!listOfbelong.to.OneOfTheDisks,]
n<- nrow(observationSitesWithMissingDisk)    #if full grid then n=n1*n1
> n
[1] 57989
### the following function will be used in the "PCG via FFT" and in
### precomputation of the preconditioning
###################### expand to the full grid (0 at missing) ##########
expand.to.fullGrid  <- function(n1, z, listOfbelong.to.OneOfTheDisks) {
	## it should be checked that length(z)+ (nb of true booleans) = n1*n1
	zWithMissingEqNA  <- matrix(0,nr=n1,nc=n1)
	zWithMissingEqNA[!(listOfbelong.to.OneOfTheDisks)]  <- z
	matrix(zWithMissingEqNA, nr=n1) }
#######
### 1 example of parameters setting:
bTrue<- 1.0    #variance (or "power") of the field
# noise level =0
nu<-1./2    #this “smoothness index” will be assumed known
#
trueTheta=10.
trueRange<-1./trueTheta
########
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
utilisateur     système      écoulé
      1.369       0.072       1.531
set.seed(321)  # so that it is reproducible #
system.time(look<- sim.rf( obj) )
########### plot of this realization
x1 <- x2 <- seq(1, n1, 1)
trueZ <-look[1:n1,1:n1]
image(x1,x2,trueZ)
utilisateur     système      écoulé
      0.548       0.033       0.694
 ##### only observed outside the disks:
prov <- c(trueZ)
z  <-   prov[!listOfbelong.to.OneOfTheDisks]
length(z)
sum(z**2)/n
[1] 1.330449
### precomputation of the preconditioning-matrix sparseG :
source("~/Documents/testFSAIinR/precondSettingWithMissingDisks.R")
# Message d'avis :
# In pdist(observationSitesWithMissingDisk[indexOfPixelInZ, ], observationSitesWithMissingDisk[max(1,  :
  # Y is the same as X, did you mean to use dist instead?
# > ut
   # user  system elapsed
# 352.688  35.782 388.394
G  <-  newMoreMoreSpeeded.G.efficientComput
transpOfG <-  t(G)

candidateThetas.Grid<- (10.) * 10**seq(-1.1,1.1,,15)
source("~/Documents/testFSAIinR/conjugate.gradient.R")
#
source("~/Documents/testFSAIinR/fsai11Precond.GEevalOnThetaGridfinalFinal.R")
out <- fsai11Precond.GEevalOnThetaGrid(n1,z, candidateThetas.Grid,nu, listOfbelong.to.OneOfTheDisks , G, transpOfG,tolPGC=1e-03)
out
# > out
# $values
            # [,1]
 # [1,] 12.5005024
 # [2,]  8.7055279
 # [3,]  6.0626384
 # [4,]  4.2221937
 # [5,]  2.9406180
 # [6,]  2.0482494
 # [7,]  1.4269924
 # [8,]  0.9946361
 # [9,]  0.6939372
# [10,]  0.4852058
# [11,]  0.3409401
# [12,]  0.2424073
# [13,]  0.1773076
# [14,]  0.1386134
# [15,]  0.1245319

# $niterForY
 # [1]  45  24  37  37  34  33  29  18  28  18  23  37  58 104 201
#
