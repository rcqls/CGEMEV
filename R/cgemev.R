# sites
grid <- function(dim) {
  # TODO
}

# matrix (only called when size od$f sites is reasonable)
cor.matern.matrix <- function(sites,theta,b) { #notation ref à l'article de Didier

}

# model
matern.cov.field <- function(nu,theta,b=1) { # contructor
  obj <- list(nu=nu,theta=theta,b=b)
  class(obj) <- "matern.cov.field"
  obj
}

# simulation

simulate.matern.cov.field <- function(n,sites=grid(c(100,100)),missing.sites) {
## z=b^(.5)L^Tw (ou version sans calcul de R ni de L => circulante embedding but only for grid !!!!)
## return z
}

## Example
# ex.mat <- matern.cov.model(.5,1)
# simulate(ex.mat,500000)

## Estimation of theta and b

cgemev <- function(z,missing.sites,tol.bissection=1e-4,tol.pcg=1e-4) {
## theta_r+1=theta_r (z^T R(theta_r)^(-1) z/ z^Tz)^(1/2*nu)
## preconditionning
## G^T R(theta) Gw=G^T z => x=Gw
## with FSAI: see algo
}


code.didier <- function() {
  library(fields)
#  n is the length of the observed vector z
# the observation sites are the pixels FALSE in the array
#       of booleans array.listOfPixels.in.OneOfTheDisks
#
#  arrayOfIndicesOfObservationInVectorZ gives for each pixel i1,i2, the index (in z)
#   of the associated observation,  or the value 0 if
#   i1,i2 is not observed,
#
# the precomputed array
#       of booleans array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks
#  gives all the pixels where the inversion of the small sub-matrix of R would give
#  the same result So this result is precomputed and stored in the following vectors
preComputedlistPi <- c(2,   3,   4, 129, 130, 131, 132, 133, 257, 258, 259)
preComputedGi <- c(0.4588033,  0.6016878,  0.7916430,  0.5423615, -0.2983703, -3.8459627,-2.3750820, -1.0472091, -0.2293937, -5.2853303, 10.6979353)
length.precomputedListPi  <- length(preComputedlistPi)
n  <- length(which(! array.listOfPixels.in.OneOfTheDisks))
ut<-system.time(
{
junk.G <-  spam(0, n, n)
for(indexColumnInGrid in 1:n1){
for(indexRowInGrid in 1:n1)  {
	indexOfPixelInZ <- arrayOfIndicesOfObservationInVectorZ[indexRowInGrid, indexColumnInGrid]
#
	if  (! array.listOfPixels.in.OneOfTheDisks[ indexRowInGrid, indexColumnInGrid ])
    #	(! indexOfPixelInZ ==0)   marcherait aussi
    {
    	if (array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks[ indexRowInGrid, indexColumnInGrid ])
      {#print(c("row=",indexRowInGrid,"comumn=", indexColumnInGrid))
      	j  <- listOfIndicesOfObservationInVectorZ[preComputedlistPi +  ( indexRowInGrid  - 3) +
      	       n1*(indexColumnInGrid - 3)]
  	    #length.listPi <- length(j)
  	    i  <- rep(indexOfPixelInZ, length.precomputedListPi)
  	    gi <- preComputedGi}
     else
			  { 	#print(indexOfPixelInZ)
				# listPi <- as.vector(pdist(	observationSitesWithMissingDisk[indexOfPixelInZ,],
	# #			observationSitesWithMissingDisk[max(1, indexOfPixelInZ - 2*n1 - 2):min(n, indexOfPixelInZ + 2*n1 + 2)])@dist)<(2.5/n1)
				# observationSitesWithMissingDisk[max(1, indexOfPixelInZ - 2*n1 - 2): indexOfPixelInZ,])@dist)<(2.5/n1) #!!!!!!! magic !!!!
				# #
			  # #if (indexOfPixelInZ <n)  listPi[(indexOfPixelInZ +1):n] <- FALSE
			  # j <-  which(listPi)
			  # j <-  j  +  max(1, indexOfPixelInZ - 2*n1 - 2) -1
			  	# if   (  length(j)==1  )
			      # {i <- indexOfPixelInZ
				     # gi <-  1 }
			    # else
			      # {  #subR <-  exp(- rdist(observationSitesWithMissingDisk[j,]))
				     # subR <-  Matern(rdist(observationSitesWithMissingDisk[j,]), nu=0.5, range=1)
			       # length.listPi <- length(j)
			       # i  <- rep(indexOfPixelInZ, length.listPi)
			       # vectorEmi <- c(rep(0,length.listPi -1),1)
			       # gTildai <-  solve(subR, vectorEmi)
			       # gi <-   gTildai/sqrt(gTildai[length.listPi])  }
			   }
			junk.G [cbind(i, j)] <-  gi
}
}
}})

}
