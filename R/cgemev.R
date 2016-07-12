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

# listOfPixels.in.OneOfTheDisks  <- listOfbelong.to.OneOfTheDisks
# #
# prov <-  listOfPixels.in.OneOfTheDisks
# arrayOfIndicesOfObservationInVectorZ  <-
# expand.to.fullGrid(n1,1:length(prov[!prov]), prov)
# #
# array.listOfPixels.in.OneOfTheDisks   <- matrix(listOfPixels.in.OneOfTheDisks,ncol=n1)
# #
# array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks   <- matrix(listOfDistantEnoughPixels.from.BoundaryAndTheDisks,ncol=n1)
# listOfIndicesOfObservationInVectorZ <- c(arrayOfIndicesOfObservationInVectorZ)
# #
################# precomputations : ##########
init.preconditioning <- function(arrayOfIndicesOfObservationInVectorZ,indexColumnInGrid,indexRowInGrid,grid.size=c(256,256),nu=.5,range=1,precond.bandwidth=2.5)  {
  ## TODO: extend grid.size[1] grid.size[2]
  library(pdist)
	indexOfPixelInZ <- arrayOfIndicesOfObservationInVectorZ[indexRowInGrid, indexColumnInGrid]
#
	if (! array.listOfPixels.in.OneOfTheDisks[ indexRowInGrid, indexColumnInGrid ]) {
#ceci marche aussi :  	(! indexOfPixelInZ ==0)
#(! listOfPixels.in.OneOfTheDisks[  ])
		 	#print(indexOfPixelInZ)
			listPi <- as.vector(pdist(	observationSitesWithMissingDisk[indexOfPixelInZ,],
			observationSitesWithMissingDisk)@dist)<(precond.bandwidth/grid.size[1])    #!!!!!!! magic !!!!
		 if (indexOfPixelInZ <n)  listPi[(indexOfPixelInZ +1):n] <- FALSE
		 j <-  which(listPi)
			#subR <-  Matern(rdist(j),nu=0.5,range=1)
			#
		 if   (  length(j)==1)
			{i <- indexOfPixelInZ
				gi <-  1
				}
		 else
			{  #subR <-  exp(- rdist(observationSitesWithMissingDisk[j,]))
				subR <-  Matern(rdist(observationSitesWithMissingDisk[j,]), nu=nu, range=range)
			length.listPi <- length(j)
			vectorEmi <- c(rep(0,length.listPi -1),1)
			gTildai <-  solve(subR, vectorEmi)
			 gi <-   gTildai/sqrt(gTildai[length.listPi])
			  }
		#i  <- rep(indexOfPixelInZ, length.listPi)
		#new.speededG [cbind(i, j)] <-  gi
}
list(i=indexOfPixelInZ,j=j,g=gi,l=length.listPi)
}


preconditioning <- function(arrayOfIndicesOfObservationInVectorZ,grid.size=c(256,256),nu=.5,range=1,precond.bandwidth=2.5)  {
  newMoreMoreSpeeded.G.efficientComput <-  spam(0, n, n)
  entriesRaw<-c()
  colindicesRaw<-c()
  rowpointersRaw <-c(1)
  n1 <- grid.size[1]
  first.precond <- NULL #found when first condition satisfied

  for(indexColumnInGrid in 1:n1) for(indexRowInGrid in 1:n1)  {
  	indexOfPixelInZ <- arrayOfIndicesOfObservationInVectorZ[indexRowInGrid, indexColumnInGrid]
  #
    if  (! array.listOfPixels.in.OneOfTheDisks[ indexRowInGrid, indexColumnInGrid ]) {
      #	(! indexOfPixelInZ ==0)
      #(! listOfPixels.in.OneOfTheDisks[  ])
      if (array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks[ indexRowInGrid, indexColumnInGrid ]) {
        #print(c("row=",indexRowInGrid,"comumn=", indexColumnInGrid))
        if(is.null(first.precond)) first.precond <- init.preconditioning(arrayOfIndicesOfObservationInVectorZ,indexRowInGrid, indexColumnInGrid,grid.size,nu,range,precond.bandwidth)

        j  <- listOfIndicesOfObservationInVectorZ[first.precond$j +  rep(( indexRowInGrid  - 3) +
        	       n1*(indexColumnInGrid - 3),  first.precond$l)]
    	  length.listPi <- length(j)
    	  i  <- rep(indexOfPixelInZ, length.listPi)
    	  gi <- first.precond$g
      } else { 	#print(indexOfPixelInZ)
  				listPi <- as.vector(pdist(	observationSitesWithMissingDisk[indexOfPixelInZ,],
  	#			observationSitesWithMissingDisk[max(1, indexOfPixelInZ - 2*n1 - 2):min(n, indexOfPixelInZ + 2*n1 + 2)])@dist)<(2.5/n1)
  				observationSitesWithMissingDisk[max(1, indexOfPixelInZ - 2*n1 - 2): indexOfPixelInZ,])@dist)<(2.5/n1)     #!!!!!!! magic !!!!
  				#
  			  #if (indexOfPixelInZ <n)  listPi[(indexOfPixelInZ +1):n] <- FALSE
  			  j <-  which(listPi)
  			  j <-  j  +  max(1, indexOfPixelInZ - 2*n1 - 2) -1
  			  if   (  length(j)==1  ) {
            i <- indexOfPixelInZ
  				  gi <-  1
  				  length.listPi <- 1
          } else {  #subR <-  exp(- rdist(observationSitesWithMissingDisk[j,]))
  				      subR <-  Matern(rdist(observationSitesWithMissingDisk[j,]), nu=0.5, range=1)
                length.listPi <- length(j)
                vectorEmi <- c(rep(0,length.listPi -1),1)
                gTildai <-  solve(subR, vectorEmi)
                gi <-   gTildai/sqrt(gTildai[length.listPi])
          }
  		}
      entriesRaw <-c(entriesRaw,gi)
      colindicesRaw <-c(colindicesRaw,j)
      rowpointersRaw <-c(rowpointersRaw, rowpointersRaw[length(rowpointersRaw)] +length.listPi)
    }
  }

  newMoreMoreSpeeded.G.efficientComput@entries <- entriesRaw
  newMoreMoreSpeeded.G.efficientComput@colindices <- as.integer(colindicesRaw)
  newMoreMoreSpeeded.G.efficientComput@rowpointers <- as.integer(rowpointersRaw)
#			newMoreMoreSpeeded.G.efficientComput[cbind(as.integer(colindicesRaw), as.integer(rowpointersRaw))] <-  entriesRaw
}

cgemev <- function(z,missing.sites,tol.bissection=1e-4,tol.pcg=1e-4) {
## theta_r+1=theta_r (z^T R(theta_r)^(-1) z/ z^Tz)^(1/2*nu)
## preconditionning
## G^T R(theta) Gw=G^T z => x=Gw
## with FSAI: see algo
}
