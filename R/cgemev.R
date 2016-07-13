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


################# precomputations : ##########



#### INPUTS
## arrayOfIndicesOfObservationInVectorZ -> obj$z
## array.listOfPixels.in.OneOfTheDisks -> obj$missing.sites
## array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks -> obj$distant.sites

preconditioning.grid.domain <- function(obj,nu=.5,range=1,precond.bandwidth=2.5)  {
  library(pdist)
  library(fields)

  # function with no arguments since all variables inside inherited from the parent environment
  # only for code more readable
  first.preconditioning <- function()  {
    iz <- obj$z[irow, icol]

  	listPi <- as.vector(pdist(	obj$non.missing.coords[iz,],
  	obj$non.missing.coords)@dist)<(precond.bandwidth/grid.size[1])    #!!!!!!! magic !!!!
    if (iz <n)  listPi[(iz +1):n] <- FALSE
    j <-  which(listPi)
    #subR <-  Matern(rdist(j),nu=0.5,range=1)
    #
    if(length(j)==1) {
      i <- iz
    	gi <-  1
    } else {
      #subR <-  exp(- rdist(obj$non.missing.coords[j,]))
    	subR <-  Matern(rdist(obj$non.missing.coords[j,]), nu=nu, range=range)
      length.listPi <- length(j)
      vectorEmi <- c(rep(0,length.listPi -1),1)
      gTildai <-  solve(subR, vectorEmi)
      gi <-   gTildai/sqrt(gTildai[length.listPi])
    }
    list(i=iz,j=j,g=gi,l=length.listPi)
  }

  ## init sparse matrix
  sparseG <-  spam(0, obj$non.missing.number, obj$non.missing.number)
  entriesRaw<-c()
  colindicesRaw<-c()
  rowpointersRaw <-c(1)

  first.precond <- NULL #found when first condition satisfied

  for(icol in 1:obj$grid.size[2]) for(irow in 1:obj$grid.size[1])  {
  	iz <- obj$z[irow, icol]

    if  (! obj$missing.sites[ irow, icol ]) {
      #	(! iz ==0)
      #(! listOfPixels.in.OneOfTheDisks[  ])
      if (obj$distant.sites[ irow, icol ]) {
        #print(c("row=",irow,"comumn=", icol))
        if(is.null(first.precond)) first.precond <- first.preconditioning()

        j  <- obj$z[first.precond$j +  rep(( irow  - 3) + obj$n1*(icol - 3),  first.precond$l)]
    	  length.listPi <- length(j)
    	  i  <- rep(iz, length.listPi)
    	  gi <- first.precond$g
      } else { 	#print(iz)
  				listPi <- as.vector(pdist(	obj$non.missing.coords[iz,],
  	#			obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2):min(n, iz + 2*obj$n1 + 2)])@dist)<(2.5/obj$n1)
  				obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2): iz,])@dist)<(precond.bandwidth/obj$n1)     #!!!!!!! magic !!!!
  				#
  			  #if (iz <n)  listPi[(iz +1):n] <- FALSE
  			  j <-  which(listPi)
  			  j <-  j  +  max(1, iz - 2*obj$n1 - 2) -1
  			  if   (  length(j)==1  ) {
            i <- iz
  				  gi <-  1
  				  length.listPi <- 1
          } else {  #subR <-  exp(- rdist(obj$non.missing.coords[j,]))
  				      subR <-  Matern(rdist(obj$non.missing.coords[j,]), nu=nu, range=range)
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

  sparseG@entries <- entriesRaw
  sparseG@colindices <- as.integer(colindicesRaw)
  sparseG@rowpointers <- as.integer(rowpointersRaw)

  # Alternative:
  # If rowindicesRaw instead of rowpointersRaw
  # sparseG[cbind(as.integer(colindicesRaw), as.integer(rowindicesRaw))] <-  entriesRaw
  sparseG
}

cgemev <- function(z,missing.sites,tol.bissection=1e-4,tol.pcg=1e-4) {
## theta_r+1=theta_r (z^T R(theta_r)^(-1) z/ z^Tz)^(1/2*nu)
## preconditionning
## G^T R(theta) Gw=G^T z => x=Gw
## with FSAI: see algo
}
