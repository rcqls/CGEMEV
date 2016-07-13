grid.domain  <- function(missing.domains,grid.size=256,oversampling.size=2.5,dim=2) {
	obj <- new.env()
	obj$missing.domains <- missing.domains # disk only now! TODO: extend to any form of domain
	obj$n.domains <- length(missing.domains)

	obj$dim <- dim
	if(length(grid.size)==1) grid.size <- rep(grid.size,dim)
	obj$grid.size <- grid.size
	obj$n1 <- obj$grid.size[1] #TODO: extend to other domain
	obj$coords <- (1./(obj$n1-1)) * matrix( c(rep(1:obj$n1,obj$n1),rep(1:obj$n1,each=obj$n1)), ncol=2) # dim or 2
	obj$oversampling.size <- oversampling.size
	class(obj) <- "grid.missing.domain"
	update(obj)
	obj
}

update.grid.domain <-function(obj) {
	missing.sites.grid.domain(obj)
	distant.sites.grid.domain(obj)
	non.missing.coords.grid.domain(obj)
	non.missing.number.grid.domain(obj)
}

missing.sites.grid.domain <- function(obj) {
	n1Xn1  <- obj$n1*obj$n1
	missing.sites <- rep(FALSE,  n1Xn1)
	for(indexPixel in 1:n1Xn1) for(id in 1:obj$n.domains) missing.sites[indexPixel] <- missing.sites[indexPixel] | any(belong.to.disk(obj$coords[indexPixel,] ,  obj$missing.domains[[id]]$half.center,obj$missing.domains[[id]]$radius^2 ))
	# OLD: listOfbelong.to.OneOfTheDisks
	obj$missing.sites <- missing.sites
}

distant.sites.grid.domain <- function(obj) {
	n1Xn1  <- obj$n1*obj$n1

	roundHalfWidth  <-  round(obj$oversampling.size/obj$n1)
	halfWidthNeighborood  <-  obj$oversampling.size/obj$n1
	augmentedRadius  <-  sapply(obj$missing.domains,function(md) md$half.center + halfWidthNeighborood)

	distant.sites <- rep(FALSE,  n1Xn1)
	for(i in (roundHalfWidth+1):(n1-roundHalfWidth)) for(j in (roundHalfWidth+1):(n1-roundHalfWidth))  {
		indexPixel <- j + (i-1)*obj$n1
		for(id in 1:obj$n.domains) distant.sites[indexPixel] <- distant.sites[indexPixel] | !any(belong.to.disk(obj$coords[indexPixel,] ,  obj$missing.domains[[id]]$half.center,augmentedRadius[id]^2 ))
	}
	# OLD: listOfDistantEnoughPixels.from.BoundaryAndTheDisks
	obj$distant.sites <- distant.sites
}

# OLD: observationSitesWithMissingDisk
non.missing.coords.grid.domain <- function(obj) obj$non.missing.coords <- obj$coords[!obj$missing.sites,]

# OLD: n
non.missing.number.grid.domain <- function(obj) obj$non.missing.number <- sum(!obj$missing.sites) #if full grid then ==obj$n1*obj$n1

expand.to.full.grid.domain  <- function(obj,z) {
	## it should be checked that length(z)+ (nb of true booleans) = n1*n1
	zWithMissingEqNA  <- matrix(0,nrow=obj$grid.size[1],ncol=obj$grid.size[2])
	zWithMissingEqNA[!obj$missing.sites]  <- z
  # arrayOfIndicesOfObservationInVectorZ
  obj$z.obs <- matrix(zWithMissingEqNA,nrow=obj$grid.size[1],ncol=obj$grid.size[2])
}
