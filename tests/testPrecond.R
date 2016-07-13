library(CGEMEV)

# md=missing.domains
ex1.md <- list(
	list(half.center=c(0.67841,0.67841),radius=0.17841),
	list(half.center=c(0.278412, 0.228412),radius=0.071365)
)

ex2.md  <- c(ex1.md,list(half.center=c(0.8, 0.25),radius=0.03568))

# gd=grid.domain
ex1.gd <- grid.domain(ex1.md,64)

## Maybe shorten the names of the following vars as in the previous functions
#listOfbelong.to.OneOfTheDisks <- ex1.gd$
#listOfDistantEnoughPixels.from.BoundaryAndTheDisks <- generate.distant.sites(ex1.missing.domains)

### 1 example of parameters setting:
########


cat("gaussian matern creation\n")
gm <- gaussian.matern()
cat("-> done\n")
set.seed(321)  # so that it is reproducible #

cat("simulation and plot\n")
simulate(gm)
plot(gm)
cat("-> done\n")

# only observed outside the disks:
prov <- gm$look[1:gm$n1,1:gm$n1]
z  <-   prov[!ex1.gd$missing.sites]
length(z)

n <- ex1.gd$missing.sites.number
sum(z**2)/n

### precomputation of the preconditioning-matrix sparseG :
listOfPixels.in.OneOfTheDisks  <- listOfbelong.to.OneOfTheDisks
prov <-  listOfPixels.in.OneOfTheDisks
arrayOfIndicesOfObservationInVectorZ  <- expand.to.fullGrid(gm$n1,1:sum(!prov), prov)
array.listOfPixels.in.OneOfTheDisks   <- matrix(listOfPixels.in.OneOfTheDisks,ncol=gm$n1)
array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks   <- matrix(listOfDistantEnoughPixels.from.BoundaryAndTheDisks,ncol=gm$n1)


ut<-system.time({
	G <- preconditioning(arrayOfIndicesOfObservationInVectorZ,array.listOfPixels.in.OneOfTheDisks,array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks,c(64,64))
})

transpOfG <-  t(G)

candidateThetas.Grid <- (10.) * 10**seq(-1.1,1.1,,15)
#
out <- fsai11Precond.GEevalOnThetaGrid(n1,z, candidateThetas.Grid,nu, listOfbelong.to.OneOfTheDisks , G, transpOfG,tolPGC=1e-03)
out
