

fsai11Precond.GEevalOnThetaGrid <- function(n1,z,
#w,  #noisy case
candidateThetas1DGrid,nu=0.5, listOf.belongTo.OneOfTheDisks, 
sparseG,transOfSparseG,tolPGC=1e-03)
{
GEvaluesOnthegrid<- matrix(NA,length(candidateThetas1DGrid))
#trace.term<- matrix(NA,length(candidateThetas1DGrid))
niter <- matrix(NA,length(candidateThetas1DGrid),2)
n=length(z)
DeltaTy <-   sparseG  %*%  z
##########################
# the matrice-vector product required by conjugate.gradient:
bEV<-sum(z**2)/n 
viaFFTwithMissings.prod.DeltaTCorrelDelta.Timesx <- function( x) {
   x <- transOfSparseG   %*%  x
   xprovFull<- expand.to.fullGrid(n1,  x, listOf.belongTo.OneOfTheDisks)
   result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
   result<-  sparseG %*%  result[! listOf.belongTo.OneOfTheDisks]
   result }
##########################
coefProvFory<- z
# loop over ranges  :
for(k in 1:length(candidateThetas1DGrid)){
   tethaCand <- candidateThetas1DGrid[k]
   #XXXXXX attention en 2016: il faut utiliser theta au lieu de range et matern 
   #XXXXXX ds le nom de la fct
   cov.obj<- matern.image.cov( setup=TRUE, 
            grid=list(x=1:n1,y=1:n1), 
            theta= (1/tethaCand)*(n1-1),smoothness=nu)
   #
   out2<- conjugate.gradient(DeltaTy,  multAx=
            viaFFTwithMissings.prod.DeltaTCorrelDelta.Timesx,  
            start= coefProvFory,
            tol= tolPGC,    kmax=200,verbose=FALSE)   
   niter[k,1]<-out2$conv$niter
   coefProvFory<- out2$x  #save the sol for "start" for the next theta
   #
   GEvaluesOnthegrid[k] <- sum(DeltaTy *  coefProvFory) /n       
   }
list(values = GEvaluesOnthegrid,  
            niterForY=niter[,1]   
            #, niterForW=niter[,2],trace.term=trace.term #noisy case
            )
}
########### end of fsai11Precond.GEevalOnThetaGrid ##################
#####################################################################