# #cd Documents/testFSAIinR-on-perturbed-grid-subSampled
# #R
# library(fields)
# library(spatstat)
# #setting for simulation of a Gaussian Matern field on 10000 locations from a perturbed grid 200*200
# n1<-200
# #each realization of the field will be observed at the sites
# #    (x_1,x_2)(i), i=1,..,n1*n1     equispaced on the unit square, next perturbed:
# xx<- (1./(n1-1))*matrix( c(rep(0: (n1-1), n1),rep(0: (n1-1),each= n1)), ncol=2)
# xx <- xx + (1./(n1-1))* matrix( runif(2*n1*n1,min=-0.4,max=0.4), ncol=2)
# ##  on change l'ordre
# n=10000
# randOrder <- sample(1:(n1*n1),n,replace=FALSE) # sample(n1*n1,n,replace=FALSE)
# x <- xx[randOrder,]
# n<- nrow(x)
# x1 <- x[,1]
# x2  <- x[,2]
#
# points <- ppp(x=x1,
#          y=x2,
#          window=owin(xrange=c(-1./(n1-1),1+1./(n1-1)), yrange=c(-1./(n1-1),1+1./(n1-1))))
# plot(points)
# # true values for the 4 parameters :
# bTrue<- 1.0    #variance (or power) of the field
# # noise level =0
# trueRange<-1./16
# nu<-1./2    #this “smoothness index” will be assumed known
# #
# lags<-(1./(n1-1))*(0:(n1-1))
# autocorrelationFunct<- Matern( lags, range= trueRange,
# 	smoothness=nu)
# plot( lags, autocorrelationFunct, type="l")
#
# system.time( distMatrix <- rdist(x) )
# ##   user  system elapsed
# ##   2.269   0.479   2.747
# sharededge <- function(X) {
#   verifyclass(X, "ppp")
#   Y <- X[as.rectangle(X)]
#   dX <- deldir(Y)
#   DS <- dX$dirsgs
#   xyxy <- DS[,1:4]
#   names(xyxy) <- c("x0","y0","x1","y1")
#   sX <- as.psp(xyxy,window=dX$rw)
#   marks(sX) <- 1:nobjects(sX)
#   sX <- sX[as.owin(X)]
#   tX <- tapply(lengths.psp(sX), marks(sX), sum)
#   jj <- as.integer(names(tX))
#   ans <- data.frame(ind1=DS[jj,5],
#                    ind2=DS[jj,6],
#                    leng=as.numeric(tX))
#   return(ans)
# }
# system.time(shared_edge_lengths <- sharededge(points))
# ##   user  system elapsed
# ##  2.531   0.078   2.609
# incidenceMatrix <- matrix(rep(0,n*n),ncol=n)
#
# system.time(incidenceMatrix[cbind(shared_edge_lengths[[1]], shared_edge_lengths[[2]])] <-  1 )
# system.time(incidenceMatrix <- incidenceMatrix +t(incidenceMatrix) +diag(n) )
# ##   user  system elapsed
# ##  2.951   1.155   4.105
#
#
# ### tentative de "re-gridder" les locations (x1,x2)_i  ,  i=1,...,n  par la numerotation "reverse Cuthill-McKee (RCM)"
#   system.time( sparseVersionOfincidenceMatrix <- as.spam( incidenceMatrix)  )
# ##   user  system elapsed
# ##  2.364   0.581   2.945
#  sparseVersionOfPower3OfIncidenceMatrix <- sparseVersionOfincidenceMatrix%*%
#  (sparseVersionOfincidenceMatrix %*% sparseVersionOfincidenceMatrix)
# system.time(U <- chol(as.spam( sparseVersionOfPower3OfIncidenceMatrix +100*diag(n)), pivot = "RCM"))
# ##   user  system elapsed
# ##  5.442   2.691  12.224
#   ord <- ordering(U)
#   ## verification:
#   first.subset <- ord[1:100]
#   points <- ppp(x=x[first.subset,1],
#          y=x[first.subset ,2],
#          window=owin(xrange=c(-1./(n1-1),1+1./(n1-1)), yrange=c(-1./(n1-1),1+1./(n1-1))))
# plot(points)
#   second.subset <- ord[101:300]
#   points <- ppp(x=x[second.subset,1],
#          y=x[ second.subset ,2],
#          window=owin(xrange=c(-1./(n1-1),1+1./(n1-1)), yrange=c(-1./(n1-1),1+1./(n1-1))))
# plot(points)
#
#
# system.time(orderedDistanceMatrix <- distMatrix[ord, ord])
# ##  user  system elapsed
# ## 6.023   0.325   6.346
# ###
# ######
# ############ sparness pattern from sparseVersionOfPower3OfIncidenceMatrix
# ######### calcul de G à partir de RforRangeeq1000
# system.time(fullVersionOfincidencePower3 <- as.matrix(sparseVersionOfPower3OfIncidenceMatrix[ord,ord]))
# ##  user  system elapsed
# ##  5.225   0.487   5.711
# system.time(ordered.corr.matrixForRange1000 <- exp(-orderedDistanceMatrix/1000) )
# ##   user  system elapsed
# ##  1.902   1.105   3.019
#
#
#  ## calcul du "sparse factor of an approximate inverse" of ordered.corr.matrixForRange1000
#  ## algorithme décrit ds Janna et al. 2015 (Journal of the ACM)
#   ut<-system.time(
# {
# 	 G <-  spam(0, n, n)
#  for(indexRow in 1:n){
# listPi <- fullVersionOfincidencePower3[indexRow,]>0 #!!!!! also magic !!!!
# if (indexRow <n)  listPi[(indexRow +1):n] <- FALSE
#  j <-  which(listPi)
# subR <-  ordered.corr.matrixForRange1000[listPi, listPi]
# length.listPi <- length(j)
# i  <- rep(indexRow, length.listPi)
# vectorEmi <- c(rep(0,length.listPi -1),1)
# gTildai <-  solve(subR, vectorEmi)
#  gi <-  gTildai/sqrt(gTildai[length.listPi])
#  G [cbind(i, j)] <-  gi
# }
# })
# ut
#  ##  user  system elapsed
#  ## 68.051  46.873 115.126
# summary(G)
# # Matrix object of class 'spam' of dimension 10000x10000,
#     # with 212702 (row-wise) nonzero elements.
#     # Density of the matrix is 0.213%.
#     ## Class 'spam'
# display(G)
# ### donc environ 20 coef par lignes    C'est très sparse !!!
#
# ###### regardons si ce nouveau G améliore le nbre d'iterations
# # simulons un z
# ut<-system.time(
# {
# Correl.mat<-   exp(- distMatrix/trueRange)
# })
# ut
#    # user  system elapsed
#   # 2.296   0.978   7.097
# ut<-system.time(
# {
# A<- chol( Correl.mat)
# })
# ut
#    # user  system elapsed
# # 164.200   0.418 164.594
# ########## si on utilise la librairie optimisée (pour les multicoeurs sous MacOS X)
#    # # # user  system elapsed
#  # # # 35.152   0.604  10.238
#
#  #one simulates 1 realization
# ut<-system.time(z<- t(A) %*% rnorm(n) )
# ut
# ########## si on utilise la librairie optimisée (pour les multicoeurs sous MacOS X)
#       # # # user  system elapsed
#      # # # 2.054   0.330   2.229
#
#
# transpOfA <- t(A)
# ut<-system.time(z<- transpOfA %*% rnorm(n) )
# ut
# ########## si on utilise la librairie optimisée (pour les multicoeurs sous MacOS X)
#      # # # user  system elapsed
#      # # #  0.450   0.001   0.208
#           # # #  donc c'est la transposition qui coute cher
# sum(z**2)/n
# #[1] 1.036223
#
#
# ### appliquons GEEV matching:
# ##source("conjugate.gradient.R")
# approx.prod.DeltaTCorrelDelta.Timesx <- function( x) {
# 	result<- t(G)   %*%  x
# 	   result<- RmatrixCandidate  %*%  result
#    result<-  G %*%  result
# #   result<- bEV*result + x
#    result
# }
# #  ATENTIOn : ne pas oublier de se placer ds le probleme avec reordering:
# orderedZ <- z[ord]
# system.time(orderedDistanceMatrix <- distMatrix[ord,ord] )
#    # user  system elapsed
#   # 6.160   0.403   8.777
# system.time( {
#  	bEV<-sum(orderedZ**2)/n
#  	DeltaTy<-   G  %*%  orderedZ
#
# #varSignal<-bEV
# nINNERiter <- matrix(NA,51)
# rangeSuccessiveIterates <- matrix(NA,51)
#  	 theta <- 1./1000
#  	 rangeSuccessiveIterates[ 1] <-  1/theta
#  	 factCorrectif <-  2.;  coefProvFory<- DeltaTy;nbIter <- 0
#  	while(abs(factCorrectif - 1) > 0.001)
#  		{nbIter <- nbIter +1;
#  		if (nbIter ==51) break;
#  		RmatrixCandidate <- exp(-theta * orderedDistanceMatrix)
# 		   out2<- conjugate.gradient(DeltaTy,
#             approx.prod.DeltaTCorrelDelta.Timesx,
#             start= coefProvFory,
#             tol=1e-03,    kmax=200,verbose=FALSE)
#    nINNERiter[nbIter]<-out2$conv$niter
#    # save for the next CG-solver :
#    coefProvFory<- out2$x
#  #
# 		ge.quadratic.term <-  sum(out2$x * DeltaTy )/n;
# 		  			factCorrectif <-  ge.quadratic.term/bEV;
#     		theta <-theta * factCorrectif;
#     		rangeSuccessiveIterates[ nbIter +1] <-  1/theta;}
#  #   to improve accuracy of b theta one does'nt apply the final corection on theta:
#   theta <-theta / factCorrectif
# #	bEVthetaCGEMcCGEMEVnbIters[, indexReplcitate] <-  c(bEV,theta,bEV*theta,nbIter)
# })
# #    user  system elapsed
# #  25.178   3.409 108.427
# # nINNERiter[1:12]
# #  [1] 22 69  3 NA NA NA NA NA NA NA NA NA
# # > rangeSuccessiveIterates[1:6]
# # [1] 1.000000e+03 5.663412e-02 5.544583e-02 5.539854e-02           NA
# # [6]           NA
# # > bEV/rangeSuccessiveIterates[1:6]
# # [1]  0.00089496 15.80248836 16.14115815 16.15493651          NA          NA
