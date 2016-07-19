#include <Rcpp.h>
using namespace Rcpp;

double matern(NumericVector d, range = 1, nu = 0.5) {

    d  = d / range; //iterator
    d[d == 0] = 1e-10;//+test abs < 1e-10

    con = pow(2,nu - 1) * R::gamma(nu);
    con = 1./con;
    return con * power(d,nu) * R::bessel_k(d, nu, 1);
}

// [[Rcpp::export]]
List preconditioning_info(environment obj,Function first_preconditioning, double nu, double range, double precond.bandwidth) {
  std::vector<double> entriesRaw;
  std::vector<int> colindicesRaw,rowpointersRaw;

  int is_first_precond=true;
  List first_precond;

  int iz,length_listPi;
  IntegerVector j;
  double gi;

  int icol,
  for(icol=0;icol<obj[["grid.size"]][1];icol++) for(irow=0;irow<obj[["grid.size"]][0];irow++) {
    iz=obj[["z"]](irow, icol);
    if (! obj[["missing.sites"]](irow, icol)) {
    	if (obj[["distant.sites"]](irow, icol)) {
    		if(is_first_precond)) {
          is_first_precond=false;
          first_precond=first_preconditioning(irow,icol);
        }
                   
          length.listPi=first_precond[["l"]];
          //TODO!!!
          j <- obj[["z"]][first_precond[[lag.coords.listPi]][0] + 
          matrix(c(rep(irow,length.listPi),rep(icol,length.listPi)),ncol=2)]

          gi=first_precond[["g"]];
    		} else {
          //TO REDO!!!
          listPi = as.vector(pdist( obj$non.missing.coords[iz,],
          obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2): iz,])@dist)<(precond.bandwidth/obj$n1)     
          j <- which(listPi)
          j <- j + max(1, iz - 2*obj$n1 - 2) - 1
          //END REDO

          if ( length(j)==1 ) {
            gi=1
            length_listPi = 1
    		} else {
          int length_listPi=length(j);
          //compute all Euclidian distance
          NumericMatrix distSites(length_listPi,length_listPi);
          // rdist(obj$non.missing.coords[j,])
          NumericMatrix subR(length_listPi,length_listPi);
          subR=matern(distSites,range,nu);
          NumericVector vectorEmi(length.listPi);
          vectorEmi[length(vectorEmi)-1]=1;
          //Armadillo
          NumericVector gTildai=solve(subR, vectorEmi);
          //To normalize
          gi <- gTildai/Math.sqrt(gTildai[length_listPi-1]);

    		}

    	}
      //TODO!!
    	entriesRaw <-c(entriesRaw,gi)
    	colindicesRaw <-c(colindicesRaw,j)
      //cumlength.listPi <- cumlength.listPi + length.listPi
      rowpointersRaw <-c(rowpointersRaw, rowpointersRaw[length(rowpointersRaw)] +length_listPi)
    	//rowpointersRaw[iz+1] <- cumlength.listPi

    }
  }
  return List::create(_["entriesRaw"]=entriesRaw,_["colindicesRaw"]=colindicesRaw,_["rowpointersRaw"]=rowpointersRaw);
}
// preconditioning.info <- function(obj,first.precond,first.preconditioning,nu,range,precond.bandwidth) {
//   entriesRaw<-c()
//   colindicesRaw<-c()
//   rowpointersRaw <- c(1)
//
//   for(icol in 1:obj$grid.size[2]) for(irow in 1:obj$grid.size[1]) {
//     iz <- obj$z[irow, icol]
//     if (! obj$missing.sites[ irow, icol ]) {
//     	# (! iz ==0)
//     	if (obj$distant.sites[ irow, icol ]) {
//     		if(is.null(first.precond)) first.precond <- first.preconditioning(irow,icol)
//                             #
//           length.listPi <- first.precond$l
//           j <- obj$z[first.precond$lag.coords.listPi + 
//           matrix(c(rep(irow,length.listPi),rep(icol,length.listPi)),ncol=2)]
//           gi <- first.precond$g
//     		} else { #print(iz)
//           listPi <- as.vector(pdist( obj$non.missing.coords[iz,],
//           obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2): iz,])@dist)<(precond.bandwidth/obj$n1) 
//                                   #
//           j <- which(listPi)
//           j <- j + max(1, iz - 2*obj$n1 - 2) -1
//           if ( length(j)==1 ) {
//           gi <- 1
//           length.listPi <- 1
//     		} else {  
//           subR <- Matern(rdist(obj$non.missing.coords[j,]), nu=nu, range=range)
//           length.listPi <- length(j)
//           vectorEmi <- c(rep(0,length.listPi -1),1)
//           gTildai <- solve(subR, vectorEmi)
//           gi <- gTildai/sqrt(gTildai[length.listPi])
//
//     		}
//
//     	}
//     	entriesRaw <-c(entriesRaw,gi)
//     	colindicesRaw <-c(colindicesRaw,j)
//       #cumlength.listPi <- cumlength.listPi + length.listPi
//       rowpointersRaw <-c(rowpointersRaw, rowpointersRaw[length(rowpointersRaw)] +length.listPi)
//     	#rowpointersRaw[iz+1] <- cumlength.listPi
//
//     }
//   }
//   return(list(entriesRaw,colindicesRaw,rowpointersRaw))
// }
