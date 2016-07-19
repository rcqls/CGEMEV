//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double matern(double d, double range=1, double nu=.5) {

    d  = d / range; //iterator
    if(fabs(d) < 1e-10) d=1e-10;//+test abs < 1e-10

    double con = pow(2,nu - 1) * gamma(nu);
    con = 1./con;
    return con * pow(d,nu) * R::bessel_k(d, nu, 1);
}

// [[Rcpp::export]]
List preconditioning_info(Environment obj,Function first_preconditioning, double nu, double range, double precond_bandwidth) {
  std::vector<double> entriesRaw;
  std::vector<int> colindicesRaw,rowpointersRaw,j;

  rowpointersRaw.push_back(1);
  int is_first_precond=true;
  List first_precond;
  NumericVector first_precond_g;
  NumericMatrix first_precond_lags;
  int first_precond_lags_size;

  int iz,length_listPi;
  std::vector<double> gi;

  int icol,irow;
  NumericVector grid_size=obj["grid.size"];
  NumericMatrix z=obj["z"];
  LogicalMatrix missing_sites=obj["missing.sites"],distant_sites=obj["distant.sites"];
  NumericMatrix non_missing_coords=obj["non.missing.coords"];

  int n1=obj["n1"];

  double dmax=precond_bandwidth/(double)n1;

  for(icol=0;icol<grid_size[1];icol++) for(irow=0;irow<grid_size[0];irow++) {
    iz=z(irow, icol)-1;
    //printf("iz=%d\n",iz);
    //printf("missing_sites(irow, icol)=%d\n",missing_sites(irow, icol));
    //printf("distant_sites(irow, icol)=%d\n",distant_sites(irow, icol));
    if (!missing_sites(irow, icol)) {
    	if (distant_sites(irow, icol)) {
        //printf("iciiiiiiiiii\n");
    		if(is_first_precond) {
          is_first_precond=false;
          first_precond=first_preconditioning(irow+1,icol+1);
          first_precond_g=first_precond["g"];
          first_precond_lags=as<NumericMatrix>(first_precond["lag.coords.listPi"]);
          first_precond_lags_size = first_precond_lags.nrow();
        }
         
        length_listPi=first_precond["l"];
        //printf("length_listPi=%d\n",length_listPi);

        //TODO ==> DONE
        // j <- z[first_precond[[lag.coords.listPi]][0] + 
        // matrix(c(rep(irow,length.listPi),rep(icol,length.listPi)),ncol=2)];
        //normally == length_listPi
        j.clear();
        //printf("vals=");
        for(int k=0;k<first_precond_lags_size;k++) {
          int val=z(irow+first_precond_lags(k,0),icol+first_precond_lags(k,1));
          //printf(" %d",val);
          j.push_back(val);
        }
        //printf("\n");
        gi.resize(length_listPi);
        gi.clear();
        gi.insert(gi.end(),first_precond_g.begin(),first_precond_g.end());
  		} else {
          //TO REDO!!!
          // listPi = as.vector(pdist( obj$non.missing.coords[iz,],
          // obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2): iz,])@dist)<(precond.bandwidth/obj$n1)     
          // j <- which(listPi)
          // j <- j + max(1, iz - 2*obj$n1 - 2) - 1
          //END REDO

          //max(1, iz - 2*obj$n1 - 2)
          int iz_from=iz-2*n1-2;
          if(iz_from<0) iz_from=0;
          //printf("iz_from=%d\n",iz_from);

          j.clear();
          //printf("dmax=%lf,rdist=",dmax);
          for(int iz_cur=iz_from;iz_cur<=iz;iz_cur++) {
            double d=sqrt(sum(pow(non_missing_coords(iz,_)-non_missing_coords(iz_cur,_),2)));
            //printf(" %lf",d);
            if(d<dmax) j.push_back(iz_cur+1);//R (index)
          }
          //printf("\n");

          if ( j.size()==1 ) {
            gi.resize(1);
            gi[0]=1;
            length_listPi = 1;
    		  } else {
            length_listPi=j.size();
            //compute all Euclidian distance
            //NumericMatrix distSites(length_listPi,length_listPi);
          // rdist(obj$non.missing.coords[j,])
          // subR=matern(distSites,range,nu);

          arma::mat subR(length_listPi,length_listPi);
          printf("length_listPi=%d,dist=",length_listPi);
          for(int k1=0;k1<length_listPi;k1++) {
            for(int k2=0;k2<k1;k2++) {
              double d=sqrt(sum(pow(non_missing_coords(j[k1]-1,_)-non_missing_coords(j[k2]-1,_),2)));
              printf("(%d,%d)=(%lf,",j[k1]-1,j[k2]-1,d);
              double val=matern(d,range,nu);
              printf("%lf) ",val);
              subR(k1,k2) = val;
              subR(k2,k1) = val;
            }
            subR(k1,k1)=1; //If I am not mistaken!
          }
          printf("\n");

          arma::colvec vectorEmi(length_listPi);
          vectorEmi[length_listPi-1]=1;
          //Armadillo
          //arma::vec gTildai=arma::inv(subR) * vectorEmi;
          arma::vec gTildai=arma::solve(subR,vectorEmi);
          gi.resize(length_listPi);
          double sqrtGTildaiLast=sqrt(gTildai[length_listPi-1]);
          for(int k=0;k<length_listPi;k++) gi[k]=gTildai(k)/sqrtGTildaiLast;
          // //To normalize
          // gi <- gTildai/Math.sqrt(gTildai[length_listPi-1]);
    		}

    	}
      //TODO!!
    	// entriesRaw <-c(entriesRaw,gi)
    	// colindicesRaw <-c(colindicesRaw,j)
      // //cumlength.listPi <- cumlength.listPi + length.listPi
      // rowpointersRaw <-c(rowpointersRaw, rowpointersRaw[length(rowpointersRaw)] +length_listPi)
    	// //rowpointersRaw[iz+1] <- cumlength.listPi
      entriesRaw.insert(entriesRaw.end(), gi.begin(), gi.end());
      colindicesRaw.insert(colindicesRaw.end(),j.begin(),j.end());
      rowpointersRaw.push_back(rowpointersRaw[rowpointersRaw.size()-1] + length_listPi);
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
