# sites
grid <- function(dim) {
  # TODO
}

# model
cor.matern.matrix <- function(theta,b) { #notation ref à l'article de Didier

}

matern.cov.field <- function(nu,theta,b=1) { # contructor
  obj <- list(nu=nu,theta=theta,b=b)
  class(obj) <- "matern.cov.field"
  obj
}

# simulation

simulate.matern.cov.field <- function(n,sites=grid(c(100,100)),missing.sites) {
## z=b^(.5)L^Tw (ou version sans calcul de R ni de L => circulante embedding....)
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
