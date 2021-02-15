\name{brik}
\alias{brik}
\title{ Computation of Initial Seeds and Kmeans Results }
\description{ 
  \code{brik} computes appropriate seeds --based on bootstrap and the MBD depth-- to initialise k-means, which is then run.
}

\usage{
brik(x, k, method="Ward", nstart=1, B=10, J = 2, ...)
}

\arguments{
  \item{x}{ a data matrix containing \code{N} observations (individuals) by rows and \code{d} variables (features) by columns }
  \item{k}{ number of clusters } 
  \item{method}{ clustering algorithm used to cluster the cluster centres from the bootstrapped replicates; \code{Ward}, by default. Currently, only \code{pam} and randomly initialised \code{kmeans} are implemented }
  \item{nstart}{ number of random initialisations when using the \code{kmeans} method to cluster the cluster centres }
  \item{B}{ number of bootstrap replicates to be generated }
  \item{J}{ number of observations used to build the bands for the MBD computation. Currently, only the value J=2 can be used }
  \item{...}{ additional arguments to be passed to the \code{kmeans} function for the final clustering; at this stage \code{nstart} is set to 1, as the initial seeds are fixed }
}

\value{
  \item{seeds}{ a matrix of size \code{k x d} containing the initial seeds obtained with the BRIk algorithm }
  \item{km}{ an object of class \code{kmeans} corresponding to the run of kmeans on \code{x} with starting points \code{seeds} }
}

\details{ The brik algorithm is a simple, computationally feasible method, which provides k-means with a set of initial seeds to cluster datasets of arbitrary dimensions. It consists of two stages: first, a set of cluster centers is obtained by applying k-means to bootstrap replications of the original data to be, next, clustered; the deepest point in each assembled cluster is returned as initial seeds for k-means. }

\author{ Javier Albert Smet \email{javas@kth.se} and 
Aurora Torrente \email{etorrent@est-econ.uc3m.es}}

\references{
Torrente, A. and Romo, J. (2020). Initializing k-means Clustering by Bootstrap and Data Depth. \emph{J Classif} (2020). https://doi.org/10.1007/s00357-020-09372-3.
}
\examples{
## brik algorithm 
    ## simulated data
    set.seed(0)
    g1 <- matrix(rnorm(200,0,3), 25, 8) ; g1[,1]<-g1[,1]+4;
    g2 <- matrix(rnorm(200,0,3), 25, 8) ; g2[,1]<-g2[,1]+4; g2[,3]<-g2[,3]-4
    g3 <- matrix(rnorm(200,0,3), 25, 8) ; g3[,1]<-g3[,1]+4; g3[,3]<-g3[,3]+4

    x <- rbind(g1,g2,g3)
    labels <-c(rep(1,25),rep(2,25),rep(3,25))

    C1 <- kmeans(x,3)
    C2 <- brik(x,3,B=25)
  
    table(C1$cluster, labels)
    table(C2$km$cluster, labels)    

}
\keyword{kmeans}
\keyword{bootstrap}
\keyword{MBD}
