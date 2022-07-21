\name{fabrik}
\alias{fabrik}
\title{ Computation of Initial Seeds for Kmeans and Clustering of Functional Data }
\description{ 
  \code{fabrik} fits splines to the multivariate dataset and runs the BRIk algorithm on the smoothed data. For functional data, this is just a straight forward application of BRIk to the k-means algorithm; for multivariate data, the result corresponds to an alternative clustering method where the objective function is not necessarily minimised, but better allocations are obtained in general.
}

\usage{
fabrik(x, k, method="Ward", nstart=1, B = 10, J = 2, x.coord = NULL, OSF = 1, 
    vect = NULL, intercept = TRUE, degPolyn = 3, degFr = 5, knots = NULL, ...)
}

\arguments{
  \item{x}{ a data matrix containing \code{N} observations (individuals) by rows and \code{d} variables (features) by columns }
  \item{k}{ number of clusters } 
  \item{method}{ clustering algorithm used to cluster the cluster centres from the bootstrapped replicates; \code{Ward}, by default. Currently, only \code{pam} and randomly initialised \code{kmeans} are implemented }
  \item{nstart}{ number of random initialisations when using the \code{kmeans} method to cluster the cluster centres }
  \item{B}{ number of bootstrap replicates to be generated }
  \item{J}{ number of observations used to build the bands for the MBD computation. Currently, only the value J=2 can be used }
  \item{x.coord}{ initial x coordinates (time points) where the functional data is observed; if not provided, it is assumed to be \code{1:d} }
  \item{OSF}{ oversampling factor for the smoothed data; an OSF of m means that the number of (equally spaced) time points observed in the approximated function is m times the number of original number of features, \code{d} }
  \item{vect}{ optional collection of x coordinates (time points) where to assess the smoothed data; if provided, it ignores the OSF }
  \item{intercept}{ if \code{TRUE}, an intercept is included in the basis; default is \code{FALSE} }
  \item{degPolyn}{ degree of the piecewise polynomial; 3 by default (cubic splines) }
  \item{degFr}{ degrees of freedom, as in the \code{bs} function }
  \item{knots}{ the internal breakpoints that define the spline }
  \item{...}{ additional arguments to be passed to the \code{kmeans} function for the final clustering; at this stage \code{nstart} is set to 1, as the initial seeds are fixed }
}

\value{
  \item{seeds}{ a matrix of size \code{k x D}, where \code{D} is either \code{m x d} or the length of \code{vect} . It contains the initial smoothed seeds obtained with the BRIk algorithm }
  \item{km}{ an object of class \code{kmeans} corresponding to the run of kmeans on the smoothed data, with starting points \code{seeds} }
}

\details{ The FABRIk algorithm extends the BRIk algorithm to the case of longitudinal functional data by adding a step that includes B-splines fitting and evaluation of the curve at specific x coordinates. Thus, it allows handling issues such as noisy or missing data. It identifies smoothed initial seeds that are used as starting points of kmeans on the smoothed data. The resulting clustering does not optimise the distortion (sum of squared distances of each data point to its nearest centre) in the original data space but it provides in general a better allocation of datapoints to real groups. }

\author{ Javier Albert Smet \email{javas@kth.se} and 
Aurora Torrente \email{etorrent@est-econ.uc3m.es}}

\references{ 
Torrente, A. and Romo, J. (2020). Initializing Kmeans Clustering by Bootstrap and Data Depth. \emph{J Classif} (2020). https://doi.org/10.1007/s00357-020-09372-3.
Albert-Smet, J., Torrente, A. and Romo J. (2021). Modified Band Depth Based Initialization of Kmeans for Functional Data Clustering. Submitted to Computational Statistics and Data Analysis.
}
\examples{
## fabrik algorithm 
    ## simulated data
    set.seed(1)
    x.coord = seq(0,1,0.01)
    x <- matrix(ncol = length(x.coord), nrow = 100)
    labels <- matrix(ncol = 100, nrow = 1)
  
    centers <-  matrix(ncol = length(x.coord), nrow = 4)
    centers[1, ] <- abs(x.coord)-0.5
    centers[2, ] <- (abs(x.coord-0.5))^2 - 0.8
    centers[3, ] <- -(abs(x.coord-0.5))^2 + 0.7
    centers[4, ] <- 0.75*sin(8*pi*abs(x.coord))
  
    for(i in 1:4){
        for(j in 1:25){
            labels[25*(i-1) + j] <- i  
            if(i == 1){x[25*(i-1) + j, ] <- abs(x.coord)-0.5 + 
                rnorm(length(x.coord),0,1.5)}
            if(i == 2){x[25*(i-1) + j, ] <- (abs(x.coord-0.5))^2 - 0.8 + 
                rnorm(length(x.coord),0,1.5)}
            if(i == 3){x[25*(i-1) + j, ] <- -(abs(x.coord-0.5))^2 + 0.7 + 
                rnorm(length(x.coord),0,1.5)}
            if(i == 4){x[25*(i-1) + j, ] <- 0.75*sin(8*pi*abs(x.coord)) + 
                rnorm(length(x.coord),0,1.5)}
            }
        }

    C1 <- kmeans(x,4)
    C2 <- fabrik(x,4,B=25)
  
    table(C1$cluster, labels)
    table(C2$km$cluster, labels)    
}

\keyword{kmeans}
\keyword{bootstrap}
\keyword{cluster} 
\keyword{functional data}
\keyword{MBD}
