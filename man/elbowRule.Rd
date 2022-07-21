\name{elbowRule}
\alias{elbowRule}
\title{ Selection of Appropriate DF Parameter Based on an Elbow Rule for the Distortion }
\description{ 
  \code{elbowRule} runs the FABRIk algorithm for different degrees of freedom (DF) and suggests the best of such values as the one where the minimum distortion is obtained. An optional visualization of the computed values allows the choice of alternative suitable DF values based on an elbow-like rule.
}

\usage{
elbowRule(x, k, method="Ward", nstart=1, B = 10, J = 2, x.coord = NULL, OSF = 1, 
    vect = NULL, intercept = TRUE, degPolyn = 3, degFr = 4:20, knots = NULL, 
    plot = FALSE, ...)
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
  \item{degFr}{ a vector containing tentative values of the degrees of freedom, to be tested }
  \item{knots}{ the internal breakpoints that define the spline }
  \item{plot}{ a Boolean parameter; it allows plotting the distortion against the degrees of freedom. Set to \code{FALSE} by default }
  \item{...}{ additional arguments to be passed to the \code{kmeans} function for the final clustering; at this stage \code{nstart} is set to 1, as the initial seeds are fixed }
}

\value{
  \item{df}{ the original vector of DF values to be tested }
  \item{tot.withinss}{ a vector containing the distortion obtained for each tested DF value }
  \item{optimal}{ DF value producing the smallest distortion among the tested \code{df} }
}

\details{ The function implements a simple elbow-like rule that allows selecting an appropriate value for the DF parameter among the tested ones. It computes the distortion obtained for each of these values and returns the one yielding to the smallest distortion. By setting the parameter \code{plot} to \code{TRUE} the distortion is plotted against the degrees of freedom and elbows or minima can be visually detected. }

\author{ Javier Albert Smet \email{javas@kth.se} and 
Aurora Torrente \email{etorrent@est-econ.uc3m.es}}

\references{ 
Torrente, A. and Romo, J. (2020). Initializing Kmeans Clustering by Bootstrap and Data Depth. \emph{J Classif} (2020). https://doi.org/10.1007/s00357-020-09372-3.
Albert-Smet, J., Torrente, A. and Romo J. (2021). Modified Band Depth Based Initialization of Kmeans for Functional Data Clustering. Submitted to Computational Statistics and Data Analysis.
}

\examples{
    ## simulated data
    set.seed(1)
    x.coord = seq(0,1,0.01)
    x <- matrix(ncol = length(x.coord), nrow = 80)
    labels <- matrix(ncol = 100, nrow = 1)
  
    centers <-  matrix(ncol = length(x.coord), nrow = 4)
    centers[1, ] <- abs(x.coord)-0.5
    centers[2, ] <- (abs(x.coord-0.5))^2 - 0.8
    centers[3, ] <- -(abs(x.coord-0.5))^2 + 0.7
    centers[4, ] <- 0.75*sin(8*pi*abs(x.coord))
  
    for(i in 1:4){
        for(j in 1:20){
            labels[20*(i-1) + j] <- i  
            if(i == 1){x[20*(i-1) + j, ] <- abs(x.coord)-0.5 + 
                rnorm(length(x.coord),0,1.5)}
            if(i == 2){x[20*(i-1) + j, ] <- (abs(x.coord-0.5))^2 - 0.8 + 
                rnorm(length(x.coord),0,1.5)}
            if(i == 3){x[20*(i-1) + j, ] <- -(abs(x.coord-0.5))^2 + 0.7 + 
                rnorm(length(x.coord),0,1.5)}
            if(i == 4){x[20*(i-1) + j, ] <- 0.75*sin(8*pi*abs(x.coord)) + 
                rnorm(length(x.coord),0,1.5)}
            }
        }

    # ER <- elbowRule(x, 4, B=25, degFr = 5:12, plot=FALSE)
    ER <- elbowRule(x, 4, B=25, degFr = 5:12, plot=TRUE)
  
}
\keyword{elbow rule}
\keyword{kmeans}
\keyword{bootstrap}
\keyword{cluster} 
\keyword{functional data}
\keyword{MBD}
