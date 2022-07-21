\name{plotKmeansClustering}
\alias{plotKmeansClustering}
\title{ Kmeans Clustering Plot }
\description{ 
  \code{plotKmeansClustering} represents, in different subpanels, each of the clusters obtained after running k-means. The corresponding centroid is highlighted.
}

\usage{
plotKmeansClustering(x, kmeansObj, col=c(8,2), lty=c(2,1), x.coord = NULL, 
    no.ticks = 5, ...)
}

\arguments{
  \item{x}{ a data matrix containing \code{N} observations (individuals) by rows and \code{d} variables (features) by columns }
  \item{kmeansObj}{ an object of class \code{kmeans}, containing the cluster labels output by kmeans } 
  \item{col}{ a vector containing colors for the elements in \code{x} and for the centroid. The last one is used for the centroid, whereas the previous ones are recycled }
  \item{lty}{ a vector containing the line type for the elements in \code{x} and for the centroid. The last one is used for the centroid, whereas the previous ones are recycled }
  \item{x.coord}{ initial x coordinates (time points) where the functional data is observed; if not provided, it is assumed to be \code{1:d} }
  \item{no.ticks}{ number of ticks to be displayed in the X axis }
  \item{...}{ additional arguments to be passed to the \code{plot} function }
}

\value{ the function returns invisibly a list with the following components:
  \item{clusters}{ a list containing one cluster per component; observations are given by rows }
  \item{centroids}{ a list with the centroid of each cluster }
}

\details{ The function creates a suitable grid where to plot the different clusters independently. In the i-th cell of the grid, the data points corresponding to the i-th cluster are represented in parallel coordinates and the final centroid is highlighted. }

\author{ Javier Albert Smet \email{javas@kth.se} and 
Aurora Torrente \email{etorrent@est-econ.uc3m.es}}

\examples{
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

    plotKmeansClustering(x, kmeans(x,4))
    plotKmeansClustering(x, brik(x,4)$km)
    plotKmeansClustering(x, fabrik(x,4)$km)
    plotKmeansClustering(x, fabrik(x,4,degFr=10)$km)
}
\keyword{kmeans}
\keyword{bootstrap}
\keyword{MBD}
