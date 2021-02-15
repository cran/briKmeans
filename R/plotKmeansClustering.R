plotKmeansClustering <- function (x, kmeansObj, col=c(8,2), lty=c(2,1),
    x.coord = NULL, no.ticks = 5, ...)
{
    if (class(kmeansObj)!="kmeans") {
        stop("Requires a kmeans object...")}
    if (nrow(x)!=length(kmeansObj$cluster)) {
        stop("Clustering does not correspond to data...")}
    if (is.null(x.coord)) {x.coord <- 1:ncol(x)}
    else {if (length(x.coord) != ncol(x)) {
        stop("Wrong number of coordinates...")
        }
    }
    
    noClust <- length(kmeansObj$size)
    gridR <- ceiling(sqrt(noClust))
    gridC <- ceiling(noClust/gridR)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    plot.new()
    par(mfrow=c(gridR,gridC))
    par(mar = c(2,2,0,0))
    
    clusters <- list()
    centroids <- list()
    for (i in 1:noClust) {
        # identify clusters
        index <- which(kmeansObj$cluster==i)
        iCluster <- x[index,]
        clusters[[length(clusters)+1]] <- iCluster
        # plot cluster
        matplot(t(iCluster), type="l", col=col[1:(length(col)-1)], 
            lty=lty[1:(length(lty)-1)], xaxt="n", ...)
        positions <- seq(1, ncol(x),length = no.ticks)
        t <- seq(x.coord[1], x.coord[ncol(x)], length=no.ticks)
        axis(side=1, at=positions, labels = t)
        legend("topleft", legend=paste("Cluster",i))
        # plot centroid
        matlines(kmeansObj$centers[i,], col=col[length(col)], 
            lty=lty[length(lty)], ...)
        centroids[[length(centroids)+1]] <- kmeansObj$centers[i,]
    }
    return(invisible(list(clusters=clusters, centroids=centroids)))
}



