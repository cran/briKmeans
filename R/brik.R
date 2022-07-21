brik <- function (x, k, method="Ward", nstart=1,
    B=10, J = 2, ...)
{
    # implemented methods
    impl.methods <- c("Ward", "pam", "kmeans")
    if (!(method %in% impl.methods)) {
        stop("Select an appropriate clustering method.")
    }
    
    # format data as matrix
    x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
    
    #  function to retrieve centroids computed by kmeans
    kmean.centers <- function(x, w) {
        centers <- kmeans(x[w, ], k)$centers
        return(centers)
    }

    # obtain the centroid for each bootstrap replicate
    if (requireNamespace("boot", quietly = TRUE)) {
        centers <- boot::boot(x, kmean.centers, B)$t
    } else {
        stop("'boot' package required!")
    }
    
    # store the centroids correctly arranged
    new.centers <- matrix(0, B * k, d)
    for (i in 1:d) {
        new.centers[, i] <- centers[, ((i - 1) * k + 1):(i * 
            k)]
    }

    # cluster the centroids with 'method'
    if (method=="Ward"){
        distances <- dist(new.centers)
        partition <- cutree(tree = hclust(distances,
            method='ward.D2'), k=k)
    }

    if (method=="PAM"){
        if (requireNamespace("cluster", quietly = TRUE)) {
            partition <- cluster::pam(new.centers, k)$cluster
        } else {
            stop("'cluster' package required!")
        }
    }
    if (method=="kmeans"){
        partition <- kmeans(new.centers, k, nstart=nstart,
            ...)$cluster
    }

     # find the deepest point within each cluster
    deepest <- matrix(0, k, d)
    for (i in 1:k) {
        iCluster <- new.centers[partition == i,]
        if (length(iCluster)==d) iCluster <- matrix(iCluster,1,ncol=d)
        if (requireNamespace("depthTools", quietly = TRUE)) {
            depth <- depthTools::MBD(iCluster, plotting=FALSE)
        } else {
            stop("'depthTools' package required!")
        }
        
        deepest[i, ] <- iCluster[which(depth$ordering==1),]
    }
    
    ## run kmeans by initialising centers with 'deepest'
    final <- kmeans(x, deepest, nstart=1, ...)

    ## output the quality indices
    return(list(seeds = deepest, km = final))
}
