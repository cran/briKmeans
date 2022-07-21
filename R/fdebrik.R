fdebrik <- function (x, k, method="Ward", nstart=1, B = 10, J = 2,
    x.coord = NULL, functionalDist="d0.pearson", 
    OSF = 1, vect = NULL, intercept = TRUE, degPolyn = 3,
    degFr = 5, knots = NULL, ...)
{
    # initial checks
    x <- as.matrix(x) # original data
    d <- ncol(x)      # dimension of original data
    n <- nrow(x)
    
    if (is.null(x.coord)) {x.coord <- 1:d} # define abscisa
    else {if (length(x.coord) != d) {
        stop("Wrong number of coordinates")
        }
    }
  
    # implemented methods
    impl.methods <- c("Ward", "pam", "kmeans")
    if (!(method %in% impl.methods)) {
        stop("Select an appropriate clustering method.")
    }
    
    # check functional distance
    if (functionalDist=="d0.pearson") {derivative <- FALSE}
    else {
        if (functionalDist=="d1.pearson") {derivative <- TRUE}
        else {stop("Not implemented functional distance.")}
    }
    
    # define oversampled abscisa
    if(is.null(vect)){
        xx.coord <- seq( min(x.coord), max(x.coord),
            length.out = length(x.coord) * OSF )
    }else{
        xx.coord <- vect
    }
    D <- length(xx.coord)
    
    ## fit a curve to each multivariate data point
    xx <- matrix(ncol=length(xx.coord), nrow=n)
    # construct the basis
    if (requireNamespace("splines2", quietly = TRUE)) {
        iBasis <- splines2::bSpline(x=x.coord, degree = degPolyn, 
            df=degFr, intercept = intercept, knots = knots)
        fBasis <- splines2::bSpline(x=xx.coord, degree = degPolyn, 
            df = degFr, intercept = intercept, knots = knots)
    } else {
        stop("'splines2' package required!")
    }

    if (derivative) {
        fBasisDer <- deriv(fBasis)
        yy0 <- yy1 <- matrix(ncol=D, nrow=n)
        # fit a b-spline to each multivariate data
        myCoeffsMatrix <- matrix(ncol=ncol(fBasis), nrow=n)
        for (i in 1:n){
            y <- x[i,]
            myCoeffs <- lm(y ~ iBasis+0)
            myCoeffsMatrix <- matrix(myCoeffs$coefficients)
            myCoeffsMatrix[is.na(myCoeffsMatrix)] <- 0
            # evaluate each fitted function and its derivative 
            yy0[i,] <- fBasis %*% myCoeffsMatrix
            yy1[i,] <- fBasisDer %*% myCoeffsMatrix
            } 
        yy0yy1 <- cbind(yy0,yy1)
        
        # define the function to compute bootstrap centers
        kma.centers1 <- function(xx0xx1, w) {
            d <- ncol(xx0xx1)/2
            Y <-xx0xx1[w, 1:d]
            YDer <- xx0xx1[w, (d+1):(2*d)]
            obtainedCenters <- matrix(nrow=k, ncol=D)
            centers <- kma(xx.coord, y0=Y, y1=YDer, n.clust=k,
                n.out=D, warping.method="NOalignment",
                similarity.method=functionalDist)$y0.centers.final
            obtainedCenters[1:nrow(centers),] <- centers
            return(obtainedCenters)
            }                  
        }
    else {
        fBasisDer <- NULL
        yy0 <- matrix(ncol=D, nrow=n)
        # fit a b-spline to each multivariate data
        myCoeffsMatrix <- matrix(ncol=ncol(fBasis), nrow=n)
        for (i in 1:n){
            y <- x[i,]
            myCoeffs <- lm(y ~ iBasis+0)
            myCoeffsMatrix <- matrix(myCoeffs$coefficients)
            myCoeffsMatrix[is.na(myCoeffsMatrix)] <- 0
            # evaluate each fitted function to get the input for kma
            yy0[i,] <- fBasis %*% myCoeffsMatrix
            yy1 = NULL
            }
        yy0yy1 <- yy0
        # define the function to compute bootstrap centers
        kma.centers0 <- function(xx0, w) {
            d <- ncol(xx0)
            obtainedCenters <- matrix(nrow=k, ncol=d)
            centers <- kma(x=xx.coord, y0=xx0[w,], n.clust=k,
                n.out=d, warping.method="NOalignment",
                similarity.method= functionalDist)$y0.centers.final
            obtainedCenters[1:nrow(centers),] <- centers
            return(obtainedCenters)
            }
        }
    
    # bootstrap step
    if (requireNamespace("boot", quietly = TRUE)) {
        if (derivative) {results <- boot::boot(yy0yy1,kma.centers1,B)}
        else {results <- boot::boot(yy0yy1, kma.centers0, B)}
        centers <- as.vector(results$t)
    }
    else {
        stop("'boot' package required!")
    }
    
    new.centers <- matrix(centers, ncol=D)
    index <- apply(new.centers,1,function(x) all(is.na(x)))
    new.centers <- new.centers[!index,]
    
    if (method == "Ward") {
        distances <- dist(new.centers)
        partition <- cutree(tree = hclust(distances, method = "ward.D2"),
        k = k)
    }
    
    if (method == "pam") {
        if (requireNamespace("cluster", quietly = TRUE)) {
            partition <- cluster::pam(new.centers, k)$cluster
        }
        else {
            stop("'cluster' package required!")
        }
    }
    if (method == "kmeans") {
        partition <- kmeans(new.centers, k, nstart = nstart)$cluster
    }
    
    deepest <- matrix(0, k, D)
    for (i in 1:k) {
        iCluster <- new.centers[partition == i, ]
        if (length(iCluster) == D)
        iCluster <- matrix(iCluster, 1, ncol = D)
        if (requireNamespace("depthTools", quietly = TRUE)) {
            depth <- depthTools::MBD(iCluster, plotting = FALSE)
        }
        else {
            stop("'depthTools' package required!")
        }
        deepest[i, ] <- iCluster[which(depth$ordering == 1),
        ]
    }
    
    final <- kmeans(x=yy0, centers=deepest, nstart = 1)
    return(list(seeds = deepest, km = final))
}
