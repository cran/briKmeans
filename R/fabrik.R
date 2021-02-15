fabrik <- function (x, k, method="Ward", nstart=1, B = 10, J = 2,
    x.coord = NULL, OSF = 1, vect = NULL, intercept = TRUE, degPolyn = 3,
    degFr = 5, knots = NULL, ...)
{
    # initial checks
    x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
    
    if (is.null(x.coord)) {x.coord <- 1:d}
    else {if (length(x.coord) != d) {
        stop("Wrong number of coordinates")
        }
    }
    
    # implemented methods
    impl.methods <- c("Ward", "pam", "kmeans")
    if (!(method %in% impl.methods)) {
        stop("Select an appropriate clustering method.")
    }
  
    # fit a curve to each multivariate data point
    if(is.null(vect)){
        xx.coord <- seq( min(x.coord), max(x.coord),
            length.out = length(x.coord) * OSF )
    }else{
        xx.coord <- vect
    }
    xx <- matrix(ncol=length(xx.coord), nrow=n)

    if (requireNamespace("splines", quietly = TRUE)) {
        iBasis <- splines::bs(x=x.coord, degree = degPolyn, df=degFr, intercept = intercept,
            knots = knots)
        fBasis <- splines::bs(x=xx.coord, degree = degPolyn, df = degFr, intercept = intercept,
            knots = knots)
    } else {
        stop("'splines' package required!")
    }
    
    # fit appropriate spline to x to construct xx
    for (i in 1:n){
        y <- x[i,]
        myCoeffs <- lm(y ~ iBasis+0) # + 0 (no intercept)
    
        # formats the coefficient data into matrix form
        myCoeffsMatrix <- matrix(myCoeffs$coefficients)
        myCoeffsMatrix[is.na(myCoeffsMatrix)] <- 0 # replaces NA's by zeros
        # formats the basis into matrix form
        fBasisMatrix <- matrix(fBasis,nrow = nrow(fBasis), ncol = ncol(fBasis))
    
        # get the smoothed data
        xx[i,] <- fBasisMatrix %*% myCoeffsMatrix
    }

    ## run kmeans by initialising centers with 'deepest'
    final <- brik(x = xx, k = k, method = method, nstart=nstart,
        B = B, J = J, ...)
    
    return(list(seeds = final$seeds, km = final$km ))
}
