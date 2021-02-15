elbowRule <- function (x, k, method="Ward", nstart=1, B = 10, J = 2,
    x.coord = NULL, OSF = 1, vect = NULL, intercept = TRUE, degPolyn = 3,
    degFr = 4:20, knots = NULL, plot = FALSE, ...){

distortion <- function (dataMatrix, labels) {
    s <- c()  
    for(i in 1:max(labels)){
        pos <- which(labels %in% i)
        if (length(pos) > 1) { a <- apply(dataMatrix[pos,],2,"mean")
            } else { a <- dataMatrix[pos,] }
        if (length(pos)<=1) {s <- rbind(s, rep(0,length(a)))
            } else { s <- rbind(s, sweep(dataMatrix[which(labels %in% i),],2,a,"-")) }
        }
    return(sum(s^2))
}

rel.distortion <- c()
for (DF in degFr){
    km <- fabrik(x=x, k=k, method=method, nstart=nstart, 
    B = B, J = 2, x.coord=x.coord, OSF = OSF, vect = vect, intercept = intercept, 
    degPolyn = degPolyn, degFr = DF, knots = NULL, ...)$km
    rel.distortion <- c(rel.distortion, distortion(x,km$cluster))
    }

optimal.df <- degFr[which.min(rel.distortion)]

if (plot) {
    plot(degFr, rel.distortion, type="b", xlab = "Degrees of Freedom", ylab = "Distortion", pch=21)
    points(optimal.df, rel.distortion[which.min(rel.distortion)], pch=22, cex=1.5, col=2)
    }
return(list(df=degFr, tot.withinss=rel.distortion, optimal=optimal.df))
}

