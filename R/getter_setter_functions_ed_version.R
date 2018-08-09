
getx <- function(ods){
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    b <- colMeans(x0)
    x <- t(t(x0) - b)
    
    return(x)
}

setD <- function(ods, D){
    if(!is.matrix(D)){
        D <- matrix(D, nrow=nrow(ods))
    }
    metadata(ods)[['D']] <- D
    return(ods)
}

getD <- function(ods){
    return(metadata(ods)[['D']])
}

setE <- function(ods, E){
    if(!is.matrix(E)){
        E <- matrix(E, nrow=nrow(ods))
    }
    metadata(ods)[['E']] <- E
    return(ods)
}

getE <- function(ods){
    return(metadata(ods)[['E']])
}

setb <- function(ods, b){
    mcols(ods)[['b']] <- b
    return(ods)
}

getb <- function(ods){
    return(mcols(ods)[['b']])
}

getw <- function(ods){
    return(c(as.vector(getE(ods)), as.vector(getD(ods)), getb(ods)))
}


replaceCountsByZscore <- function(ods, zScore=15){
    lk <- log2(counts(ods) + 1)
    z <- (lk - rowMeans(lk)) / rowSds(lk)
    k2replace <- abs(z) > zScore
    
    print(table(k2replace))
    counts(ods)[k2replace] <- matrix(rowMeans(counts(ods)), )
}
