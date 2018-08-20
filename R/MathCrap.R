n <- 5
A <- matrix(rnorm(n*n), n, n)

# Minor and cofactor
minor <- function(A, i, j) det( A[-i,-j] )
cofactor <- function(A, i, j) (-1)^(i+j) * minor(A,i,j)

# With a loop
adjoint1 <- function(A) {
    n <- nrow(A)
    B <- matrix(NA, n, n)
    for( i in 1:n )
        for( j in 1:n )
            B[j,i] <- cofactor(A, i, j)
    B
}

# With `outer`
adjoint <- function(A) {
    n <- nrow(A)
    t(outer(1:n, 1:n, Vectorize(
        function(i,j) cofactor(A,i,j)
    )))
}

# Check the result: these should be equal
det(A) * diag(nrow(A))
A %*% adjoint1(A)
A %*% adjoint2(A)

if(FALSE){
microbenchmark::microbenchmark(A %*% adjoint1(A),
                               A %*% adjoint2(A))

n <- 50
A <- matrix(rnorm(n*n), n, n)

microbenchmark::microbenchmark(A %*% adjoint1(A),
                               A %*% adjoint2(A))
}