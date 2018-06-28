context("Test loss and gradient functions: ")

# w : weights
# k : read counts matrix (sample x genes)
# s : size factor (scalar for now)
# xbar: mean gene expression levels

## simulation of data
n = 100 # samples
p = 20 # genes
q = 5 # latent space dimension
s=rnorm(n,mean=1, sd = 0.1)
theta = 25

h_true <- matrix(rnorm(n*q), nrow=n, ncol=q)
Wd_true <- matrix(rnorm(p*(q+1)), nrow=p, ncol=q+1)
y_true <- Wd_true%*%t(cbind(rep(1,n),h_true))
k <- apply(
    y_true,
    1,
    function(yi)
        rnbinom(n, mu=s*exp(yi), size=theta)
)

# compute x
x0 <- log((1+k)/s)
xbar <- colMeans(x0)
x <- t(t(x0)-xbar)

# ## value for random weight init
w <- c(rnorm(p*q, sd=1/p*q), rep(0,p))

# expect a real positiv number.
test_that("loss is real number", {
    l <- loss(w, k, x, s, xbar, theta)
    expect_that(l ,is.numeric)
    expect_that(length(l) == 1)
})

# numeric Gradient
numericLossGrad <- function(fn, epsilon, w,...){
    grad <- numeric(length(w))
    for(i in seq_along(w)){
        eps <- integer(length(w))
        eps[i] <- epsilon
        grad[i] <- (fn(w + eps, ...) - fn(w -eps, ...))/(2*epsilon)
    }
    return(grad)
}

# testing the gradient
#plot(numericLossGrad(loss, 1E-8, w, k=k, x=x, s=s, xbar=xbar, theta=25),
#     lossGrad(w, k, x, s, xbar, theta));abline(0,1)

test_that("Analytic gradient ~= numeric gradient", {
    for(i in 1:5){
        w <- c(rnorm(p*q, sd=1/p*q), rep(0,p))
        expect_equal(lossGrad(w, k, x, s, xbar, theta), 
            matrix(numericLossGrad(
                    loss, 1E-8, w, k=k, x=x, s=s, xbar=xbar, theta=theta
                        ), ncol=q+1), tol=0.0001)
    }
})