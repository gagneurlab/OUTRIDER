
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

## check that the simulated counts seem to make sense
# hist(log10(1+k))

# compute x
x0 <- log((1+k)/s)
xbar <- colMeans(x0)
x <- t(t(x0)-xbar)

# ## value for random weight init
w <- c(rnorm(p*q, sd=1/p*q), rep(0,p))
# expect a real positiv number.
loss(w, k, x, s, xbar, theta)

# numeric Gradient
numericLossGrad <- function(fn, epsilon, w,...){
    grad <- list()
    for(i in seq_along(w)){
        eps <- rep(0, length(w))
        eps[i] <- epsilon
        grad[i] <- (fn(w + eps, ...) - fn(w -eps, ...))/2*epsilon
    }
    return(grad)
}

# testing the gradient
plot(numericLossGrad(loss, 1E-5, w, k=k, x=x, s=s, xbar=xbar, theta=25),
     lossGrad(w, k, x, s, xbar, theta))

for(i in 1:5){
    w <- c(rnorm(p*q, sd=1/p*q), rep(0,p))
    print(cor(as.numeric(numericLossGrad(loss, 1E-3, w, k=k, x=x, s=s, xbar=xbar, theta=theta)),
     lossGrad(w, k, x, s, xbar, theta)) > 0.99)
}
