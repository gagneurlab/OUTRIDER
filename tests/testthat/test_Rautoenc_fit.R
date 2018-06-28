context("Autoencoder: ")

set.seed(42)

## simulation of data
n = 100 # samples
p = 20 # genes
q = 6 # latent space dimension
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
randLoss <- loss(w, k, x, s, xbar, theta)

## init with pca on centered log counts
library(pcaMethods)
x  <- t(t(x0)- xbar)
pca <- pca(x, nPcs = q) ## could also use ppca
pc  <- loadings(pca)
w_guess <- as.vector(pc)
## indeed much better guess
pcaLoss <- loss(c(w_guess, rep(0,p)), k, x, s, xbar, theta)

## init with pca on centered log counts
## we'd need the gradient for this one to go fast
## It might be a good idea to restrict to Wd = t(We).
## I feel this is too much over-fitted otherwise

fit <- optim(c(w_guess, rep(0,p)), loss, k=k, s=s, x=x, xbar=xbar, theta=theta, 
        method="L-BFGS-B", control = list(maxiter=300))
fitLoss <- loss(fit$par, k, x, s, xbar, theta)
fit2 <- optim(c(w_guess, rep(0,p)), loss, gr = lossGrad, k=k, x=x, s=s, 
        xbar=xbar, theta=25, method="L-BFGS-B", control = list(maxiter=300))
fitGradLoss <- loss(fit2$par, k, x, s, xbar, theta)

test_that("Dimensions match.", {
    expect_equal(length(w_guess), p*q)
    expect_equal(length(fit$par), p*(q+1))
    expect_equal(length(lossGrad(fit$par, k, x, s, xbar, theta)), p*(q+1))
})

test_that("Losses match (fitLoss ~= fitGradLoss)",{
    expect_equal(fitLoss, fitGradLoss, tol=0.1)
})

test_that("randLoss > pcaLoss > fitGradLoss",{ 
    expect(randLoss >= pcaLoss)
    expect(pcaLoss >= fitGradLoss)
})

