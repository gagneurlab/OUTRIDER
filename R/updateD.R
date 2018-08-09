
samples <- 80
q<- 2

#Hidden Space is a sample time q matrix. 
H <- matrix(c(rep(c(-1,1), each=samples/2), 
              rep(c(-2,2), samples/2)), ncol=2)

H
D_true <- rnorm(q)
y_true <- H %*% D_true + 3
mu_true <- 0.01 + exp(y_true)

k <- rnbinom(length(mu_true), mu = mu_true, size=25)

#library(MASS)
#glm.nb(k~H)

init<-c(mean(log(k+1)),0,0)
#b <- 3
#d <- D_true

lossD <- function(d, k, H, s=1, theta){
    b <- d[1]
    d <- d[-1]
    y <- H %*% d + b
    mu <- 0.01 + exp(y)
    -mean(dnbinom(k, mu=mu, size=theta, log=TRUE))
}


gradD <- function(d, k, H, s=1, theta){
    b <- d[1]
    d <- d[-1]
    y <- c(H %*% d + b)
    mu <- 0.01 + s*exp(y)
    t1<- colMeans(c(k*s*exp(y)/mu) * H)
    
    kt <- (k + theta)*s*exp(y)/(s*exp(y)+theta)
    t2 <- colMeans(c(kt*s*exp(y)/mu) * H) 
    
    dd<- t2-t1
    db<- mean(kt - k*s*exp(y)/mu)
    c(db, dd)
}

fit <- optim(init, fn=lossD, gr=gradD, k=k, H=H, s=1, theta=25, method='L-BFGS')
fit$par

D_true
lossD(k, init, H, 25)
lossD(k,c(3, D_true), H,  25)
