## -----------------------------------------------------------------------------
set.seed(22035)
pareto <- function(a,b) {
  # Generate random numbers with cdf F(x)
  u <- runif(10000)
  x <- b*(1-u)^(-1/a)
  
  # Draw the histogram of random numbers generated
  hist(x, prob = TRUE, main = paste('Pareto(',a,',',b,')'))
  
  # Draw the density function f(x)
  y <- seq(0, max(x), 0.1)
  lines(y, a*b^a/(y^(a+1)))
}

pareto(2, 2)

## -----------------------------------------------------------------------------
set.seed(22035)
beta <- function(a,b) {
  # Calculate constant c
  x0 <- (a-1)/(a+b-2)
  c <- x0^(a-1)*(1-x0)^(b-1)  # constant in pdf can be ignored
  
  # Generate random numbers with pdf f(x)
  n <- 10000
  k <- 0
  y <- numeric(n)
  while (k < n) {
    u <- runif(1)
    x <- runif(1) # random variate from g(x)
    if (x^(a-1)*(1-x)^(b-1) / c > u) {
      # accept x
      k <- k + 1
      y[k] <- x
    }
  }
  
  # Draw the histogram of random numbers generated
  hist(y, prob = TRUE, main = paste('Beta(',a,',',b,')'), xlab = "x")
  
  # Draw the density function f(x)
  z <- seq(0, 1, 0.01)
  lines(z, z^(a-1)*(1-z)^(b-1)*gamma(a+b)/gamma(a)/gamma(b))
}

beta(3, 2)

## -----------------------------------------------------------------------------
set.seed(22035)
expgamma <- function(r, beta) {
  # Generate random numbers from the mixture
  n <- 1000
  x <- rgamma(n, r, beta)
  y <- rexp(n, x)
  return(y)
}

r <- 4; beta <- 2
rnd = expgamma(r, beta)

## -----------------------------------------------------------------------------
# Draw the histogram of random numbers generated
hist(rnd, prob = TRUE, main = paste('Pareto(',r,',',beta,')'), xlab = "y")

# Draw the density function f(y)
y <- seq(0, max(rnd), 0.01)
lines(y, r*beta^r/(beta+y)^(r+1))

## -----------------------------------------------------------------------------
set.seed(22035)
# This part is copied from bb
quick_sort <- function(x){
  num <- length(x)
  if(num==0||num==1){return(x)
  }else{
    a <- x[1]
    y <- x[-1]
    lower <- y[y<a]
    upper <- y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}#form a loop
}


test<-sample(1:1e4)
system.time(quick_sort(test))[1]
test <- quick_sort(test)
# show the result of fast sort algorithm
test[1:10]
test[9991:10000]

## -----------------------------------------------------------------------------
set.seed(22035)
n <- c(1e4, 2e4, 4e4, 6e4, 8e4)
computation_time <- function(n){
  t <- numeric(100)
  set.seed(22035)
  for(i in 1:100){
    test <- sample(1:n)
    t[i] <- system.time(quick_sort(test))[1]
  }
  t_mean <- mean(t)
  return(t_mean)
}


an <- c(computation_time(n[1]),computation_time(n[2]),computation_time(n[3]),
       computation_time(n[4]),computation_time(n[5]))
an

## -----------------------------------------------------------------------------
tn <- n*log(n)
mylm <- lm(an~tn)
x <- seq(0,1e6,length.out=100)
b <- coefficients(mylm)
plot(tn, an, main="Regression line")
lines(x, b[1]+b[2]*x, col="red")

## -----------------------------------------------------------------------------
set.seed(22035)
m <- 1e4
U <- runif(m)
theta1 <- mean(exp(U))                # simple MC estimator
theta2 <- mean((exp(U)+exp(1-U))/2)   # antithetic variables estimator
var1 <- var(exp(U))                   # sample variance of simple MC
var2 <- var((exp(U)+exp(1-U))/2)      # sample variance of antithetic variables
theta1
theta2
100*(var1-var2)/var1      # empirical estimator of percent reduction of variance

## -----------------------------------------------------------------------------
set.seed(22035)
m <- 1e4
x <- rnorm(m)

g <- function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
}
f1 <- function(x) dnorm(x)

theta.hat1 <- mean(g(x)/f1(x))
var1 <- var(g(x)/f1(x))
cbind(theta.hat1, var1)

## -----------------------------------------------------------------------------
set.seed(22035)
y <- rgamma(m,shape = 3,rate = 1)
f2 <- function(x) dgamma(x, shape = 3, rate = 1)

theta.hat2 <- mean(g(y)/f2(y))
var2 <- var(g(y)/f2(y))
cbind(theta.hat2, var2)

## -----------------------------------------------------------------------------
d <- seq(1, 5, 0.05)
gs <- c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)),
        expression(f[1](x)==e^{-x^2/2}/sqrt(2*pi)),
        expression(f[2](x)==x^2*e^{-x}/2))
par(mfrow=c(1,2))
#figure (a)
plot(d, g(d), type = "l", ylab = "", ylim = c(0,0.5),
     lwd = 2,col=1,main='(A)')
lines(d, f1(d), lty = 2, lwd = 2,col=2)
lines(d, f2(d), lty = 3, lwd = 2,col=3)
legend("topright", legend = gs, lty = 1:3,
       lwd = 2, inset = 0.02,col=1:3)
#figure (b)
plot(d, g(d)/f1(d), type = "l", ylab = "", ylim = c(0,3),
     lwd = 2, lty = 2, col = 2, main = "(B)")
lines(d, g(d)/f2(d), lty = 3, lwd = 2, col = 3)
legend("topright", legend = gs[-1], lty = 2:3, lwd = 2,
       inset = 0.02, col = 2:3)

## -----------------------------------------------------------------------------
a <- numeric(5)
a[1] <- -log(0.8+exp(-1)/5)
a[5] <- 1
for(i in 2:4){
  a[i] <- -log(exp(-a[i-1])-0.2+exp(-1)/5)
}
a

## -----------------------------------------------------------------------------
set.seed(22035)
M <- 1e4
U <- runif(M)

g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

# density function on each subinterval
f <- function(x){
  5*exp(-x)/(1-exp(-1))
}

# inverse of distribution functions
F1_inv <- function(x){
  -log(1-(1-exp(-1))*x/5)
}
F2_inv <- function(x){
  -log(exp(-a[1])-(1-exp(-1))*x/5)
}
F3_inv <- function(x){
  -log(exp(-a[2])-(1-exp(-1))*x/5)
}
F4_inv <- function(x){
  -log(exp(-a[3])-(1-exp(-1))*x/5)
}
F5_inv <- function(x){
  -log(exp(-a[4])-(1-exp(-1))*x/5)
}

# samples generated by inverse transform method
x <- matrix(0, nrow = 2000, ncol = 5)
x[,1] <- F1_inv(U[1:2000])
x[,2] <- F2_inv(U[2001:4000])
x[,3] <- F3_inv(U[4001:6000])
x[,4] <- F4_inv(U[6001:8000])
x[,5] <- F5_inv(U[8001:10000])

# estimator of mean and variance on each subinterval
theta.hat <- numeric(5)
sigma2.hat <- numeric(5)
for(i in 1:5){
  theta.hat[i] <- mean(g(x[,i])/f(x[,i]))
  sigma2.hat[i] <- var(g(x[,i])/f(x[,i]))
}

# show the result
theta <- sum(theta.hat)
sigma2 <- sum(sigma2.hat)
se <- sqrt(sigma2)
cbind(theta, sigma2,se)

## -----------------------------------------------------------------------------
# sample generation function
sample_gen <- function(n, mu=0, sigma=1){
  x <- rlnorm(n=n, meanlog = mu, sdlog = sigma)
  return(x)
}

# data analysis function (constuct a confidence interval with level alpha)
CI <- function(x, alpha=0.05){
  n <- length(x)
  y <- log(x)
  mu.hat <- mean(y)
  sigma2.hat <- var(y)
  lower <- mu.hat+qt(alpha/2,df=n-1)*sqrt(sigma2.hat/n)
  upper <- mu.hat+qt(1-alpha/2,df=n-1)*sqrt(sigma2.hat/n)
  return(c("lower.bound"=lower,"upper.bound"=upper))
}

## -----------------------------------------------------------------------------
set.seed(22035)
m <- 1e4
lower <- upper <- numeric(m)

for(i in 1:m){
  Sample <- sample_gen(n=10, mu=0, sigma=1)
  lower[i] <- CI(x=Sample)[1]
  upper[i] <- CI(x=Sample)[2]
}

CP <- mean((lower<0)&(upper>0))
cat("CP =",CP)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# The functions of "Count Five" test is copied from the book
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

F.test <- function(x, y, alpha=0.05){
  S1 <- var(x)
  S2 <- var(y)
  m <- length(x)
  n <- length(y)
  f <- S2/S1
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(f>qf(1-alpha/2,df1 = n-1,df2 = m-1)||
                           f<qf(alpha/2,df1 = n-1,df2 = m-1)))
}

## -----------------------------------------------------------------------------
power_count5test <- function(m, n1, n2, sigma1, sigma2){
  mean(replicate(m, expr={
    x <- rnorm(n1, 0, sigma1)
    y <- rnorm(n2, 0, sigma2)
    count5test(x, y)
  }))
}

power_F.test <- function(m, n1, n2, sigma1, sigma2){
  mean(replicate(m, expr = {
    x <- rnorm(n1, 0, sigma1)
    y <- rnorm(n2, 0, sigma2)
    F.test(x, y, alpha = 0.055)
  }))
}

## -----------------------------------------------------------------------------
set.seed(22035)
m <- 1e4
# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
result1 <- numeric(3)
result2 <- numeric(3)
n <- c(20,100,1000)

for(i in 1:3){
  result1[i] <- power_count5test(m, n1=n[i], n2=n[i], sigma1, sigma2)
  result2[i] <- power_F.test(m, n1=n[i], n2=n[i], sigma1, sigma2)
}


knitr::kable(data.frame("size"=c(20,100,200),"count five test"=result1,
                          "F test"=result2),align="c")

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
rm(list = ls())

library(boot)

Sample1 <- function(x){
  samp <- aircondit[x]
  samp
}

Rate1 <- function(samp, i) {
  rat <- 1/mean(as.matrix(samp[i, ]))
  rat
}

Result1 <- function(samp,func,Rr){
  bo <- boot(samp, statistic = func, R = Rr)
  print(bo)
}

set.seed(22035)
samp <- Sample1(1)
resu <- Result1(samp,Rate1,2000)

detach(package:boot)

rm(list = ls())

## -----------------------------------------------------------------------------
rm(list = ls())

library(boot)

Sample2 <- function(x){
  samp <- aircondit[x]
  samp
}

Meant2 <- function(x, i) {
  mea <- mean(as.matrix(x[i, ]))
  mea
}

Result2 <- function(samp,func,Rr){
  bo <- boot(samp, statistic = func, R = Rr)
  re <- boot.ci(bo, type = c("norm", "perc", "basic", "bca"))
  print(bo)
  print(re)
  hist(bo$t, prob = TRUE, main = " ")
  points(bo$t0, 0, cex = 2, pch = 16)
  bo
}

set.seed(22035)
samp <- Sample2(1)
resu <- Result2(samp,Meant2,2000)

detach(package:boot)

rm(list = ls())


## -----------------------------------------------------------------------------
rm(list = ls())

skewness <- function(x,i) {
  #computes the sample skewness coeff.
  x_bar <- mean(x[i])
  x_bar
}

Sample3 <- function(n, mea, sd){
  samp <- rnorm(n, mea, sd)
  samp
}

Analysis3 <- function(m, func, Rr, n, mea, sd){
  library(boot)
  nornorm <- matrix(0, m, 2)
  norbasi <- matrix(0, m, 2)
  norperc <- matrix(0, m, 2)
  for (i in 1:m) {
    Samp <- Sample3(n, mea, sd)
    Skew <- boot(Samp, statistic = func, R=Rr)
    Nor <- boot.ci(Skew, type=c("norm","basic","perc"))
    nornorm[i,] <- Nor$norm[2:3]
    norbasi[i,] <- Nor$basic[4:5]
    norperc[i,] <- Nor$percent[4:5]
  }
  #Calculate the coverage probability of a normal distribution
  norm <- mean(nornorm[,1] <= s & nornorm[,2] >= s)
  basi <- mean(norbasi[,1] <= s & norbasi[,2] >= s)
  perc <- mean(norperc[,1] <= s & norperc[,2] >= s)
  #Calculate the probability of the left side of the normal distribution
  normleft <- mean(nornorm[,1] >= s )
  basileft <- mean(norbasi[,1] >= s )
  percleft <- mean(norperc[,1] >= s )
  #Calculate the right side probability of a normal distribution
  normright <- mean(nornorm[,2] <= s )
  basiright <- mean(norbasi[,2] <= s )
  percright <- mean(norperc[,2] <= s )
  analyresu <- c(norm, basi, perc, normleft, basileft, percleft, normright, basiright, percright)
  analyresu
}

Result3 <- function(sd, analyresu){
  dnam <- paste("N ( 0 ,", as.character(sd^2),")",seq="")
  Distribution <- c(dnam)
  Type <- c("basic", "norm", "perc")
  Left <- analyresu[4:6]
  Right <- analyresu[7:9]
  P.coverage <- analyresu[1:3]
  result <- data.frame(Distribution, Type, Left, Right, P.coverage)
  result
}

s <- 0
n <- 20
m <- 1000
R <- 1000

mea <- 0
sd <- 3 

# We can set n, m, R, mea, sd any way we want.

set.seed(22035)
library(boot)

Analyresu <- Analysis3(m, skewness, R, n, mea, sd)
Resu <- Result3(sd, Analyresu)

knitr::kable (Resu, align="c")

rm(list = ls())

## -----------------------------------------------------------------------------
library(bootstrap)
x <- as.matrix(scor)
n <- nrow(x)
lambda <- eigen(cov(x))$values
theta.hat <- max(lambda/sum(lambda)) # original estimate
theta.jack <- numeric(n)
for(i in 1:n){
  lambda.jack <- eigen(cov(x[-i, ]))$values
  theta.jack[i] <- max(lambda.jack/sum(lambda.jack)) # jackknife estimate
}
# the jackknife estimates of bias
bias.jack <- (n-1)*(mean(theta.jack) - theta.hat) 
# the jackknife estimates of standard error
se.jack <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2)) 
c(est=theta.hat, bias=bias.jack, se=se.jack)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic)  
N <- n*(n-1)/2 # all possible combination
e1 <- e2 <- e3 <- e4 <- numeric(n)
h<-1
# leave-two-out cross validation
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    k<-c(i,j)
    y <- magnetic[-k]
    x <- chemical[-k]
    # Model 1: Linear
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2]*chemical[k]
    e1[h] <- sum((magnetic[k] - yhat1)^2)
    # Model 2：Quadratic
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2]*chemical[k] +J2$coef[3]*chemical[k]^2
    e2[h] <- sum((magnetic[k] - yhat2)^2)
    # Model 3: Exponential
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2]*chemical[k]
    yhat3 <- exp(logyhat3)
    e3[h] <- sum((magnetic[k] - yhat3)^2)
    # Model 4: Log-Log
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4 <- exp(logyhat4)
    e4[h] <- sum((magnetic[k] - yhat4)^2)
    h<-h+1
  }
}
# the average squared prediction error by leave-two-out cross validation
c(Linear=sum(e1), Quadratic=sum(e2), Exponential=sum(e3), LogLog=sum(e4))/(2*N)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# define a function to conduct permutation test
spear.perm <- function(x, y, R=1e3){
  rho0<-cor.test(x, y, method = "spearman")$estimate # the original test estimate
  n<-length(y)
  reps <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:n)
    reps[i] <- cor.test(x, y[k], method = "spearman")$estimate
    }
  p <- mean(c(rho0, reps) >= rho0) # p-value of the permutation test
  return(p) # return p-value
}

## -----------------------------------------------------------------------------
library(MASS)
# generate data from multi-normal distribution
mu <- c(0, 0); sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
n<-30
set.seed(22035)
xy <- mvrnorm(n, mu, sigma)
x<-xy[, 1]; y<-xy[, 2]
# the p-value reported by cor.test on the same samples
(p0<-cor.test(x, y, method = "spearman")$p.value)
# the achieved signiﬁcance level of the permutation test
(p.perm<- spear.perm(x,y))
# compare the two p-value
round(c(p0=p0,perm=p.perm),4)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22035)

rl.metropolis <- function(sigma, x0, N) {
  # sigma: sd of proposal distribution N(xt,sigma^2)
  # x0: initial value
  # N: length of chain
  
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0  # to calculate acceptance rate
  for (t in 2:N) {
    y <- rnorm(1, x[t-1], sigma)
    if (u[t] <= exp(abs(x[t-1]) - abs(y))) { x[t] <- y; k <- k + 1 }
    else { x[t] <- x[t-1] }
  }
  return(list(mc = x, acc.prob = k / N))
}

N <- 10000
b <- 1000
k <- 4
sigma <- c(0.5, 1, 4, 16)
x0 <- c(-5, -2, 2, 5)
X <- matrix(nrow = k, ncol = N)
acc.prob <- numeric(k)
for (i in 1:k) {
  rl <- rl.metropolis(sigma[i], x0[i], N)
  X[i, ] <- rl$mc
  acc.prob[i] <- rl$acc.prob
}
acc.prob

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
for (i in 1:k) {
  plot(X[i,], type = "l", xlab = bquote(sigma == .(sigma[i])),
       ylab = "X", ylim = range(X[i,]))
}

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
x <- seq(-6, 6, 0.01)
fx <- exp(-abs(x)) / 2
for (i in 1:k) {
  hist(X[i, -(1:b)], breaks = "Scott", freq = FALSE, main = "",
       xlab = bquote(sigma == .(sigma[i])), xlim = c(-6, 6), ylim = c(0, 0.5),)
  lines(x, fx, col = 2, lty = 2)
}

## -----------------------------------------------------------------------------
z <- rexp(100, 1)
z <- c(-rev(z), z) # generate laplace random numbers
p <- c(0.05, seq(0.1, 0.9, 0.1), 0.95)
Q <- quantile(z, p)
mc <- X[, -(1:b)]
Qmc <- apply(mc, 1, function(x) quantile(x, p))
QQ <- data.frame(round(cbind(Q, Qmc), 3))
names(QQ) <- c('True', 'sigma=0.5', 'sigma=1', 'sigma=4', 'sigma=16')
knitr::kable(QQ)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# ergodic mean plot
phi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(phi)) {
  phi[i,] <- phi[i,] / (1:ncol(phi))
}
for (i in 1:k) {
  if (i == 1) {
    plot((b+1):N, phi[i, (b+1):N], ylim = c(-0.5, 0.5),
         type = "l", xlab = 'Index', ylab = bquote(phi))
  } else { lines(phi[i, (b+1):N], col = i) }
}

# plot of R_hat
rhat <- rep(0, N)
for (j in (b+1):N) {
  rhat[j] <- Gelman.Rubin(phi[, 1:j])
}
plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
# clear memory and set seed
#rm(list = ls())
set.seed(22035)

rbn.metropolis <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XY <- rbn.metropolis(mu, sigma, rho, mu, N)
X <- XY$X[-(1:b)]; Y <- XY$Y[-(1:b)]
plot(X, Y, xlab = bquote(X[t]), ylab = bquote(Y[t]),
     main = "", cex = 0.5, ylim = range(Y))
cov(cbind(X, Y))

## ----fig.height=4, fig.width=9------------------------------------------------
k <- 4
x0 <- matrix(c(2,2,-2,-2,4,-4,-4,4), nrow = 2, ncol = k)
Xmc <- Ymc <- XYmc <- matrix(0, nrow = k, ncol = N)
for (i in 1:k) {
  XY <- rbn.metropolis(mu, sigma, rho, x0[,i], N)
  Xmc[i,] <- XY$X; Ymc[i,] <- XY$Y
  XYmc[i,] <- Xmc[i,] * Ymc[i,]
}

# ergodic mean plot
cal_phi <- function(X) {
  phi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(phi)) {
    phi[i,] <- phi[i,] / (1:ncol(phi))
  }
  return(phi)
}
phiX <- cal_phi(Xmc)
phiY <- cal_phi(Ymc)
phiXY <- cal_phi(XYmc)

plot_erg_mean <- function(phi, rg) {
  for (i in 1:k) {
    if (i == 1) {
      plot((b+1):N, phi[i, (b+1):N], type = "l", ylim = rg,
           xlab = "Index", ylab = bquote(phi))
    }
    else { lines(phi[i, (b+1):N], col = i) }
  }
}
par(mfrow = c(1, 3))
plot_erg_mean(phiX, rg = c(-0.5, 0.5))
plot_erg_mean(phiY, rg = c(-0.5, 0.5))
plot_erg_mean(phiXY, rg = c(0.7, 1.1))

## ----fig.height=4, fig.width=9------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# plot of R_hat
plot_R_hat <- function(phi) {
  rhat <- rep(0, N)
  for (j in (b+1):N) {
    rhat[j] <- Gelman.Rubin(phi[, 1:j])
  }
  plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R", ylim = c(1, 1.25))
  abline(h = 1.2, lty = 2)
}
par(mfrow = c(1, 3))
plot_R_hat(phiX)
plot_R_hat(phiY)
plot_R_hat(phiXY)

## ----comment = ''-------------------------------------------------------------
lm.fit <- lm(Y ~ X)
summary(lm.fit)

## ----fig.height=5, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
e <- lm.fit$residuals
qx <- seq(-2, 2, 0.01)
hist(e, breaks = "Scott", freq = FALSE, main = "", xlim = c(-2, 2), ylim = c(0, 1))
lines(qx, dnorm(qx, 0, sqrt(0.19)), col = 2, lwd = 1.5)
qqnorm(e)
qqline(e, col = 2, lwd = 2, lty = 2)

## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22034)

rl.metropolis <- function(sigma, x0, N) {
  # sigma: sd of proposal distribution N(xt,sigma^2)
  # x0: initial value
  # N: length of chain
  
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0  # to calculate acceptance rate
  for (t in 2:N) {
    y <- rnorm(1, x[t-1], sigma)
    if (u[t] <= exp(abs(x[t-1]) - abs(y))) { x[t] <- y; k <- k + 1 }
    else { x[t] <- x[t-1] }
  }
  return(list(mc = x, acc.prob = k / N))
}

N <- 10000
b <- 1000
k <- 4
sigma <- c(0.5, 1, 4, 16)
x0 <- c(-5, -2, 2, 5)
X <- matrix(nrow = k, ncol = N)
acc.prob <- numeric(k)
for (i in 1:k) {
  rl <- rl.metropolis(sigma[i], x0[i], N)
  X[i, ] <- rl$mc
  acc.prob[i] <- rl$acc.prob
}
acc.prob

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
for (i in 1:k) {
  plot(X[i,], type = "l", xlab = bquote(sigma == .(sigma[i])),
       ylab = "X", ylim = range(X[i,]))
}

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
x <- seq(-6, 6, 0.01)
fx <- exp(-abs(x)) / 2
for (i in 1:k) {
  hist(X[i, -(1:b)], breaks = "Scott", freq = FALSE, main = "",
       xlab = bquote(sigma == .(sigma[i])), xlim = c(-6, 6), ylim = c(0, 0.5),)
  lines(x, fx, col = 2, lty = 2)
}

## -----------------------------------------------------------------------------
z <- rexp(100, 1)
z <- c(-rev(z), z) # generate laplace random numbers
p <- c(0.05, seq(0.1, 0.9, 0.1), 0.95)
Q <- quantile(z, p)
mc <- X[, -(1:b)]
Qmc <- apply(mc, 1, function(x) quantile(x, p))
QQ <- data.frame(round(cbind(Q, Qmc), 3))
names(QQ) <- c('True', 'sigma=0.5', 'sigma=1', 'sigma=4', 'sigma=16')
knitr::kable(QQ)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# ergodic mean plot
phi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(phi)) {
  phi[i,] <- phi[i,] / (1:ncol(phi))
}
for (i in 1:k) {
  if (i == 1) {
    plot((b+1):N, phi[i, (b+1):N], ylim = c(-0.5, 0.5),
         type = "l", xlab = 'Index', ylab = bquote(phi))
  } else { lines(phi[i, (b+1):N], col = i) }
}

# plot of R_hat
rhat <- rep(0, N)
for (j in (b+1):N) {
  rhat[j] <- Gelman.Rubin(phi[, 1:j])
}
plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
# clear memory and set seed
#rm(list = ls())
set.seed(22034)

rbn.metropolis <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XY <- rbn.metropolis(mu, sigma, rho, mu, N)
X <- XY$X[-(1:b)]; Y <- XY$Y[-(1:b)]
plot(X, Y, xlab = bquote(X[t]), ylab = bquote(Y[t]),
     main = "", cex = 0.5, ylim = range(Y))
cov(cbind(X, Y))

## ----fig.height=4, fig.width=9------------------------------------------------
k <- 4
x0 <- matrix(c(2,2,-2,-2,4,-4,-4,4), nrow = 2, ncol = k)
Xmc <- Ymc <- XYmc <- matrix(0, nrow = k, ncol = N)
for (i in 1:k) {
  XY <- rbn.metropolis(mu, sigma, rho, x0[,i], N)
  Xmc[i,] <- XY$X; Ymc[i,] <- XY$Y
  XYmc[i,] <- Xmc[i,] * Ymc[i,]
}

# ergodic mean plot
cal_phi <- function(X) {
  phi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(phi)) {
    phi[i,] <- phi[i,] / (1:ncol(phi))
  }
  return(phi)
}
phiX <- cal_phi(Xmc)
phiY <- cal_phi(Ymc)
phiXY <- cal_phi(XYmc)

plot_erg_mean <- function(phi, rg) {
  for (i in 1:k) {
    if (i == 1) {
      plot((b+1):N, phi[i, (b+1):N], type = "l", ylim = rg,
           xlab = "Index", ylab = bquote(phi))
    }
    else { lines(phi[i, (b+1):N], col = i) }
  }
}
par(mfrow = c(1, 3))
plot_erg_mean(phiX, rg = c(-0.5, 0.5))
plot_erg_mean(phiY, rg = c(-0.5, 0.5))
plot_erg_mean(phiXY, rg = c(0.7, 1.1))

## ----fig.height=4, fig.width=9------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# plot of R_hat
plot_R_hat <- function(phi) {
  rhat <- rep(0, N)
  for (j in (b+1):N) {
    rhat[j] <- Gelman.Rubin(phi[, 1:j])
  }
  plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R", ylim = c(1, 1.25))
  abline(h = 1.2, lty = 2)
}
par(mfrow = c(1, 3))
plot_R_hat(phiX)
plot_R_hat(phiY)
plot_R_hat(phiXY)

## ----comment = ''-------------------------------------------------------------
lm.fit <- lm(Y ~ X)
summary(lm.fit)

## ----fig.height=5, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
e <- lm.fit$residuals
qx <- seq(-2, 2, 0.01)
hist(e, breaks = "Scott", freq = FALSE, main = "", xlim = c(-2, 2), ylim = c(0, 1))
lines(qx, dnorm(qx, 0, sqrt(0.19)), col = 2, lwd = 1.5)
qqnorm(e)
qqline(e, col = 2, lwd = 2, lty = 2)

## -----------------------------------------------------------------------------
#data generation
model <- function(alpha,beta,gamma,N){
  e_M <- rnorm(N,0,1)
  e_Y <- rnorm(N,0,1)
  x <- rnorm(N,0,1)
  m <- alpha*x+e_M
  y <- beta*m+gamma*x+e_Y
  return(data.frame(x,m,y))
}

## -----------------------------------------------------------------------------
t <- function(x,m,y){              #Test statistics
  data <- data.frame(x,m,y)
  f1 <- lm(m~x,data=data)
  f2 <- lm(y~m+x,data=data)
  s_alpha <- summary(f1)$coefficients[2,2]
  s_beta <- summary(f2)$coefficients[2,2]
  t <- f1$coef[2]*f2$coef[2]/sqrt(f1$coef[2]^2*s_alpha^2+f2$coef[2]^2*s_beta^2)
  return(as.numeric(t))
}

## -----------------------------------------------------------------------------
N <- 50
R <- 499
K <- 1:N
test1 <- function(data){
  reps <- numeric(R)
for (i in 1:R) {
  k <- sample(K,size = N,replace = FALSE)
  x <- data[k,1]
  reps[i] <- t(x,data$m,data$y)
}
p_hat <- (sum(abs(reps)>abs(t(data$x,data$m,data$y)))+1)/(1+R)
return(p_hat<0.05)
}
test2 <- function(data){
  reps <- numeric(R)
for (i in 1:R) {
  k <- sample(K,size = N,replace = FALSE)
  y <- data[k,3]
  reps[i] <- t(data$x,data$m,y)
}
  p_hat <- (sum(abs(reps)>abs(t(data$x,data$m,data$y)))+1)/(1+R)
return(p_hat<0.05)
}
test3 <- function(data){
  reps <- numeric(R)
for (i in 1:R) {
  k <- sample(K,size = N,replace = FALSE)
  m <- data[k,2]
  reps[i] <- t(data$x,m,data$y)
}
  p_hat <- (sum(abs(reps)>abs(t(data$x,data$m,data$y)))+1)/(1+R)
return(p_hat<0.05)
}

## -----------------------------------------------------------------------------
stimulation <- function(alpha,beta,gamma=1,N=50){
type_1_error <- matrix(0,nrow = 3,ncol = 200)
for (i in 1:200) {
  data <- model(alpha,beta,gamma,N)
  type_1_error[1,i] <- test1(data)
  type_1_error[2,i] <- test2(data)
  type_1_error[3,i] <- test3(data)
}
return(c(mean(type_1_error[1,]),mean(type_1_error[2,]),mean(type_1_error[3,])))
}

## -----------------------------------------------------------------------------
s1 <- stimulation(alpha=0,beta=0)
s2 <- stimulation(alpha=0,beta=1)
s3 <- stimulation(alpha=1,beta=0)

## -----------------------------------------------------------------------------
result <- data.frame(model1=s1,model2=s2,model3=s3)
row.names(result) <- c("test1","test2","test3")
result

## -----------------------------------------------------------------------------
f <- function(f0,N=10^6,b1=0,b2=1,b3=-1){
x1 <- rpois(N,1)
x2 <- rexp(N,1)
x3 <- rbeta(N,1,0.5)
g <- function(alpha){
tmp <- exp(alpha+b1*x1+b2*x2+b3*x3)
p <- 1/(1+tmp)
mean(p) - f0
}
solution <- uniroot(g,c(-100,100))
alpha <- solution$root
return(alpha)
}

## -----------------------------------------------------------------------------
f0 <- c(0.1,0.01,0.001,0.0001)
alpha <- numeric(4)
for (i in 1:4) {
  alpha[i] <- f(f0[i])
}

## -----------------------------------------------------------------------------
plot(-log(f0),alpha)

## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------

#'*lower bound*
u <- c(11,8,27,13,16,0,23,10,24,2)
#'*upper bound*
v <- c(12,9,28,14,17,1,24,11,25,3)


#'*Observed data likelihood*
o.likelihood <- function(lambda){
  sum((v*exp(-lambda*v)-u*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v)))
}



solution <- uniroot(o.likelihood,interval = c(0.05,0.1),extendInt = "yes")

k <- round(unlist(solution),5)[1:3]


MLE <- k[1]

#'*EM algorithm*

lambda.old <- 0.0000000001
N <- 1e5

tol <- .Machine$double.eps

options(digits=10)
for(j in 1:N) {
  
  lambda <- length(u)/(sum((u*exp(-lambda.old*u)-v*exp(-lambda.old*v))/(exp(-lambda.old*u)-exp(-lambda.old*v)))+length(u)/lambda.old)

  if ((abs(lambda - lambda.old)/lambda.old) < tol) break
  lambda.old <- lambda
}




## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------
x <- list(c(1, 2), list(3, 4))
str(x)
str(unlist(x))
str(as.vector(x))

## -----------------------------------------------------------------------------
dim(c(1, 2, 3)) # atomic vector
dim(list(1, 2, list(3))) # list

## -----------------------------------------------------------------------------
x <- matrix(1:6, nrow = 2, ncol = 3)
c(is.matrix(x), is.array(x))

y <- array(1:12, c(2, 3, 2))
c(is.matrix(y), is.array(y))

## -----------------------------------------------------------------------------
x <- data.frame(V1 = c(1, 2, 3),
                V2 = c("a", "b", "c"),
                V3 = c(TRUE, FALSE, FALSE),
                row.names = c("X1", "X2", "X3"))
x
attributes(x)
dim(x)

## -----------------------------------------------------------------------------
x <- data.frame(
  V1 = c(1L, 2L),
  V2 = c(FALSE, TRUE),
  V3 = c("a", "b")
)
as.matrix(x)

y <- data.frame(
  V1 = c(1.5, 2.0),
  V2 = c(FALSE, TRUE)
)
as.matrix(y)

## -----------------------------------------------------------------------------
# 0 rows
x <- data.frame(V1 = numeric())
c(nrow(x), ncol(x))

# 0 columns
y <- data.frame(row.names = "X1")
c(nrow(y), ncol(y))

# 0 rows, 0 columns
z <- data.frame()
c(nrow(z), ncol(z))

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
x <- data.frame(x1 = c(1.5, 2.5, 3.5, 4.5), x2 = rnorm(4, 4, 4))
str(x)

## -----------------------------------------------------------------------------
res1 <- data.frame(lapply(x, scale01))
res1

## -----------------------------------------------------------------------------
# add a non-numeric column
x$x3 = c(rep("A", 2), rep("B", 2))

res2 <- data.frame(lapply(x, function(x) if (is.numeric(x)) scale01(x) else x))
res2

## -----------------------------------------------------------------------------
rm(list = ls())
x <- data.frame(x1 = c(0.6, 1.3, 7.6, 2.4), x2 = rnorm(4, 2, 2))
str(x)

## -----------------------------------------------------------------------------
res1 <- vapply(x, sd, 1)
res1

## -----------------------------------------------------------------------------
# add a non-numeric column
x$x3 = c(rep("A", 2), rep("B", 2))

res2 <- vapply(x[vapply(x, is.numeric, TRUE)], sd, 1)
res2

## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22034)

gibbsR <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

## ----fig.height=4, fig.width=8------------------------------------------------
# import package
library(Rcpp)
library(StatComp22035)

# generate chains
N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XYR <- gibbsR(mu, sigma, rho, mu, N)
XR <- XYR$X[-(1:b)]; YR <- XYR$Y[-(1:b)]
XYC <- gibbsC(mu, sigma, rho, mu, N)
XC <- XYC[-(1:b), 1]; YC <- XYC[-(1:b), 2]

par(mfrow = c(1, 2))
qqplot(XR, XC, plot.it = TRUE)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
qqplot(YR, YC, plot.it = TRUE)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)

## ----fig.height=4, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
plot(XR, YR, cex = 0.5)
plot(XC, YC, cex = 0.5)

## -----------------------------------------------------------------------------
# import package
library(microbenchmark)
ts <- microbenchmark(gibbsR = gibbsR(mu, sigma, rho, mu, N),
                     gibbsC = gibbsC(mu, sigma, rho, mu, N))
summary(ts)[, c(1, 3, 5, 6)]

