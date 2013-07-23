### Power analysis for IES Assessment Studies
### Compute in two ways: simple t-test and multiple regression.
### Output should be graphics showing three lines (three different effect sizes).
### One graphic will be power vs. sample size
### Another graphic will be sample size vs. effect size

library(pwr)

### Set some constants.
### Sample size will be on horizontal axis, will range from 10 to 2000 in plots.
ssize <- seq(50,1500,50)
### Set alpha level to 0.05. This means there is 5% chance of comitting a Type I error (rejecting the null when in fact the null is true)
alpha <- 0.05
### Set three effect sizes: small, medium, and large
effsize <- c(0.2, 0.5, 0.8)
### Let power be revealed on the vertical axis. Power, or beta, is 1 minus the probability of making a Type II error (failing to reject the null when in fact it is false)

## T-tests
## Power vs sample size
plot(ssize,pwr.t.test(n=ssize,d=effsize[1],sig.level=alpha)$power, type="l", xlab="Sample Size", ylab="Power",axes=F,main=expression(paste("Power Analysis for T-Test, ", alpha==0.05, " (Two-Tailed)")))
lines(ssize, pwr.t.test(n=ssize,d=effsize[2],sig.level=alpha)$power, lty=2, col="red")
lines(ssize, pwr.t.test(n=ssize,d=effsize[3],sig.level=alpha)$power, lty=3, col="blue")
axis(side=1, at=seq(50,1500,50))
axis(side=2, at=seq(0,1,.2))
abline(v=0, h=seq(0,1,.2), lty=2, col="grey89")
abline(v=seq(50,1500,50), h=0, lty=2, col="grey89")
legend("bottomright", title="Effect Size", as.character(effsize), lty=1:3, col=c("black", "red", "blue"))

### Sample size vs. Effect Size
effsize <- seq(.05, .95, .05)
pwr <- c(0.7, 0.8, 0.9)
plot(effsize, sapply(effsize, function(x) pwr.t.test(d=x,sig.level=alpha, power=pwr[1])$n), type="l", xlab="Effect Size", ylab="Sample Size",axes=F)
title(expression(paste("Power Analysis for T-Test, ", alpha==0.05, " (Two-Tailed)")))
lines(effsize, sapply(effsize, function(x) pwr.t.test(d=x,sig.level=alpha, power=pwr[2])$n), lty=2, col="red")
lines(effsize, sapply(effsize, function(x) pwr.t.test(d=x,sig.level=alpha, power=pwr[3])$n), lty=3, col="blue")
axis(side=1, at=effsize)
axis(2, at=seq(0, 5000, 500))
abline(v=0, h=seq(0, 5000, 500), lty=2, col="grey89")
abline(v=effsize, h=0, lty=2, col="grey89")
legend("topright", title="Power", as.character(pwr), lty=1:3, col=c("black", "red", "blue"))


## Multiple Regression
### Need more assumptions for multiple regression (u)
ssize <- seq(50,1500,50)
alpha <- 0.05
effsize <- c(0.02, 0.15, 0.35) #different values here for multiple regression
## u is equal to the number of continuous predictors plus the number of dummy variables minus one.
## Here assume that we have 8 continuous vars and 8 dummy vars, so 16 total -1 = 15
u <- 15
## v will be our denominator degrees of freedom. so u+v+1 is our sample size
## Power vs sample size
plot(ssize,pwr.f2.test(u=u,v=ssize, f2=effsize[1],sig.level=alpha)$power, type="l", xlab="Sample Size", ylab="Power",axes=F,main=expression(paste("Power Analysis for Regression, ", alpha==0.05, ", 15 Controls")))
lines(ssize, pwr.f2.test(u=u,v=ssize, f2=effsize[2],sig.level=alpha)$power, lty=2, col="red")
lines(ssize, pwr.f2.test(u=u,v=ssize, f2=effsize[3],sig.level=alpha)$power, lty=3, col="blue")
axis(side=1, at=seq(50,1500,50))
axis(side=2, at=seq(0,1,.2))
abline(v=0, h=seq(0,1,.2), lty=2, col="grey89")
abline(v=seq(50,1500,50), h=0, lty=2, col="grey89")
legend("bottomright", title="Effect Size", as.character(effsize), lty=1:3, col=c("black", "red", "blue"))


### Sample size vs. Effect Size
effsize <- seq(.02, .5, .02)
pwr <- c(0.7, 0.8, 0.9)
u <-15
plot(effsize, sapply(effsize, function(x) pwr.f2.test(u=u, f2=x,sig.level=alpha, power=pwr[1])$v)+u+1, type="l", xlab="Effect Size", ylab="Sample Size",axes=F, main=expression(paste("Power Analysis for Regression, ", alpha==0.05, ", 15 Controls")))
lines(effsize, sapply(effsize, function(x) pwr.f2.test(u=u, f2=x,sig.level=alpha, power=pwr[2])$v)+u+1, lty=2, col="red")
lines(effsize, sapply(effsize, function(x) pwr.f2.test(u=u, f2=x,sig.level=alpha, power=pwr[3])$v)+u+1, lty=3, col="blue")
axis(side=1, at=effsize)
axis(2, at=seq(0, 800, 100))
abline(v=0, h=seq(0, 800, 100), lty=2, col="grey89")
abline(v=effsize, h=0, lty=2, col="grey89")
legend("topright", title="Power", as.character(pwr), lty=1:3, col=c("black", "red", "blue"))

