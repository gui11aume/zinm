source("../Rpackage/R/nm.R")
dyn.load("../Rpackage/src/Rzinm.so")

set.seed(123)

lambda = rgamma(n=1000, shape=.8, rate=1)

x = rpois(n=1000, lambda=100*lambda)
y = rpois(n=1000, lambda=35*lambda)


target = c(1.18090772, 0.01092251, 0.73333581, 0.25574168)
stopifnot(max(abs(nm(cbind(x,y)) - target)) < 1e-6)
