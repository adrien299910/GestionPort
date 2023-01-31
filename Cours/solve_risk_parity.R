# Risk parity portfolio

library(NlcOptim)

sigma <- c(.1, .2, .3)
rho <- matrix(c(1,.8, .7,.8, 1, .6, .7, .6, 1), nrow = 3)
Sigma <- diag(sigma) %*% rho %*% diag(sigma)

obj <- function(w) {
  w <- matrix(w, ncol=1)
  Sw <- Sigma %*% w
  MRC <- w * Sw
  sum(diff(MRC)^2)*1000
}

Aeq <- matrix(rep(1,3), nrow=1)
Beq <- 1
lb <- rep(0,3)
ub <- rep(1,3)

w.0 <- rep(1/3,3)

res <- solnl(X=w.0, objfun=obj, Aeq=Aeq, Beq=Beq, lb=lb, ub=ub)
