# Merton's model: option on index
# r: interest rate
# b: cost of carry: int. rate - div. yield

GBSPrice <- function(PutCall, S, K, T, r, b, sigma) {
  d1 <- (log(S/K) + (b+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)


  if(PutCall == 'c')
    px <- S*exp((b-r)*T)*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
  else
    px <- K*exp(-r*T)*pnorm(-d2) - S*exp((b-r)*T)*pnorm(-d1)

px
}

GBSVega <- function(PutCall, S, K, T, r, b, sigma) {
  d1 <- (log(S/K) + (b+sigma^2/2)*T)/(sigma*sqrt(T))
  S*exp((b-r)*T) * dnorm(d1)
}

ImpliedVolNewton <- function(p, TypeFlag, S, X, Time, r, b,
                             sigma=NULL, maxiter=500, tol=1.e-5) {

  if(is.null(sigma))
     s <- sqrt(2*abs(log(S/(X*exp((b*T)))))/T)
  else
    s <- sigma

  not_converged <- T
  i=1
  vega <- GBSVega(TypeFlag, S, X, Time, r, b, s)
  while(not_converged & (i<maxiter)) {
    err <- (p-GBSPrice(TypeFlag, S, X, Time, r, b, s))
    s <- s + err/vega
    # print(paste('i:', i, 's:', s))
    not_converged <- (abs(err/vega) > tol)
    i <- i+1
  }
s
}

BinaryPrice <- function(PutCall, S, K, T, r, b, sigma) {
  d1 <- (log(S/K) + (b+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)


  if(PutCall == 'c')
    px <- K*exp(-r*T)*pnorm(d2)
  else
    px <- K*exp(-r*T)*pnorm(-d2)

px
}

# 1 - Volatility interpolation
##############################

T <- 1
Spot <- 100
r <- 0
b <- 0
eps <- .001
sigma <- .3

# Benchmark data: (strike, volatility)
VolData <- list(c(80, .32), c(100, .30), c(120, .315))


# Define an array of pricing functions for
# the three benchmark instruments:

C1 <- function(vol=sigma, spot=Spot) GBSPrice(PutCall='c',
      S=spot, K=VolData[[1]][1], T=T, r=r, b=b, sigma=vol)

C2 <- function(vol=sigma, spot=Spot) GBSPrice(PutCall='c',
      S=spot, K=VolData[[2]][1], T=T, r=r, b=b, sigma=vol)

C3 <- function(vol=sigma, spot=Spot) GBSPrice(PutCall='c',
      S=spot, K=VolData[[3]][1], T=T, r=r, b=b, sigma=vol)

C <- c(C1, C2, C3)

# utility functions to compute the risk indicators,
# all by finite difference:

Vega <- function(f, vol, spot=Spot) (f(vol+eps, spot)-f(vol-eps, spot))/(2*eps)

Vanna <- function(f, vol, spot=Spot) {
  (Vega(f, vol, spot+1)-Vega(f, vol, spot-1))/2
}

Volga <- function(f, vol) {
    (Vega(f,vol+eps)-Vega(f,vol-eps))/(eps)
}

# Example: To compute the vega of 3 hedge options
# B.vega <- sapply(1:3, function(i) Vega(C[[i]], sigma))

# 1 - Compute vectors of vega, vanna, volga for
#     3 hedge instruments

B.vega <- sapply(1:3, function(i) Vega(C[[i]], sigma))
B.vanna <- sapply(1:3, function(i) Vanna(C[[i]], sigma))
B.volga <- sapply(1:3, function(i) Volga(C[[i]], sigma))

# 2 - Option with new strike

O <- function(vol=sigma, spot=Spot) GBSPrice('c', S=spot,
       K=Knew, T=T, r=r, b=b, sigma=vol)

Knew <- 90

# Its Black-Scholes price
O.BS <- O()

# 2 - risk indicators for new option
####################################

O.vega <- Vega(O, sigma)
O.vanna <- Vanna(O, sigma)
O.volga <- Volga(O, sigma)

# 3 - Benchmark costs: difference between price
#     at market vol and price at ATM vol
###############################################

B.cost <- sapply(1:3, function(i) C[[i]](VolData[[i]][2]) - C[[i]](sigma))

# 4 - calculation of price adjustment
#####################################

A <- t(matrix(c(B.vega, B.vanna, B.volga),  nrow=3))
x <- matrix(c(O.vega, O.vanna, O.volga), nrow=3)
w <- solve(A, x)
CF <- t(w) %*% matrix(B.cost, nrow=3)

# adjusted price
O.Price <- O.BS + CF

# implied volatility
O.iv <- ImpliedVolNewton(O.Price, 'c', Spot, Knew, T, r, b,
                         sigma=sigma)

# 5 - TODO
# Wrap the above logic in a function to interpolate/extrapolate

## VVVol(K) computes the implied vol at strike K

VVVol <- function(K) {

## this function computes the price of a call strike K,
## as a function of vol and spot

  O <- function(vol=sigma, spot=Spot) GBSPrice('c', S=spot,
       K=K, T=T, r=r, b=b, sigma=vol)

  # Its Black-Scholes price
  O.BS <- O()

  # risk indicators for new option
  O.vega <- Vega(O, sigma)
  O.vanna <- Vanna(O, sigma)
  O.volga <- Volga(O, sigma)

  # calculation of price adjustment
  A <- t(matrix(c(B.vega, B.vanna, B.volga),  nrow=3))
  x <- matrix(c(O.vega, O.vanna, O.volga), nrow=3)
  w <- solve(A, x)
  CF <- t(w) %*% matrix(B.cost, nrow=3)

  print(paste('CF:', CF, 'BS:', O.BS, 'K:', K))

  # implied volatility
  iv <- ImpliedVolNewton(O.BS+CF, 'c', Spot, K, T, r, b,
                         sigma=sigma)
  iv
}

# the vol curve from K=70 to K=130

v <- sapply(seq(70, 130, 2), VVVol)

        plot(seq(70, 130,2), v, type='l', lwd=3, xlab='Strike',
     ylab='Implied Volatility')
points(sapply(VolData, function(v) v[1]),
          sapply(VolData, function(v) v[2]), pch=19, col='red')

# 6 - Price a binary call
#########################

# Compare
# 1 - price with ATM vol
# 2 - price with interpolated vol
# 3 - price with VV adjustment
#################################

T <- 1
Spot <- 100
r <- 0
b <- 0

# ATM BS vol
sigma <- .30

# strike
Strike <- 110

# smile function
smile <- function(X) (-(0/20)*(X-Spot) + (1/300)*(X-Spot)^2)/100

SmileVol <- function(X) sigma + smile(X)

Bin <- function(vol=sigma, spot=Spot) BinaryPrice('c', S=spot,
       K=Strike, T=T, r=r, b=b, sigma=vol)

# Its Black-Scholes price
Bin.BS <- Bin()

# Price with interpolated vol
Bin.IntVol <- Bin(vol=SmileVol(Strike))

# price and risk indicators for benchmark instruments

B.cost <- sapply(1:3,
          function(i) C[[i]](SmileVol(VolData[[i]][1])) - C[[i]](sigma))

B.vega <- sapply(1:3, function(i) Vega(C[[i]], sigma))
B.vanna <- sapply(1:3, function(i) Vanna(C[[i]], sigma))
B.volga <- sapply(1:3, function(i) Volga(C[[i]], sigma))

# 2 - risk indicators for new option
####################################

Bin.vega <- Vega(Bin, sigma)
Bin.vanna <- Vanna(Bin, sigma)
Bin.volga <- Volga(Bin, sigma)

# 4 - calculation of price adjustment
#####################################

A <- t(matrix(c(B.vega, B.vanna, B.volga),  nrow=3))
x <- matrix(c(Bin.vega, Bin.vanna, Bin.volga), nrow=3)
w <- solve(A, x)
CF <- t(w) %*% matrix(B.cost, nrow=3)

# adjusted price
Bin.VVPrice <- Bin.BS + CF

print('Binary price')
print(paste('BS:', Bin.BS))
print(paste('Interpolated Vol:', Bin.IntVol))
print(paste('VV price:', Bin.VVPrice))

## wrap calculation in a function

BinaryVVPrice <- function(Strike) {

  Bin <- function(vol=sigma, spot=Spot) BinaryPrice('c', S=spot,
       K=Strike, T=T, r=r, b=b, sigma=vol)

  # Its Black-Scholes price
  Bin.BS <- Bin()

  # risk indicators for new option
  ####################################

  Bin.vega <- Vega(Bin, sigma)
  Bin.vanna <- Vanna(Bin, sigma)
  Bin.volga <- Volga(Bin, sigma)

  # calculation of price adjustment
  #####################################

  x <- matrix(c(Bin.vega, Bin.vanna, Bin.volga), nrow=3)
  w <- solve(A, x)
  CF <- t(w) %*% matrix(B.cost, nrow=3)

  # adjusted price
  c(Bin.BS, Bin.BS + CF)
}

## at Strike 65 and 66
tmp <- BinaryVVPrice(65)
BS.Price.65 <- tmp[1]
VV.Price.65 <- tmp[2]

tmp <- BinaryVVPrice(66)
BS.Price.66 <- tmp[1]
VV.Price.66 <- tmp[2]

# BS Delta and VV Delta: see the difference

print(paste('BS delta at 65', round(BS.Price.66-BS.Price.65,3)))
print(paste('VV delta at 65', round(VV.Price.66-VV.Price.65,3)))
