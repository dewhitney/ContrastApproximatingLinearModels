# DGP1 <- function(n=10L, p=3L){
#   X <- data.frame(matrix(runif(n*p,-1,1), ncol=p))
#   A <- rbinom(n,1,plogis(X[,1]))
#   Q0 <- 1 + 2*X[,1] + X[,2] + 0.5*X[,3]
#   CATE <- 0.442*(1-X[,1]-X[,2]-0*X[,3])
#   Q <- Q0 + A*CATE
#   e <- rnorm(n)
#   Y <- Q + e
#   out <- list(X=X, A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
#   return(out)
# }
# 
# DGP1a <- function(n=300L) DGP1(n, p=3L)
# DGP1b <- function(n=300L) DGP1(n, p=10L)
# 
# DGP2 <- function(n=300L){
#   Z <- matrix(rnorm(n*4L), ncol=4L)
#   A <- rbinom(n,1,plogis(Z %*% c(-1,.5,-.25,-.1)))
#   Q0 <- 2 + Z %*% c(2,1,1,1)
#   CATE <- 1/2 + Z %*% c(2,1,-2,-1)
#   Q <- Q0 + A*CATE
#   e <- rnorm(n)
#   X <- data.frame(X1 = exp(Z[,1]/2),
#                   X2 = 1 + Z[,2]/(1 + exp(Z[,1])),
#                   X3 = (Z[,1]*Z[,3]/25 + 0.6)^3,
#                   X4 = (Z[,2] + Z[,4] + 2)^2)
#   Y <- Q + e
#   out <- list(X=X, A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
#   return(out)
# }
# 
# DGP3 <- function(n=300L){
#   Z <- matrix(rnorm(n*4L), ncol=4L)
#   X <- data.frame(X1 = exp(Z[,1]/2),
#                   X2 = 1 + Z[,2]/(1 + exp(Z[,1])),
#                   X3 = (Z[,1]*Z[,3]/25 + 0.6)^3,
#                   X4 = (Z[,2] + Z[,4] + 2)^2)
#   A <- rbinom(n,1,plogis(Z %*% c(-1,.5,-.25,-.1)))
#   Q0 <- 2 + Z %*% c(2,1,1,1)
#   CATE <- 1/2 + as.matrix(X) %*% c(2,1,-2,-1)
#   Q <- Q0 + A*CATE
#   e <- rnorm(n)
#   Y <- Q + e
#   out <- list(X=X, A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
#   return(out)
# }
# 
# DGP4 <- function(n=300L){
#   Z <- matrix(rnorm(n*4L), ncol=4L)
#   X <- data.frame(X1 = exp(Z[,1]/2),
#                   X2 = 1 + Z[,2]/(1 + exp(Z[,1])),
#                   X3 = (Z[,1]*Z[,3]/25 + 0.6)^3,
#                   X4 = (Z[,2] + Z[,4] + 2)^2)
#   A <- rbinom(n,1,plogis(as.matrix(X) %*% c(-1,.5,-.25,-.1)))
#   Q0 <- 2 + Z %*% c(2,1,1,1)
#   CATE <- 1/2 + Z %*% c(2,1,-2,-1)
#   Q <- Q0 + A*CATE
#   e <- rnorm(n)
#   Y <- Q + e
#   out <- list(X=X, A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
#   return(out)
# }

DGP_proto <- function(n=300L, linear = c(F,F,F)){
  Z <- matrix(rnorm(n*4L), ncol=4L)
  X <- data.frame(X1 = exp(Z[,1]/2),
                  X2 = 1 + Z[,2]/(1 + exp(Z[,1])),
                  X3 = (Z[,1]*Z[,3]/25 + 0.6)^3,
                  X4 = (Z[,2] + Z[,4] + 2)^2)
  if (linear[3]) {
    A <- rbinom(n,1,plogis(as.matrix(X) %*% c(-1,.5,-.25,-.1)))
  } else {
    A <- rbinom(n,1,plogis(Z %*% c(-1,.5,-.25,-.1)))
  }
  if (linear[2]) {
    Q0 <- 2 + as.matrix(X) %*% c(2,1,1,1)
  } else {
    Q0 <- 2 + Z %*% c(2,1,1,1)
  }
  if (linear[1]) {
    CATE <- 1/2 + as.matrix(X) %*% c(2,1,-2,-1)
  } else {
    CATE <- 1/2 + Z %*% c(2,1,-2,-1)
  }
  Q <- Q0 + A*CATE
  e <- rnorm(n)
  Y <- Q + e
  out <- list(X=X, A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
  return(out)
}

DGP1 <- function(n) {
  DGP_proto(n, linear = c(T,T,T))
}

DGP2 <- function(n) {
  DGP_proto(n, linear = c(T,T,F))
}

DGP3 <- function(n) {
  DGP_proto(n, linear = c(T,F,T))
}

DGP4 <- function(n) {
  DGP_proto(n, linear = c(T,F,F))
}

DGP5 <- function(n) {
  DGP_proto(n, linear = c(F,T,T))
}

DGP6 <- function(n) {
  DGP_proto(n, linear = c(F,T,F))
}

DGP7 <- function(n) {
  DGP_proto(n, linear = c(F,F,T))
}

DGP8 <- function(n) {
  DGP_proto(n, linear = c(F,F,F))
}

# Linear models in observed variables!
# DGP9 <- function(n=300L){ 
#   Z <- matrix(rnorm(n*4L), ncol=4L)
#   A <- rbinom(n,1,plogis(Z %*% c(-1,.5,-.25,-.1)))
#   Q0 <- 2 + Z %*% c(2,1,1,1)
#   CATE <- 1/2 + Z %*% c(2,1,-2,-1)
#   Q <- Q0 + A*CATE
#   e <- rnorm(n)
#   Y <- Q + e
#   X <- data.frame(Z)
#   names(X) <- paste0("X", 1:4)
#   out <- list(X=X, A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
#   return(out)
# }

DGP9 <- function(n=10L, p=10L){
  X <- data.frame(matrix(rnorm(n*p,-1,1), ncol=p))
  A <- rbinom(n,1,plogis(X[,1]))
  Q0 <- apply(X^2, 1, sum) + apply(X, 1, sum)
  CATE <- X[,1]+X[,2]-0.1
  Q <- Q0 + A*CATE
  e <- rnorm(n)
  Y <- Q + e
  out <- list(X=data.frame(X,X^2), A=A, Y=Y, Q0=Q0, CATE=CATE, e=e, dopt=1L*(CATE >= 0))
  return(out)
}


