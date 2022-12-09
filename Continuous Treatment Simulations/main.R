library(doParallel)
library(doRNG)

source("calm.R")

Test.Rules1 <- function(n = 1e4, Betas) {
  x1 <- runif(n)
  x2 <- runif(n)
  a <- (-(cbind(1,x1,x2) %*% t(Betas[, colnames(Betas) != "A2"])) / 
          t(2*replicate(n, Betas[, colnames(Betas) == "A2"])))
  a <- pmin(a, 5)
  a <- pmax(a, 0)
  p <- ncol(a)
  y <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    v <- (2 + x1 + x2 + a[,j]/5)
    unexposed.model <- 6 + 3*x1 + 2*x1^2 + 3*x2^2
    contrast.model <- 2*a[,j] + a[,j]*x1 + a[,j]*x2 - a[,j]^2
    y[,j] <- unexposed.model + contrast.model + rnorm(n, 0, sd = sqrt(v))
  }
  apply(y, 2, mean)
}

do.One <- function(n, SL.lib = c("SL.mean", "SL.glm.interaction", "SL.nnet", 
                                 "SL.rpart", "SL.earth")) {
  x1 <- runif(n)
  x2 <- runif(n)
  a <- 5*exp(-rnorm(n, plogis(-0.5 + 0.2*x1 + 0.2*x2))^2) * 
    (1-rbinom(n, 1, plogis(-0.5 + 0.2*x1 + 0.2*x2)))
  v <- (2 + x1 + x2 + a/5)
  unexposed.model <- 6 + 3*x1 + 2*x1^2 + 3*x2^2
  contrast.model <- 2*a + a*x1 + a*x2 - a^2
  y <- unexposed.model + contrast.model + rnorm(n, 0, sd = sqrt(v))
  
  Betas <- calm(baseline = data.frame(X1 = x1, X2 = x2), 
                exposure = data.frame(A1 = a, A2 = a^2), 
                outcome = y,
                contrast.models = list( ~ 1 + . , ~ 1 ),
                outcome.learners = c(SL.lib), 
                exposure.learners = c(SL.lib), 
                variance.learners = c(SL.lib),
                zero.learners = c(SL.lib))
  
  Values <- Test.Rules1(1e4, Betas)
  names(Values) <- rownames(Betas)
  return(Values)
}


Test.Rules2 <- function(n = 1e4, Betas) {
  x1 <- runif(n)
  x2 <- runif(n)
  a <- (-(cbind(1,x1,x2) %*% t(Betas[, colnames(Betas) != "A2"])) / 
          t(2*replicate(n, Betas[, colnames(Betas) == "A2"])))
  a <- pmin(a, 5)
  a <- pmax(a, 0)
  p <- ncol(a)
  y <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    v <- (2 + x1 + x2 + a[,j]/5)
    unexposed.model <- 6 + 3*x1 + 2*x1^2 + 3*x2^2
    contrast.model <- 2*a[,j] + a[,j]*x1^2 + a[,j]*x2^2 - a[,j]^2
    y[,j] <- unexposed.model + contrast.model + rnorm(n, 0, sd = sqrt(v))
  }
  apply(y, 2, mean)
}

do.Two <- function(n, SL.lib = c("SL.mean", "SL.glm.interaction", "SL.nnet", 
                                 "SL.rpart", "SL.earth")) {
  x1 <- runif(n)
  x2 <- runif(n)
  a <- 5*exp(-rnorm(n, plogis(-0.5 + 0.2*x1 + 0.2*x2))^2) * 
    (1-rbinom(n, 1, plogis(-0.5 + 0.2*x1 + 0.2*x2)))
  v <- (2 + x1 + x2 + a/5)
  unexposed.model <- 6 + 3*x1 + 2*x1^2 + 3*x2^2
  contrast.model <- exp(2*a + a*x1^2 + a*x2^2 - a^2)
  y <- unexposed.model + contrast.model + rnorm(n, 0, sd = sqrt(v))
  
  Betas <- calm(baseline = data.frame(X1 = x1, X2 = x2), 
                exposure = data.frame(A1 = a, A2 = a^2), 
                outcome = y,
                contrast.models = list( ~ 1 + . , ~ 1 ),
                outcome.learners = c(SL.lib), 
                exposure.learners = c(SL.lib), 
                variance.learners = c(SL.lib),
                zero.learners = c(SL.lib))
  
  Values <- Test.Rules2(1e4, Betas)
  names(Values) <- rownames(Betas)
  return(Values)
}


num_to_parallelize <- detectCores(all.tests = FALSE, logical = TRUE)-1 #find out how many cores are on the machine
cl <- makeCluster(num_to_parallelize)
registerDoParallel(cl)

registerDoRNG(912)
start.time <- Sys.time()
result <- foreach (i = 1:200, .combine = rbind, 
                   .packages = c('SuperLearner')) %dopar% {do.One(1000, "SL.ranger")}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

result |> boxplot()
result |> summary()

registerDoRNG(912)
start.time <- Sys.time()
result2 <- foreach (i = 1:200, .combine = rbind, 
                   .packages = c('SuperLearner')) %dopar% {do.Two(1000, "SL.ranger")}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

result2 |> boxplot()
result2 |> summary()
