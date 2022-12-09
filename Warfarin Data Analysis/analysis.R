
library(SuperLearner)
library(MASS)
library(dplyr)
library(doParallel)
library(foreach)
library(doRNG)

source('clean_data.R')

# W-M weights based on quantile binned doses (Schulz & Moodie 2021, Supplement)
w.m.fcn<-function(treat.outcome,treat.mod,data,m){
  
  A<-model.response(model.frame(treat.outcome,data))
  
  Xalpha<-model.matrix(treat.mod[[1]],model.frame(treat.mod[[1]],data))
  
  a.binned<-ntile(A,m)
  a.cat<-as.factor(a.binned)
  
  if(ncol(Xalpha)==1){
    pom<-polr(a.cat~1)
  }
  else{
    pom<-polr(a.cat~Xalpha[,-1])
  }
  
  a.pom.prob<-pom$fitted.values[cbind(seq_along(a.binned), a.binned)]
  
  weight.vec<-(1-pom$fitted.values[,m])/(a.pom.prob)
  
  return(weight.vec)
  
}

# SL wrapper for glmnet to incorporate quadratic dose term and interactions
SL.glmnet.blip <- function(Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
                           nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + Dose*. + I(Dose^2)*., X)
    newX <- model.matrix(~-1 + Dose*. + I(Dose^2)*., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

# Pharmacogenetic dosing rule
do.one <- function(Obs, SL.lib = c("SL.mean", "SL.glm"), n.train = 800) {
  # sample size
  N <- dim(Obs)[1]
  
  # choose IDs for training (n = 800) and testing datasets
  train.id <- sample(1:N,n.train)
  test.id <- setdiff(c(1:N),train.id)
  
  # Extract outcome, dose, and covariate data frames
  Y <- -abs(Obs$INR - 2.5)
  A0 <- as.numeric(Obs$Dose == 35)
  A1 <- Obs$Dose
  A2 <- Obs$Dose^2
  X <- Obs[,setdiff(names(Obs), c("Dose","INR"))]
  AX <- Obs[,setdiff(names(Obs), c("INR"))]
  
  # Test/train set outcomes
  test.Y <- Y[test.id]
  train.Y <- Y[train.id]
  
  # Test/train set dose variables
  test.A0 <- A0[test.id]
  train.A0 <- A0[train.id]
  test.A1 <- A1[test.id]
  train.A1 <- A1[train.id]
  test.A2 <- A2[test.id]
  train.A2 <- A2[train.id]
  
  # Test/train set covariates
  test.X <- X[test.id,]
  train.X <- X[train.id,]
  test.AX <- AX[test.id,]
  train.AX <- AX[train.id,]
  
  # Nuisance regressions
  fit.Y.AX <- SuperLearner(train.Y, train.AX, SL.library = c("SL.glmnet.blip", SL.lib))
  fit.Y.0X <- SuperLearner(train.Y[train.A0==1], train.X[train.A0==1,],
                           newX = train.X, SL.library = SL.lib)
  fit.Y.X <- SuperLearner(train.Y, train.X, SL.library = SL.lib)
  fit.A0.X <- SuperLearner(train.A0, train.X, SL.library = SL.lib)
  fit.A1.X <- SuperLearner(train.A1, train.X, SL.library = SL.lib)
  fit.A2.X <- SuperLearner(train.A2, train.X, SL.library = SL.lib)
  
  # Robinson estimator
  Y.X <- train.Y - fit.Y.X$SL.predict[,1]
  A1.X <- train.A1 - fit.A1.X$SL.predict[,1]
  A2.X <- train.A2 - fit.A2.X$SL.predict[,1]
  A1.mat <- A1.X * cbind(1, train.X)
  names(A1.mat) <- paste0("A1.", names(A1.mat))
  A2.mat <- A2.X
  names(A2.mat) <- paste0("A2.", names(A2.mat))
  A.mat <- cbind(A1.mat, A2.mat)
  
  b.R <- coef(lm(Y.X ~ -1 + ., data = A.mat))
  a.R <- -0.5*(b.R[1]+as.matrix(test.X)%*%cbind(b.R[2:14]))/b.R[15]
  a.R <- pmin(a.R, 95)
  a.R <- pmax(a.R, 5.81)

  yhat.R <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.R,test.X))[1]$pred)
  yhat.R <- pmin(yhat.R, 0)
  yhat.R <- pmax(yhat.R, -1.2)
  v.R <- mean(yhat.R)
  
  # Blip-to-35 estimator
  Y.0 <- train.Y - fit.Y.0X$SL.predict[,1]
  A1.X <- train.A1 - 35
  A2.X <- train.A2 - 35^2
  A1.mat <- A1.X * cbind(1, train.X)
  names(A1.mat) <- paste0("A1.", names(A1.mat))
  A2.mat <- A2.X
  names(A2.mat) <- paste0("A2.", names(A2.mat))
  A.mat <- cbind(A1.mat, A2.mat)
  
  b.35 <- coef(lm(Y.0 ~ -1 + ., data = A.mat)) #plug-in
  
  wts <- train.A0 / fit.A0.X$SL.predict[,1]
  a1.pred <- fit.A1.X$SL.predict[,1]
  a2.pred <- fit.A2.X$SL.predict[,1]
  A.X.mat <- as.matrix(cbind(a1.pred, a1.pred*train.X, a2.pred))
  wtsY.0 <- cbind(wts*Y.0)
  one.step <- -c(chol2inv(chol(t(as.matrix(A.mat)) %*% as.matrix(A.mat))) %*% t(A.X.mat) %*% wtsY.0)
  b.35 <- b.35 + one.step #one-step
  
  a.35 <- -0.5*(b.35[1]+as.matrix(test.X)%*%cbind(b.35[2:14]))/b.35[15]
  a.35 <- pmin(a.35, 95)
  a.35 <- pmax(a.35, 5.81)

  yhat.35 <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.35,test.X))[1]$pred)
  yhat.35 <- pmin(yhat.35, 0)
  yhat.35 <- pmax(yhat.35, -1.2)
  v.35 <- mean(yhat.35)
  
  # Q-learning (linear regression)
  b.Q <- coef(lm(train.Y ~ Dose*. + I(Dose^2), data=train.AX))
  a.Q <- -0.5*(b.Q[2]+as.matrix(test.X)%*%cbind(b.Q[17:29]))/b.Q[16]
  a.Q <- pmin(a.Q, 95)
  a.Q <- pmax(a.Q, 5.81)
  
  yhat.Q <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.Q,test.X))[1]$pred)
  yhat.Q <- pmin(yhat.Q, 0)
  yhat.Q <- pmax(yhat.Q, -1.2)
  v.Q <- mean(yhat.Q)

  # DWOLS (IPW, binned quantiles w/ordered logistic regression)
  w.D <- w.m.fcn(train.A1~1, list(train.A1 ~ .), train.X, 20) # Weights
  b.D <- coef(lm(train.Y ~ Dose*. + I(Dose^2), data=train.AX, weights = w.D))
  a.D <- -0.5*(b.D[2]+as.matrix(test.X)%*%cbind(b.D[17:29]))/b.D[16]
  a.D <- pmin(a.D, 95)
  a.D <- pmax(a.D, 5.81)
  yhat.D <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.D,test.X))[1]$pred)
  yhat.D <- pmin(yhat.D, 0)
  yhat.D <- pmax(yhat.D, -1.2)
  v.D <- mean(yhat.D)
  
  c(Robinson = v.R, BlipTo35 = v.35, Qlearning = v.Q, DWOLS = v.D, weights = fit.Y.AX$coef)
}

# Clinical dosing rule
do.two <- function(Obs, SL.lib = c("SL.mean", "SL.glm"), n.train = 800) {
  # sample size
  N <- dim(Obs)[1]
  
  # choose IDs for training (n = 800) and testing datasets
  train.id <- sample(1:N,n.train)
  test.id <- setdiff(c(1:N),train.id)
  
  # Extract outcome, dose, and covariate data frames
  Y <- -abs(Obs$INR - 2.5)
  A0 <- as.numeric(Obs$Dose == 35)
  A1 <- Obs$Dose
  A2 <- Obs$Dose^2
  X <- Obs[,setdiff(names(Obs), c("Dose","INR"))]
  AX <- Obs[,setdiff(names(Obs), c("INR"))]
  
  clin <- Obs[,c("Age", "Weight", "Height", "Enzyme", "Amiodarone", "Gender", 
                   "Black", "Asian")]
  
  # Test/train set outcomes
  test.Y <- Y[test.id]
  train.Y <- Y[train.id]
  
  # Test/train set dose variables
  test.A0 <- A0[test.id]
  train.A0 <- A0[train.id]
  test.A1 <- A1[test.id]
  train.A1 <- A1[train.id]
  test.A2 <- A2[test.id]
  train.A2 <- A2[train.id]
  
  # Test/train set covariates
  test.X <- X[test.id,]
  train.X <- X[train.id,]
  test.AX <- AX[test.id,]
  train.AX <- AX[train.id,]
  train.clin <- clin[train.id,]
  test.clin <- clin[test.id,]
  
  # Nuisance regressions
  fit.Y.AX <- SuperLearner(train.Y, train.AX, SL.library = c("SL.glmnet.blip", SL.lib))
  fit.Y.0X <- SuperLearner(train.Y[train.A0==1], train.X[train.A0==1,],
                           newX = train.X, SL.library = SL.lib)
  fit.Y.X <- SuperLearner(train.Y, train.X, SL.library = SL.lib)
  fit.A0.X <- SuperLearner(train.A0, train.X, SL.library = SL.lib)
  fit.A1.X <- SuperLearner(train.A1, train.X, SL.library = SL.lib)
  fit.A2.X <- SuperLearner(train.A2, train.X, SL.library = SL.lib)
  
  # Robinson estimator
  Y.X <- train.Y - fit.Y.X$SL.predict[,1]
  A1.X <- train.A1 - fit.A1.X$SL.predict[,1]
  A2.X <- train.A2 - fit.A2.X$SL.predict[,1]
  A1.mat <- A1.X * cbind(1, train.clin)
  names(A1.mat) <- paste0("A1.", names(A1.mat))
  A2.mat <- A2.X
  names(A2.mat) <- paste0("A2.", names(A2.mat))
  A.mat <- cbind(A1.mat, A2.mat)
  
  b.R <- coef(lm(Y.X ~ -1 + ., data = A.mat))
  a.R <- -0.5*(b.R[1]+as.matrix(test.clin)%*%cbind(b.R[2:9]))/b.R[10]
  a.R <- pmin(a.R, 95)
  a.R <- pmax(a.R, 5.81)
  
  yhat.R <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.R,test.X))[1]$pred)
  yhat.R <- pmin(yhat.R, 0)
  yhat.R <- pmax(yhat.R, -1.2)
  v.R <- mean(yhat.R)
  
  # Blip-to-35 estimator
  Y.0 <- train.Y - fit.Y.0X$SL.predict[,1]
  A1.X <- train.A1 - 35
  A2.X <- train.A2 - 35^2
  A1.mat <- A1.X * cbind(1, train.clin)
  names(A1.mat) <- paste0("A1.", names(A1.mat))
  A2.mat <- A2.X
  names(A2.mat) <- paste0("A2.", names(A2.mat))
  A.mat <- cbind(A1.mat, A2.mat)
  
  b.35 <- coef(lm(Y.0 ~ -1 + ., data = A.mat)) #plug-in
  
  wts <- train.A0 / fit.A0.X$SL.predict[,1]
  a1.pred <- fit.A1.X$SL.predict[,1]
  a2.pred <- fit.A2.X$SL.predict[,1]
  A.X.mat <- as.matrix(cbind(a1.pred, a1.pred*train.clin, a2.pred))
  wtsY.0 <- cbind(wts*Y.0)
  one.step <- -c(chol2inv(chol(t(as.matrix(A.mat)) %*% as.matrix(A.mat))) %*% t(A.X.mat) %*% wtsY.0)
  b.35 <- b.35 + one.step #one-step
  
  a.35 <- -0.5*(b.35[1]+as.matrix(test.clin)%*%cbind(b.35[2:9]))/b.35[10]
  a.35 <- pmin(a.35, 95)
  a.35 <- pmax(a.35, 5.81)
  
  yhat.35 <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.35,test.X))[1]$pred)
  yhat.35 <- pmin(yhat.35, 0)
  yhat.35 <- pmax(yhat.35, -1.2)
  v.35 <- mean(yhat.35)
  
  # Q-learning (linear regression)
  b.Q <- coef(lm(train.Y ~ . + Dose:Age + Dose:Weight + Dose:Height + 
                   Dose:Enzyme + Dose:Amiodarone + Dose:Gender + Dose:Black + 
                   Dose:Asian + I(Dose^2), data=train.AX))
  print(b.Q)
  a.Q <- -0.5*(b.Q[15]+as.matrix(test.clin)%*%cbind(b.Q[17:24]))/b.Q[16]
  a.Q <- pmin(a.Q, 95)
  a.Q <- pmax(a.Q, 5.81)
  
  yhat.Q <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.Q,test.X))[1]$pred)
  yhat.Q <- pmin(yhat.Q, 0)
  yhat.Q <- pmax(yhat.Q, -1.2)
  v.Q <- mean(yhat.Q)
  
  # DWOLS (IPW, binned quantiles w/ordered logistic regression)
  w.D <- w.m.fcn(train.A1~1, list(train.A1 ~ .), train.X, 20) # Weights
  b.D <- coef(lm(train.Y ~ . + Dose:Age + Dose:Weight + Dose:Height + 
                   Dose:Enzyme + Dose:Amiodarone + Dose:Gender + Dose:Black + 
                   Dose:Asian + I(Dose^2), data=train.AX, weights = w.D))
  a.D <- -0.5*(b.D[15]+as.matrix(test.clin)%*%cbind(b.D[17:24]))/b.D[16]
  a.D <- pmin(a.D, 95)
  a.D <- pmax(a.D, 5.81)
  yhat.D <- c(predict(fit.Y.AX,newdata=data.frame(Dose=a.D,test.X))[1]$pred)
  yhat.D <- pmin(yhat.D, 0)
  yhat.D <- pmax(yhat.D, -1.2)
  v.D <- mean(yhat.D)
  
  c(Robinson = v.R, BlipTo35 = v.35, Qlearning = v.Q, DWOLS = v.D, weights = fit.Y.AX$coef)
}

### Repeated evaluations

num_to_parallelize <- detectCores(all.tests = FALSE, logical = TRUE)-1 #find out how many cores are on the machine
cl <- makeCluster(num_to_parallelize)
registerDoParallel(cl)

registerDoRNG(2022)
start.time <- Sys.time()
result <- foreach (i = 1:200, .combine = rbind, 
                   .packages = c('SuperLearner','MASS','dplyr'), .export = c('SL.glmnet.blip')) %dopar%
  {
    myLib <- c("SL.mean","SL.glm.interaction","SL.nnet","SL.rpart",
               "SL.glmnet","SL.earth")
    do.one(traindata, SL.lib = myLib, n.train = 800)
  }
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

sprintf("%.4f", apply(result[,1:4], 2, mean))
sprintf("%.4f", apply(result[,1:4], 2, sd))



par(mar=c(4.1,4.1,1.1,1.1))
SL.wts <- apply(result[,-(1:4)],2,mean)
barplot(SL.wts[rev(order(SL.wts))], 
        names.arg = c("LASSO+ixn", "GLM", "Mean", "GLM+ixn", 
                       "MARS", "Neural net", "SVM", "RPART", 
                       "LASSO")[rev(order(SL.wts))],  las=1, col=2, ylab = "Weight")

# Run the clinical model? 

registerDoRNG(2022)
start.time <- Sys.time()
result <- foreach (i = 1:200, .combine = rbind, 
                   .packages = c('SuperLearner','MASS','dplyr'), .export = c('SL.glmnet.blip')) %dopar%
  {
    myLib <- c("SL.mean","SL.glm.interaction","SL.nnet","SL.rpart",
               "SL.glmnet","SL.earth")
    do.two(traindata, SL.lib = myLib, n.train = 800)
  }
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

sprintf("%.4f", apply(result[,1:4], 2, mean))
sprintf("%.4f", apply(result[,1:4], 2, sd))
