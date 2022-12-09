library(SuperLearner)

# EARL implementation for simulations

earl <- function(y, a, x, Q1, Q0, g, loss = "logistic", lambdas = 0, cvFolds = 0L) {
  D1 <- a*(y-Q1)/g + Q1
  D0 <- (1-a)*(y-Q0)/(1-g) + Q0
  s1 <- sign(D1)
  s0 <- sign(D0)
  X <- cbind(1, x)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(lambdas)
  
  if ( loss == "logistic" ) {
    phi <- function(t) log(1 + exp(-t))
  } else if ( loss == "exponential" ) {
    phi <- function(t) exp(-t)
  }
  
  if ( length(lambdas) == 1L ) {
    lambda <- lambdas
    fn <- function(b) {
      eta <- as.matrix(X) %*% cbind(b)
      loss1 <- abs(D1)*phi(s1*eta)
      loss0 <- abs(D0)*phi(-s0*eta)
      risk <- mean(loss0 + loss1) + lambda*sum(b^2)
      return(risk)
    }
    res <- optim(par = rep(0, p), fn = fn)
    beta <- res$par
    return(beta)
    
  } else {
    fold_id <- sample(rep_len(1:cvFolds, length.out = n))
    risks <- matrix(NA, nrow = cvFolds, ncol = k)
    for( i in 1:cvFolds ) {
      # Training data for stage 1 regressions
      X_train <- X[fold_id != i,]
      D0_train <- D0[fold_id != i]
      D1_train <- D1[fold_id != i]
      s0_train <- s0[fold_id != i]
      s1_train <- s1[fold_id != i]
      
      # Test data for stage 2 regressions
      X_test <- X[fold_id == i,]
      D0_test <- D0[fold_id == i]
      D1_test <- D1[fold_id == i]
      s0_test <- s0[fold_id == i]
      s1_test <- s1[fold_id == i]
      
      for( j in 1:k ){
        lambda <- lambdas[j]
        fn <- function(b) {
          eta <- as.matrix(X_train) %*% cbind(b)
          loss1 <- abs(D1_train)*phi(s1_train*eta)
          loss0 <- abs(D0_train)*phi(-s0_train*eta)
          risk <- mean(loss0 + loss1) + lambda*sum(b^2)
          return(risk)
        }
        res <- optim(par = rep(0, p), fn = fn)
        beta <- res$par
        eta <- as.matrix(X_test) %*% cbind(beta)
        loss1 <- abs(D1_test)*phi(s1_test*eta)
        loss0 <- abs(D0_test)*phi(-s0_test*eta)
        risks[i,j] <- mean(loss0 + loss1) + lambda*sum(beta^2)
        # print(mean(loss0 + loss1) + lambda*sum(beta^2))
      }
      # print(risks)
    }
    
    cvRisk <- colSums(risks)
    print(cvRisk)
    lambda.min <- lambdas[which.min(cvRisk)]
    # print(lambda.min)
    fn <- function(b) {
      eta <- as.matrix(X) %*% cbind(b)
      loss1 <- abs(D1)*phi(s1*eta)
      loss0 <- abs(D0)*phi(-s0*eta)
      risk <- mean(loss0 + loss1) + lambda.min*sum(b^2)
      return(risk)
    }
    # print("Made it to last optim")
    res <- optim(par = rep(0, p), fn = fn)
    beta <- res$par
    return(beta)
  }
}

crossfit.nuisances <- function(y, a, x, SL_Y = NULL, SL_A = NULL, V = 2, 
                               SL.cv = list(V = 5L)) {
  # Generate V folds for sample-split/cross-fit steps
  n <- length(y)
  folds <- 1:V
  fold_id <- sample(rep_len(folds, length.out = n))
  
  # Initialize nuisance parameter predictions
  Q0 <- Q1 <- g <- numeric(n)
  
  for( i in folds ){
    # Training data for stage 1 regressions
    x_train <- subset(x, fold_id != i)
    a_train <- a[fold_id != i]
    y_train <- y[fold_id != i]
    
    # Test data for stage 2 regressions
    x_test <- subset(x, fold_id == i)
    a_test <- a[fold_id == i]
    y_test <- y[fold_id == i]
    
    # First stage estimates of nuisance parameters
    Q0_fit <- SuperLearner(Y = y_train[a_train==0], 
                           X = subset(x_train, a_train == 0), newX = x_test, 
                           SL.library = SL_Y, cvControl = SL.cv)
    Q1_fit <- SuperLearner(Y = y_train[a_train==1], 
                           X = subset(x_train, a_train == 1), newX = x_test, 
                           SL.library = SL_Y, cvControl = SL.cv)
    g_fit <- SuperLearner(Y = a_train, X = x_train, newX = x_test, 
                          family = binomial(), SL.library = SL_A, 
                          cvControl = SL.cv)
    
    # Predicted values from stage 1 estimates
    Q0[fold_id == i] <- c(Q0_fit$SL.predict)
    Q1[fold_id == i] <- c(Q1_fit$SL.predict)
    g[fold_id == i] <- c(g_fit$SL.predict)
  }
  
  # Observed and marginal outcome regression
  Qa <- Q0 + a*(Q1 - Q0)
  m <- Q0 + g*(Q1-Q0)
  
  # Store results
  nuisances <- data.frame(Q0, Q1, Qa, g, m)
  attr(nuisances, "fold_id") <- fold_id
  
  return(nuisances)
}

SuperCATE <- function(y, a, x, SL_Y = NULL, SL_A = NULL, SL_CATE = NULL, V = 2, 
                      SL.cv = list(V = 5L)){
  
  # Split the sample, estimate nuisances, and cross-fit for predictions
  nuis <- crossfit.nuisances(y, a, x, SL_Y, SL_A, V, SL.cv)
  Q0 <- nuis$Q0
  Q1 <- nuis$Q1
  Qa <- nuis$Qa
  m <- nuis$m
  g <- nuis$g
  
  # # Generate V folds for sample-split/cross-fit steps
  # n <- length(y)
  # folds <- 1:V
  # fold_id <- sample(rep_len(folds, length.out = n))
  # 
  # 
  # 
  # 
  # # Initialize nuisance parameter predictions
  # Q0 <- Q1 <- g <- m <- Qa <- numeric(n)
  # 
  # for( i in folds ){
  #   # Training data for stage 1 regressions
  #   x_train <- subset(x, fold_id != i)
  #   a_train <- a[fold_id != i]
  #   y_train <- y[fold_id != i]
  #   
  #   # Test data for stage 2 regressions
  #   x_test <- subset(x, fold_id == i)
  #   a_test <- a[fold_id == i]
  #   y_test <- y[fold_id == i]
  #   
  #   # First stage estimates of nuisance parameters
  #   Q0_fit <- SuperLearner(Y = y_train[a_train==0], 
  #                          X = subset(x_train, a_train == 0), newX = x_test, 
  #                          SL.library = SL_Y, cvControl = SL.cv)
  #   Q1_fit <- SuperLearner(Y = y_train[a_train==1], 
  #                          X = subset(x_train, a_train == 1), newX = x_test, 
  #                          SL.library = SL_Y, cvControl = SL.cv)
  #   g_fit <- SuperLearner(Y = a_train, X = x_train, newX = x_test, 
  #                         family = binomial(), SL.library = SL_A, 
  #                         cvControl = SL.cv)
  #   
  #   # Predicted values from stage 1 estimates
  #   Q0[fold_id == i] <- c(Q0_fit$SL.predict)
  #   Q1[fold_id == i] <- c(Q1_fit$SL.predict)
  #   g[fold_id == i] <- c(g_fit$SL.predict)
  # }
  
  ### DR-Learner: regression of doubly-robust pseudo-outcomes
  # Qa <- Q0 + a*(Q1 - Q0)
  
  dr_pseudo <- (a - g)/(g*(1 - g))*(y - Qa) + Q1 - Q0
  dr_learner <- SuperLearner(Y = dr_pseudo, X = x, newX = x, 
                             SL.library = SL_CATE, cvControl = SL.cv)
  dr_learner_lm <- SuperLearner(Y = dr_pseudo, X = x, newX = x, 
                                SL.library = "SL.speedlm", 
                                cvControl = SL.cv)
  
  ### R-learner: regression of "Robinson" pseudo-outcomes with weights
  # m <- Q0 + g*(Q1-Q0)
  
  y_tilde <- y - m
  a_tilde <- a - g
  r_pseudo <- y_tilde/a_tilde
  r_weights <- c(a_tilde)^2
  r_learner <- SuperLearner(Y = r_pseudo, X = x, newX = x, 
                            obsWeights = r_weights, SL.library = SL_CATE, 
                            cvControl = SL.cv)
  
  ### CALM: regression referenced to the treatment free group
  y0 <- y - Q0
  X <- model.matrix(~., data = x)
  x0 <- a * X
  calm_inv <- chol2inv((chol(t(x0) %*% x0)))
  calm_plugin <- calm_inv %*% t(x0) %*% y0 
  eps0 <- (a == 0)*g/(1-g)*y0
  calm_bias <- calm_inv %*% t(X) %*% eps0
  calm_coefs <- calm_plugin - calm_bias
  
  ### DR-PLAM: doubly robust partially linear regression estimates
  A <- diag(a)
  A_tilde <- diag(a_tilde)
  dr_plam_coefs <- chol2inv(chol(t(X) %*% A %*% A_tilde %*% X)) %*% 
    t(X) %*% A_tilde %*% y0
  
  ### R-PLAM: Robinson's partially linear regression estimates
  r_plam_coefs <- chol2inv(chol(t(X) %*% A_tilde %*% A_tilde %*% X)) %*% 
    t(X) %*% A_tilde %*% y_tilde
  
  ### dWOLS: dynamically weighted ordinary least squares
  # Inverse probability weights (i) via Super Learner; and (ii) via GLM
  X_blip <- model.frame(~.-1, data = x)
  names(X_blip) <- paste0("a:",names(X_blip))
  df_dwols <- cbind(X,a=a,a*X_blip)
  X_dwols <- as.matrix(df_dwols) # design matrix
  
  # (i) Super Learner propensity score
  W <- diag(1/(a*g+(1-a)*(1-g))) # weight matrix
  dwols_fit <- chol2inv(chol(t(X_dwols) %*% W %*% X_dwols)) %*% t(X_dwols) %*%
    W %*% y
  dwols_coefs <- dwols_fit[grep("a", names(df_dwols))]
  
  # (ii) Logistic regression propensity score
  g_glm <- predict(glm(a ~ ., data = x, family = binomial()), 
                   x, type = "response")
  W_glm <- diag(1/(a*g_glm+(1-a)*(1-g_glm))) # weight matrix
  dwols_fit_glm <- chol2inv(chol(t(X_dwols) %*% W_glm %*% X_dwols)) %*% 
    t(X_dwols) %*% W_glm %*% y
  dwols_coefs_glm <- dwols_fit[grep("a", names(df_dwols))]
  
  # Efficient augmentation and relaxation learning (unregularized)
  EARL_logit <- earl(y, a, x, Q1, Q0, g, loss = "logistic")
  
  # Efficient augmentation and relaxation learning (regularized)
  EARL_expo <- earl(y, a, x, Q1, Q0, g, loss = "logistic", lambdas = 2^{-5:5}, cvFolds = 5L)
  
  # Estimate Conditional Average Treatment Effect
  CATE_fit <- list(DR.learner = dr_learner, 
                   DR.learner.GLM = dr_learner_lm, 
                   R.learner = r_learner, CALM = calm_coefs,
                   DR.PLAM = dr_plam_coefs, R.PLAM = r_plam_coefs, 
                   DWOLS.SL = dwols_coefs, DWOLS.GLM = dwols_coefs_glm,
                   EARL.logit = EARL_logit, EARL.expo = EARL_expo)
  
  return(CATE_fit)
}

###=|=|=| Testing |=|=|=###
# set.seed(1)
# O <- DGP3()
# SL_lib <- c("SL.mean","SL.earth", "SL.glmnet", "SL.rpart", "SL.glm", "SL.glm.interaction", "SL.nnet")
# SL_lib <- c("SL.glm")
# cate <- SuperCATE(O$Y, O$A, O$X, SL_Y=SL_lib, SL_A=SL_lib, SL_CATE=SL_lib, V = 5L, SL.cv = list(V=5L))
# cate



