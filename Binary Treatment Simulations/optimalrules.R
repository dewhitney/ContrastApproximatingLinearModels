
optimal.rules <- function(object, newdata){
  n <- dim(newdata)[1]
  Xmat <- model.matrix(~ ., newdata)
  
  tmp <- numeric(n)
  rules <- data.frame(DR.learner = tmp, DR.learner.GLM = tmp, R.learner = tmp, 
                      CALM = tmp, DR.PLAM = tmp, R.PLAM = tmp, DWOLS.SL = tmp, 
                      DWOLS.GLM = tmp)
  
  # DR-Learner with Super Learner
  DR.learner <- predict(object$DR.learner, newdata = newdata)$pred
  rules$DR.learner <- c(1L*(DR.learner >= 0))
  
  # DR-Learner with least squares projection
  DR.learner.GLM <- predict(object$DR.learner.GLM, newdata = newdata)$pred
  rules$DR.learner.GLM <- c(1L*(DR.learner.GLM >= 0))
  
  # R-Learner with Super Learner
  R.learner <- predict(object$R.learner, newdata = newdata)$pred
  rules$R.learner <- c(1L*(R.learner >= 0))
  
  # CALM
  CALM <- Xmat %*% object$CALM
  rules$CALM <- c(1L*(CALM >= 0))
  
  # DR.PLAM
  DR.PLAM <- Xmat %*% object$DR.PLAM
  rules$DR.PLAM <- c(1L*(DR.PLAM >= 0))
  
  # R.PLAM
  R.PLAM <- Xmat %*% object$R.PLAM
  rules$R.PLAM <- c(1L*(R.PLAM >= 0))
  
  # DWOLS.SL
  DWOLS.SL <- Xmat %*% object$DWOLS.SL
  rules$DWOLS.SL <- c(1L*(DWOLS.SL >= 0))
  
  # DWOLS.GLM
  DWOLS.GLM <- Xmat %*% object$DWOLS.GLM
  rules$DWOLS.GLM <- c(1L*(DWOLS.GLM >= 0))
  
  # EARL (logistic surrogate loss)
  EARL.logit <-  Xmat %*% object$EARL.logit
  rules$EARL.logit <- c(1L*(EARL.logit >= 0))
  
  # EARL (exponential surrogate loss)
  EARL.expo <-  Xmat %*% object$EARL.expo
  rules$EARL.expo <- c(1L*(EARL.expo >= 0))
  
  return(rules)
}
