library(SuperLearner)

# calm

# Description:
# Fits contrast-approximating linear models.

# Arguments:
# baseline             - a data.frame of pre-exposure measurements.
# exposure             - a data.frame of exposure measurements (may have  
#                        multiple columns).
# outcome              - a vector of outcome measurements.
# contrast.models      - a list of formula for each exposure specified in terms 
#                        of baseline variables.
# outcome.learners     - a character vector with function names for the outcome 
#                        SuperLearner.
# exposure.learners    - a character vector with function names for the exposure
#                        SuperLearner.
# variance.learners    - a character vector with function names for the variance
#                        SuperLearner.
# zero.learners        - a character vector with function names for the  
#                        probability of unexposed SuperLearner.
# include.blip.to.zero - a boolean indicating whether to include the blip to zero 
#                        approximation.
# cvControl            - a list of arguments to control cross-validation for 
#                        calls SuperLearner. See ?SuperLearner for details.

# Contrast-approximating linear models
calm <- function(baseline, exposure, outcome, contrast.models = list(~ 1), 
                 outcome.learners = "SL.glm", exposure.learners = "SL.glm", 
                 variance.learners = "SL.glm", zero.learners = "SL.glm", 
                 cvControl = list(V = 5)) {
  
  # Unweighted methods: Robinson, blip-to-zero, model-based
  
  # Outcome regression on baseline only (m)
  SL.outcome <- SuperLearner(outcome, baseline, SL.library = outcome.learners,
                             cvControl = cvControl)
  m <- predict(SL.outcome, newdata = baseline)$pred[,1]
  
  # Exposure regression on baseline (e1, e2, ...)
  n.exposures <- ncol(exposure)
  SL.exposure <- e <- vector("list", n.exposures)
  for ( i in 1:n.exposures ) {
    SL.exposure[[i]] <- SuperLearner(exposure[,i], baseline, 
                                     SL.library = outcome.learners, 
                                     cvControl = cvControl)
    e[[i]] <- predict(SL.exposure[[i]], newdata = baseline)$pred[,1]
  }
  
  # Robinson coefficients
  Y.Rob <- outcome - m
  X.model <- X.Robs <- vector("list", n.exposures)
  for ( i in 1:n.exposures ) {
    X.model[[i]] <- model.matrix(contrast.models[[i]], data = baseline)
    if (length(colnames(X.model[[i]])) > 1) {
      colnames(X.model[[i]]) <- c(names(exposure)[i], 
                                  paste(names(exposure)[i], 
                                        colnames(X.model[[i]])[-1], sep = ":"))
    } else {
      colnames(X.model[[i]]) <- names(exposure)[i]
    }
    
    X.Robs[[i]] <- (exposure[,i] - e[[i]])*X.model[[i]]
  }
  X.Rob <- do.call(cbind, X.Robs)
  Beta.Rob <- coef(lm(Y.Rob ~ 0 + ., data = data.frame(X.Rob)))
  
  # Probability of zero exposure group (g)
  zero.group <- apply(as.matrix(exposure), 1, function(x) all(x == 0))
  SL.prob.zero <- SuperLearner(1*zero.group, baseline, 
                               SL.library = zero.learners, 
                               family = "binomial", cvControl = cvControl)
  g <- predict(SL.prob.zero, newdata = baseline)$pred[,1]
  g <- pmax(g, 1e-8)
  
  # Outcome regression on baseline in zero exposure group (M0)
  SL.out.zero <- SuperLearner(outcome[zero.group], baseline[zero.group,],
                              SL.library = outcome.learners,
                              family = "gaussian", cvControl = cvControl)
  M0 <- predict(SL.out.zero, newdata = baseline)$pred[,1]
  
  # Blip-to-zero coefficients
  Y.0 <- outcome - M0
  X.0s <- X.es <- vector("list", n.exposures)
  for ( i in 1:n.exposures ) {
    X.0s[[i]] <- exposure[,i]*X.model[[i]]
    X.es[[i]] <- e[[i]]*X.model[[i]]
  }
  X.0 <- do.call(cbind, X.0s)
  X.e <- do.call(cbind, X.es)
  plugin.0 <- coef(lm(Y.0 ~ 0 + ., data = data.frame(X.0)))
  update.0 <- -chol2inv((chol(t(X.0) %*% X.0))) %*% t(X.e) %*% cbind(zero.group*Y.0/g)
  Beta.0 <- plugin.0 + update.0[,1]
  
  # Weighted methods: Weighted/re-centered Robinson (plug-in and one-step)
  
  # Exposure adjusted outcome regression on baseline (M)
  SL.out.adj <- SuperLearner(outcome, data.frame(exposure, baseline), 
                             SL.library = outcome.learners, 
                             cvControl = cvControl)
  M <- predict(SL.out.adj, newdata = data.frame(exposure, baseline))$pred[,1]
  
  # Variance-estimation
  SL.var <- SuperLearner((outcome - M)^2, data.frame(exposure, baseline), 
                         SL.library = variance.learners, cvControl = cvControl)
  
  # Construct inverse-variance weights
  Var <- predict(SL.var, newdata = data.frame(exposure, baseline))$pred[,1]
  Var <- pmax(Var, min((outcome - M)^2))
  W <- 1/Var
  
  # Expectations for inverse-variance weights (w, mw, ew1, ew2, ...)
  SL.wts <- SuperLearner(W, baseline, SL.library = variance.learners, 
                         cvControl = cvControl)
  w <- predict(SL.wts, newdata = baseline)$pred[,1]
  w <- pmax(w, min(W))
  
  SL.out.wts <- SuperLearner(outcome*W, baseline, SL.library = outcome.learners,
                             cvControl = cvControl)
  mw <- pmax(predict(SL.out.wts, newdata = baseline)$pred[,1], min(outcome*W))/w
  
  SL.exp.wts <- ew <- vector("list", n.exposures)
  for ( i in 1:n.exposures ) {
    SL.exp.wts[[i]] <- SuperLearner(exposure[,i]*W, baseline, 
                                    SL.library = exposure.learners, 
                                    cvControl = cvControl)
    ew[[i]] <- pmax(predict(SL.exp.wts[[i]], newdata = baseline)$pred[,1], 
                    min(exposure[,i]*W))/w
  }
  
  # Weighted/re-centered Robinson (plug-in)
  Y.w <- outcome - mw
  X.ws <- vector("list", n.exposures)
  for ( i in 1:n.exposures ) {
    X.ws[[i]] <- (exposure[,i] - ew[[i]])*X.model[[i]]
  }
  X.w <- do.call(cbind, X.ws)
  Beta.w1 <- coef(lm(Y.w ~ 0 + ., data = data.frame(X.w), weights = W))
  
  # Weighted/re-centered Robinson (one-step)
  update.w <- -chol2inv((chol(t(X.w) %*% diag(W) %*% X.w))) %*% t(X.w) %*% 
    cbind(W^2*((outcome - M)^2-1/W)*(M - mw - X.w %*% Beta.w1))
  Beta.w2 <- Beta.w1 + c(update.w)
  
  return(rbind(Beta.Rob, Beta.0, Beta.w1, Beta.w2))
}
