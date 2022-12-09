
library(DynTxRegime)
library(SuperLearner)
library(doParallel)
library(doRNG)

fpre <- "6Dec22"

work.dir <- "/Users/david/Documents/GitHub/calm-optimized"
out.dir <- "/Users/david/Documents/GitHub/calm-optimized/output/"

setwd(work.dir)

################################################################################

##### HELPER WRAPPER FUNCTIONS
my.SuperLearner <- function(Y, X, newX = NULL, family = gaussian(), SL.library, 
                            method = "method.NNLS", id = NULL, verbose = FALSE, 
                            control = list(), cvControl = list(), 
                            obsWeights = NULL, env = parent.frame()) {
  Y <- c(Y)
  X <- as.data.frame(X)
  out <- SuperLearner(Y, X, newX, family, SL.library, method, id, verbose, control, 
               cvControl, obsWeights, env)
  class(out) <- c("SuperLearner", "mySuperLearner")
  return(out)
}

predict.mySuperLearner <- function(object, ...) {
  predict.SuperLearner(object, ...)[1]$pred
}

#######################################
#####                             #####
#####      MAIN EARL WRAPPER      #####
#####                             #####
#######################################

my.earl <- function(Y, A, X, SL.lib) {
  namX <- names(as.data.frame(X))
  formX <- as.formula(paste("~", paste(namX, collapse = "+")))
  
  O <- data.frame(A = A, X)
  
  # propensity model
  moPropen <- buildModelObj(model = formX,
                            solver.method = 'my.SuperLearner',
                            solver.args = list('X' = 'x', 'Y' = 'y',
                                               'family'='binomial',
                                               'SL.library' = SL.lib),
                            predict.method = 'predict.mySuperLearner',
                            predict.args = NULL)
  
  # outcome model
  moMain <- buildModelObj(model = formX,
                          solver.method = 'my.SuperLearner',
                          solver.args = list('X' = 'x', 'Y' = 'y',
                                             'family'='gaussian',
                                             'SL.library' = SL.lib),
                          predict.method = 'predict.mySuperLearner',
                          predict.args = NULL)
  moCont <- buildModelObj(model = formX,
                          solver.method = 'my.SuperLearner',
                          solver.args = list('X' = 'x', 'Y' = 'y',
                                             'family'='gaussian',
                                             'SL.library' = SL.lib),
                          predict.method = 'predict.mySuperLearner',
                          predict.args = NULL)
  
  # EARL
  fitEARL <- earl(moPropen = moPropen, moMain = moMain, moCont = moCont,
                  data = O, response = Y,  txName = 'A', 
                  lambdas = 2^(-5:5), cvFolds = 10L,
                  regime = formX, surrogate = 'logit', kernel = 'linear',
                  verbose = 0)
  regimeCoef(fitEARL)
}

################################################################################

optimal.rules <- function(object, newdata){
  n <- dim(newdata)[1]
  Xmat <- model.matrix(~ ., newdata)
  
  tmp <- numeric(n)
  rules <- data.frame(EARL = tmp)
  
  # EARL
  EARL <-  Xmat %*% object
  rules$EARL <- c(1L*(EARL >= 0))
  
  return(rules)
}

earl_rules <- function(train, eval, SL_lib){
  X <- train$X
  A <- train$A
  Y <- train$Y
  Xnew <- eval$X
  cate <- my.earl(Y, A, X, SL.lib = SL_lib)
  out <- optimal.rules(cate, Xnew)
  return(out)
}

eval_value <- function(rules, eval){
  rules <- data.frame(opt = eval$dopt, rules)
  eval_df <- data.frame(Q0=eval$Q0, CATE=eval$CATE, e=eval$e)
  value <- function(rule){
    mean(eval_df$Q0 + eval_df$CATE*rule + eval_df$e)
  }
  out <- apply(rules, 2, value)
  return(out)
}

eval_frac <- function(rules, eval){
  frac <- function(rule){
    mean(eval$dopt == rule)
  }
  out <- apply(rules, 2, frac)
  return(out)
}

################################################################################

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

DGP1 <- function(n) {DGP_proto(n, linear = c(T,T,T))}
DGP2 <- function(n) {DGP_proto(n, linear = c(T,T,F))}
DGP3 <- function(n) {DGP_proto(n, linear = c(T,F,T))}
DGP4 <- function(n) {DGP_proto(n, linear = c(T,F,F))}
DGP5 <- function(n) {DGP_proto(n, linear = c(F,T,T))}
DGP6 <- function(n) {DGP_proto(n, linear = c(F,T,F))}
DGP7 <- function(n) {DGP_proto(n, linear = c(F,F,T))}
DGP8 <- function(n) {DGP_proto(n, linear = c(F,F,F))}

################################################################################

do_one <- function(n = 300L, DGP = DGP1, SL_lib = "SL.glm"){
  train <- DGP(n)
  eval <- DGP(n = 1e4)
  d <- earl_rules(train, eval, SL_lib)
  
  out_value <- eval_value(d, eval)
  out_frac <- eval_frac(d, eval)
  out <- c(value = out_value, frac = out_frac)
  return(out)
}

################################################################################

#find out how many cores are on the machine
num_to_parallelize <- detectCores(all.tests = FALSE, logical = TRUE)-1
cl <- makeCluster(num_to_parallelize)
registerDoParallel(cl)

args <- expand.grid(DGP = c(1:8), n = c(500L, 1000L, 2000L, 4000L))
n_args <- dim(args)[1]
n_studies <- 200L

clusterExport(cl, c("my.SuperLearner","predict.mySuperLearner"))

for(j in seq_len(n_args)){
  registerDoRNG(2022)
  start.time <- Sys.time()
  result <- foreach (i = 1:n_studies, .combine = rbind, 
                     .packages = c('SuperLearner','DynTxRegime','rpart'),
                     .export = c()) %dopar%
    {
      # Lib <- c("SL.mean","SL.glm.interaction","SL.nnet","SL.rpart","SL.glmnet","SL.earth")
      ## rpart leads to Error: task 1 failed - "prediction method could not be executed successfully"
      Lib <- c("SL.mean","SL.glm.interaction","SL.nnet","SL.glmnet","SL.earth")
      DGPs <- list(DGP1, DGP2, DGP3, DGP4, DGP5, DGP6, DGP7, DGP8)
      do_one(n = args$n[j], DGP = DGPs[[args$DGP[j]]], SL_lib = Lib)
    }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste("Setting",j,"of",n_args,"completed"))
  print(time.taken)
  fname <- paste0(out.dir, paste0("earl_", fpre,"_", args$n[j], "_", args$DGP[j], ".csv"))
  write.csv(t(result), fname)
}
