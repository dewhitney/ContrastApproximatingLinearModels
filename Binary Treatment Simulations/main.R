library(doParallel)
library(doRNG)

fpre <- "6Dec22"

work.dir <- "/Users/david/Documents/GitHub/calm-optimized"
out.dir <- "/Users/david/Documents/GitHub/calm-optimized/output/"

setwd(work.dir)

source("supercate.R")
source("optimalrules.R")
source("evalrules.R")
source("dgp.R")

do_one <- function(n = 300L, DGP = DGP1a, SL_lib = "SL.glm", V = 10L){
  train <- DGP(n)
  eval <- DGP(n = 1e4)
  d <- eval_rules(train, eval, SL_lib, V = V)
  
  out_value <- eval_value(d, eval)
  out_frac <- eval_frac(d, eval)
  out <- c(value = out_value, frac = out_frac)
  return(out)
}

num_to_parallelize <- detectCores(all.tests = FALSE, logical = TRUE)-1 #find out how many cores are on the machine
cl <- makeCluster(num_to_parallelize)
registerDoParallel(cl)

# args <- expand.grid(DGP = c(1:3,5:9), n = c(500L, 1000L, 2000L, 4000L))
args <- expand.grid(DGP = c(1), n = c(500L))

n_args <- dim(args)[1]
n_studies <- 200L

for(j in seq_len(n_args)){
  registerDoRNG(2022)
  start.time <- Sys.time()
  result <- foreach (i = 1:n_studies, .combine = rbind, 
                     .packages = c('SuperLearner','DynTxRegime'),
                     .export = c('my.earl', 'my.SuperLearner', 
                                 'my.predict.SuperLearner')) %dopar%
    {
      Lib <- c("SL.mean","SL.glm.interaction","SL.nnet","SL.rpart","SL.glmnet","SL.earth")
      DGPs <- list(DGP1, DGP2, DGP3, DGP4, DGP4, DGP5, DGP6, DGP7, DGP8, DGP9)
      do_one(n = args$n[j], DGP = DGPs[[args$DGP[j]]], SL_lib = Lib, V = 5L)
    }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste("Setting",j,"of",n_args,"completed"))
  print(time.taken)
  fname <- paste0(out.dir, paste0(fpre,"_", args$n[j], "_", args$DGP[j], ".csv"))
  write.csv(t(result), fname)
}

