
eval_rules <- function(train, eval, SL_lib, V = 5L){
  X <- train$X
  A <- train$A
  Y <- train$Y
  Xnew <- eval$X
  cate <- SuperCATE(Y, A, X, SL_Y=SL_lib, SL_A=SL_lib, SL_CATE=SL_lib, V = V, SL.cv = list(V=3L))
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