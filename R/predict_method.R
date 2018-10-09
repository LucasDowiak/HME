# predict hme object

predict.hme <- function(obj, newdata)
{
  
  form <- Formula::as.Formula(obj[["call_"]][["formula"]])
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = newdata)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  Z <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  
  list_nodes <- napply(obj[["gate.nodes"]], par_to_gate_paths, obj[["gate.pars"]], Z)
  prior_exp_wts <- napply(obj[["expert.nms"]], prior_weights, obj[["tree"]], list_nodes)
  preds <-  lapply(obj[["expert.pars"]], expert_pred, X, obj[["expert.type"]])
  return(c(Reduce(`+`, mapply(`*`, preds, prior_exp_wts, SIMPLIFY=FALSE))))
  
}


expert_pred <- function(beta, X, expert_type)
{
  if (expert_type == "gaussian") {
    beta <- as.matrix(beta[-length(beta)])
    return(X %*% beta)
  } else {
    mssg <- "predict.hme cannot handle expert type `%s`."
    stop(mssg, expert_type)
  }
}
