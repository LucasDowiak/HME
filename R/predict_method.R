# predict hme object
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

predict.hme <- function(obj, newdata)
{
  
  form <- Formula::as.Formula(obj[["call_"]][["formula"]])
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = newdata)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  Z <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  pred <- internal_predict(Y, X, Z, obj[["tree"]], obj[["expert.nms"]], obj[["gate.nodes"]],
                           obj[["gate.pars"]], obj[["expert.pars"]], obj[["expert.type"]])
  return(pred)
}

internal_predict <- function(X, Z, expertnames, gatenodes, gatepars,
                             expertpars, experttype)
{
  list_nodes <- napply(gatenodes, par_to_gate_paths, gatepars, Z)
  prior_exp_wts <- napply(expertnames, prior_weights, list_nodes)
  preds <-  lapply(expertpars, expert_pred, X, experttype)
  return(c(Reduce(`+`, mapply(`*`, preds, prior_exp_wts, SIMPLIFY=FALSE))))
}

