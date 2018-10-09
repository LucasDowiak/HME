logLik.hme <- function(obj)
{
  m <- c(obj[["logL"]])
  logL <- m[!is.na(m)]
  return(logL[length(logL)])
}


criterion <- function(obj, type=c("aic", "bic"))
{
  type <- match.arg(type)
  L <- logLik(obj)
  K <- obj[["no.of.pars"]]
  
  if (type == "aic") {
    penalty <- 2
  } else if (type == "bic") {
    penalty <- log(obj[["N"]])
  }
  return(penalty * K - 2 * L)
}


grab_vitals <- function(obj, expert)
{
  # gate product to expert to pass to root_posterior
  root_prior <- gate_path_product("0", expert, obj$list_priors)
  # two new expert parameters (from fitted values) (gaussian deviations with std errors)
  pars <- obj[["expert.pars"]][[expert]]
  newpars <- list(rnorm(length(pars), mean=pars, sd=.1),
                  rnorm(length(pars), mean=pars, sd=.1))
  names(newpars) <- c("0.1", "0.2")
  return(list(root_prior=root_prior, init_expert_pars=newpars))
}


grow_the_tree <- function(obj)
{
  call_ <- obj[["call_"]]
  expnms <- obj[["expert.nms"]]
  buds <- vector("list", length(expnms))
  names(buds) <- expnms
  for (ee in expnms) {
    newcall <- c(obj[["call_"]][-1], grab_vitals(obj, ee))
    buds[[ee]] <- do.call(hme, newcall)
  }
  ll <- which.max(sapply(buds, logLik))
  
  new_gate_pars <- buds[[ll]][["gate.pars"]]
  names(new_gate_pars) <- expnms[ll]
  new_gate_pars <- c(obj[["gate.pars"]], new_gate_pars)[sort(c(obj[["gate.nodes"]], expnms[ll]))]
  newtree <- c(obj[["tree"]], paste(names(ll), 1:2, sep="."))
  return(list(tree=newtree, init_gate_pars=new_gate_pars))
}


bootstrap_glm <- function(n=2, ...)
{
  object <- glm(...)
  df.r <- object$df.residual
  dispersion <- if (object$family$family %in% c("poisson", "binomial")) 
    1
  else if (df.r > 0) {
    est.disp <- TRUE
    if (any(object$weights == 0)) 
      warning("observations with zero weight not used for calculating dispersion")
    sum((object$weights * object$residuals^2)[object$weights > 
                                                0])/df.r
  }
  else {
    est.disp <- TRUE
    NaN
  }
  p <- object$rank
  p1 <- 1L:p
  Qr <- stats:::qr.lm(object)
  covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  covmat <- dispersion * covmat.unscaled
  
  replicate(n,
            rnorm(p, mean=object$coefficients, sd = sqrt(diag(covmat))),
            simplify=FALSE)
}




