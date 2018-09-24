logLik.hme <- function(obj)
{
  logL <- obj[["logL"]][!is.na(tst[["logL"]])]
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
    call2 <- c(obj[["call_"]][-1], grab_vitals(obj, ee))
    buds[[ee]] <- do.call(hme, call2)
  }
  ll <- which.max(sapply(buds, logLik))
  newtree <- c(obj[["tree"]], paste(names(ll), 1:2, sep="."))
  return(newtree)
}


