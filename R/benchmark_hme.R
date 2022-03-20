# ----------------------------------------------------------------------------
source("R/building_blocks.R")
source("R/marginal_effects.R")
source("R/optimization_functions.R")
source("R/summary_method.R")
source("R/plot_method.R")
source("R/predict_method.R")
source("R/tree_construction.R")
library(Formula)


hme <- function(tree, formula, hme_type=c("hme", "hmre"),
                expert_type=c("gaussian", "bernoulli"), data, holdout=NULL,
                root_prior=1, init_gate_pars=NULL, init_expert_pars=NULL,
                maxiter=100, tolerance=1e-4, trace=0)
{
  cl <- match.call()
  require(Formula)
  stopifnot(!missing(data))
  nullholdout <- is.null(holdout)
  
  tree <- sort(tree)
  # TODO: You should be able to tell hme_type by looking at `tree`
  hme_type <- match.arg(hme_type)
  expert_type <- match.arg(expert_type)
  
  mf <- match.call(expand.dots = FALSE)
  call_ <- as.list(mf)
  m <- match(c("formula", "data"), names(mf), 0)
  form <- Formula::as.Formula(formula)
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  Z <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  N <- length(Y)
  
  if (!nullholdout) {
    mfp <- model.frame(form, data = holdout)
    Yp <- model.response(mfp, "numeric")
    Xp <- model.matrix(terms(form, data = mfp, rhs = 1), mfp)
    Zp <- model.matrix(terms(form, data = mfp, rhs = 2), mfp)
  }
  
  tree <- sort(tree)
  expert.nodes <- tree[unlist(is_terminal(tree, tree))]
  gate.nodes <- setdiff(tree, expert.nodes)
  expert.index <- expert_index(hme_type, tree)
  
  if (is(init_gate_pars, "NULL")) {
    gate.pars <- napply(gate.nodes, init_gate_node_pars, tree, ncol(Z))
  } else {
    gate.pars <- init_gate_pars
  }
  
  if (is(init_expert_pars, "NULL")) {
    if (expert_type=="gaussian") {
      expert.pars <- bootstrap_glm(length(expert.nodes),
                                   formula=terms(form, data=mf, rhs=1),
                                   data=data,
                                   family=gaussian,
                                   model=FALSE)
      expert.pars <- napply(expert.pars, function(x) c(x, log(var(Y))))
      names(expert.pars) <- expert.nodes
    } else {
      expert.pars <- napply(expert.nodes, function(x) runif(ncol(X) , -2, 2))
    }
  } else {
    expert.pars <- init_expert_pars
  }

  NN <- length(c(unlist(gate.pars), unlist(expert.pars)))
  
  list_nodes <- napply(gate.nodes, par_to_gate_paths, gate.pars, Z)
  list_experts <- napply(expert.nodes, par_to_expert_dens, expert_type,
                          expert.pars, Y, X)

  
  logL <- errorR <- matrix(NA_real_, nrow=maxiter + 1)
  parM <- matrix(NA_real_, ncol=length(unlist(c(gate.pars, expert.pars))),
                 nrow=maxiter + 1)
  parM[1, ] <- unlist(c(gate.pars, expert.pars))
  logL[1, ] <- newLL <- -Inf
  
  for (ii in seq_len(maxiter)) {
    oldLL <- newLL
    mstep <- m_step(tree, hme_type, expert_type, Y=Y, X=X, Z=Z,
                    exp.pars=expert.pars, gat.pars=gate.pars,
                    root_prior=root_prior)
    expert.pars <- mstep[["exp.pars"]]
    gate.pars <- mstep[["gat.pars"]]
    newLL <- mstep[["loglik"]] / length(Y)
    ewmaLL <- if (is.infinite(oldLL)) {
      newLL
    } else {
      0.9 * oldLL + 0.1 * newLL
    }
    logL[ii + 1, ] <- newLL
    parM[ii + 1, ] <- unlist(c(gate.pars, expert.pars))
    if (!nullholdout) {
      yhat <- internal_predict(Xp, Zp, expert.nodes, gate.nodes,
                               gate.pars, expert.pars, expert_type)
      errorR[ii + 1, ] <- mean((Yp - yhat)**2)
    }
    if (trace > 0) {
      cat('\r', sprintf("Step: %d - Log-Likelihood: %f - Weighted LL: %f", ii, newLL, ewmaLL))
    } else {
      cat('\n', sprintf("Step: %d - Log-Likelihood: %f - Weighted LL: %f", ii, newLL, ewmaLL))
    }
    if (abs(oldLL - newLL) < tolerance)
      break
  }
  cat("\n")
  logL <- logL[!is.na(logL), , drop=FALSE]
  parM <- parM[apply(!is.na(parM), 1, FUN=all), ]
  errorR <- errorR[!is.na(errorR), , drop=FALSE]
  
  # Create the full score vector of theta = (omega + beta)
  gate_scores <- napply(gate.nodes, logistic_score, mstep$list_priors, mstep$list_density, Z)
  expt_scores <- napply(expert.nodes, gaussian_score, expert.pars, mstep$list_priors,
                        mstep$list_density, Y, X)
  scores <- do.call(cbind, c(unlist(gate_scores, recursive=FALSE), expt_scores))
  
  # gte.nms, exp.nms, ld, ln, gate.pars, exp.pars, Y, X, Z, N
  # full.vcv <- sandwich_vcov(gate.nodes,
  #                          expert.nodes,
  #                          mstep$list_density,
  #                          mstep$list_priors,
  #                          gate.pars,
  #                          expert.pars,
  #                          Y, X, Z, N)

  structure(list(tree=tree,
                 hme.type=hme_type,
                 gate.nodes=gate.nodes,
                 gate.pars=gate.pars,
                 gate.pars.nms=colnames(Z),
                 expert.type=expert_type,
                 expert.nms=expert.nodes,
                 expert.pars=expert.pars,
                 expert.pars.nms=c(colnames(X), "Dispersion"),
                 dependent.var.nme=all.vars(form)[1],
                 MSE=errorR,
                 logL=logL,
                 parM=parM,
                 list_priors=mstep$list_priors,
                 list_posteriors=mstep$list_posteriors,
                 list_density=mstep$list_density,
                 Y=Y,
                 X=X,
                 Z=Z,
                 N=N,
                 no.of.pars=NN,
                 gate.info.matrix=NA, # deprecated
                 expert.info.matrix=NA, # deprecated
                 full.vcv=NULL, # full.vcv,
                 scores=scores,
                 call_=call_),
            class="hme")
}


refactor_hme <- function(obj, how=c("loglik", "mse"), pars=NULL)
{
  stopifnot(inherits(obj, "hme"))
  how <- match.arg(how)
  
  if (is.null(pars)) {
    if (how == "mse") {
      pars <- obj[["parM"]][which.min(obj[["MSE"]]), ]
    } else {
      pars <- obj[["parM"]][which.max(obj[["logL"]]), ]
    }
  }
  depth <- max(sapply(strsplit(obj[["expert.nms"]], "\\."), length))
  lgpn <- length(obj[["gate.pars.nms"]])
  if (depth == 2 && length(obj[["expert.nms"]]) > 2) {
    lgp <- length(mod[["gate.pars"]][[1]])
    ngates <- lgp * lgpn
    gate.pars <- pars[1:ngates]
    gate.pars <- list(split(gate.pars, ((seq_along(gate.pars) - 1) %/% lgpn) + 1))
  } else {
    lgp <- length(obj[["gate.pars"]])
    ngates <- lgp * lgpn
    gate.pars <- pars[1:ngates]
    gate.pars <- split(gate.pars, ((seq_along(gate.pars) - 1) %/% lgpn) + 1)
  }
  names(gate.pars) <- names(obj[["gate.pars"]])
  lep <- length(obj[["expert.pars"]])
  lepn <- length(obj[["expert.pars.nms"]])
  expert.pars <- pars[(ngates + 1):length(pars)]
  expert.pars <- split(expert.pars, ((seq_along(expert.pars) - 1) %/% lepn))
  names(expert.pars) <- names(obj[["expert.pars"]])
  
  list_priors <- napply(names(gate.pars), par_to_gate_paths, gate.pars, obj[["Z"]])
  list_density <- napply(names(expert.pars), par_to_expert_dens, obj[["expert.type"]],
                         expert.pars, obj[["Y"]], obj[["X"]])
  list_posteriors <- napply(names(gate.pars), posterior_weights, obj[["tree"]],
                            list_priors, list_density)
  lik <- log_likelihood(obj[["tree"]], list_priors, list_posteriors, list_density,
                        1)

  gate.info.matrix <- napply(obj[["gate.nodes"]], multinomial_info_matrix,
                             obj[["tree"]], gate.pars, list_priors,
                             list_posteriors, obj[["Z"]], 1)
  
  expert.info.matrix <- napply(obj[["expert.nms"]], expert_info_matrix,
                               obj[["expert.type"]], expert.pars,
                               list_posteriors, obj[["X"]], obj[["Y"]])
  
  
  obj[["expert.pars"]] <- expert.pars
  obj[["gate.pars"]] <- gate.pars
  obj[["list_posteriors"]] <- list_posteriors
  obj[["list_density"]] <- list_density
  obj[["list_priors"]] <- list_priors
  obj[["gate.info.matrix"]] <- gate.info.matrix
  obj[["expert.info.matrix"]] <- expert.info.matrix
  obj[["logL"]] <- rbind(obj[["logL"]], lik)
  return(obj)
}
