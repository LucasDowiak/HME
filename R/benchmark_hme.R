# ----------------------------------------------------------------------------
source("building_blocks.R")
library(Formula)

"standard hme approach where each node is updated independently"
m_step <- function(tree, hme_type, Y, X, Z, exp.pars, gat.pars)
{
  gat.out <- lapply(gat.pars, function(x) list())
  exp.out <- lapply(exp.pars, function(x) list())
  
  list_priors <- napply(names(gat.pars), par_to_gate_paths, "binomial",
                         gat.pars, Z)
  list_density <- napply(names(exp.pars), par_to_expert_dens, "gaussian",
                       exp.pars, Y, X)
  list_posteriors <- napply(names(gat.pars), posterior_weights,
                             tree, list_priors, list_density)
  
  for (e in names(exp.pars)) {
    opt_blk <- optimize_block(e, tree, exp.pars, Y, X,
                              joint_posterior_weight(e, tree, list_posteriors),
                              "gaussian", hme_type)
    exp.out[[e]] <- opt_blk$par
  }
  
  for (g in names(gat.pars)) {
    # 1) wt = find joint posterior weights that represents current partition
    wt <- joint_posterior_weight(g, tree, list_posteriors)
    # 2) target = h^{a,i} using g^{a,i} parameters with weights wt
    target <- list_posteriors[[g]][, 1]
    wts=list(branch=target, joint_post=wt)
    opt_blk <- optimize_block(g, tree, gat.pars, NULL, Z, wts,
                              "binomial", hme_type)
    gat.out[[g]] <- opt_blk$par
  }
  list_priors <- napply(names(gat.pars), par_to_gate_paths, "binomial",
                        gat.pars, Z)
  list_density <- napply(names(exp.pars), par_to_expert_dens, "gaussian",
                      exp.pars, Y, X)
  list_posteriors <- napply(names(gat.pars), posterior_weights,
                            tree, list_priors, list_density)
  lik <- sum(posterior_weights("0", tree, list_priors, list_density))
  
  return(list(exp.pars=exp.out, gat.pars=gat.out, loglik=lik,
              list_priors=list_priors, list_posteriors=list_posteriors,
              list_density=list_density))
}


optimize_block <- function(node, treestr, par.list, Y=NULL, X, wts,
                           block_type=c("gaussian", "binomial"),
                           hme_type=c("hme", "hmre"), ...)
{
  block_type <- match.arg(block_type)
  hme_type <- match.arg(hme_type)
  theta <- par.list[[node]]
  # log-likelihood and gradient functionals
  D <- gradient(block_type) # needs to be function that returns a function
  Q <- Q(block_type)
  if (block_type == "binomial") {
    out <- 
    out <- optim(theta, fn=Q, gr=D, Z=X, joint_post=wts$joint_post, branch=wts$branch,
                 method="Nelder-Mead", control=list(fnscale=-1), ...)
    
  } else if (block_type == "gaussian") {
    exp.idx <- expert_index(hme_type, treestr)[node]
    #wt <- as.array(H[, exp.idx])
    out <- optim(theta, fn=Q, gr=D, Y=Y, X=X, wt=wts, method="Nelder-Mead", 
                 control=list(fnscale=-1), ...)
    
  }
  return(out)
}


# function factory
gradient <- function(grad_type=c("binomial", "gaussian"))
{
  grad_type <- match.arg(grad_type)
  if (grad_type == "binomial") {
    #f_ <- function(w, Z, wt_left, wt_right)
    #{
    #  eta <- exp(Z %*% w)
    #  g <- eta / (1 + eta)
    #  leftsplit <- wt_left * (1 - g)
    #  rightsplit <- wt_right * g
    #  contributions <- sweep(Z, 1, as.array(leftsplit - rightsplit), `*`)
    #  return(colSums(contributions))
    #}
    f_ <- function(w, Z, branch, joint_post)
    {
      eta <- exp(Z %*% w)
      g <- eta / (1 + eta)
      contributions <- sweep(Z, 1, as.array(joint_post * branch * (1 - g)), `*`)
      return(colSums(contributions))
    }
  } else if (grad_type=="gaussian") {
    f_ <- function(parm, Y, X, wt)
    {
      beta <- parm[1:(length(parm) - 1)]
      variance <- exp(parm[length(parm)])
      eps <- Y - (X %*% beta)
      g <- wt * (-0.5 * ((1 / variance) - (eps / variance)**2))
      beta_contrib <- sweep(X, 1, (wt * eps) / variance, `*`)
      return(colSums(cbind(beta_contrib, g)))
    }
  }
  return(f_)
}


# likelihood factory
Q <- function(like_type=c("binomial", "gaussian"))
{
  like_type <- match.arg(like_type)
  if (like_type == "binomial") {
    # f_ <- function(w, Z, wt_left, wt_right)
    # {
    #   eta <- Z %*% w
    #   l1p_eta <- log1p(exp(eta))
    #   leftsplit <- rowSums(H[, expert_list[[1]], drop=F]) * (eta - l1p_eta)
    #   leftsplit <- wt_left * (eta - l1p_eta)
    #   rightsplit <- wt_right * l1p_eta
    #   return(sum(leftsplit - rightsplit))
    # }
    f_ <- function(w, Z, branch, joint_post)
    {
      eta <- Z %*% w
      l1p_eta <- log1p(exp(eta))
      return(sum(joint_post * branch * (eta - l1p_eta)))
    }
  } else if (like_type == "gaussian") {
    f_ <- function(parm, Y, X, wt)
    {
      beta <- parm[1:(length(parm) - 1)]
      variance <- exp(parm[length(parm)])
      epsq <- (Y - (X %*% beta))**2
      contributions <- -0.5 * log(2 * pi * variance) - (epsq / (2 * variance))
      return(sum(wt * contributions))
    }
  }
  return(f_)
}


hme <- function(formula, tree, hme_type=c("hme", "hmre"), data, maxiter=100,
                tolerance=1e-3)
{
  cl <- match.call()
  require(Formula)
  stopifnot(!missing(data))
  
  tree <- sort(tree)
  hme_type <- match.arg(hme_type)
  mf <- match.call(expand.dots = FALSE)
  #m <- match(c("formula", "data"), names(mf), 0)
  form <- Formula::as.Formula(formula)
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  Z <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  
  tree <- sort(tree)
  expert.nodes <- tree[unlist(is_terminal(tree, tree))]
  gate.nodes <- setdiff(tree, expert.nodes)
  
  # This needs to be generalized for HMRE
  expert.index <- expert_index(hme_type, tree)
  
  gate.pars <- napply(gate.nodes, function(x) runif(ncol(Z), -2, 2))
  expert.pars <- napply(expert.nodes, function(x) c(runif(ncol(X) , -2, 2),
                                                    runif(1, 1, 5)))
  
  list_nodes <- napply(gate.nodes, par_to_gate_paths, "binomial", gate.pars, Z)
  list_experts <- napply(expert.nodes, par_to_expert_dens, "gaussian",
                          expert.pars, Y, X)

  
  logL <- matrix(NA_real_, nrow=maxiter + 1)
  parM <- matrix(NA_real_, ncol=length(unlist(c(gate.pars, expert.pars))),
                 nrow=maxiter + 1)
  # What is the likelihood
  logL[1,] <- sum(posterior_weights("0", tree, list_nodes, list_experts))
  parM[1, ] <- unlist(c(gate.pars, expert.pars))
  for (ii in seq_len(maxiter)) {
    mstep <- m_step(tree, hme_type, Y=Y, X=X, Z=Z, exp.pars=expert.pars, gat.pars=gate.pars)
    expert.pars <- mstep[["exp.pars"]]
    gate.pars <- mstep[["gat.pars"]]
    logL[ii + 1, ] <- mstep[["loglik"]]
    parM[ii + 1, ] <- unlist(c(gate.pars, expert.pars))
  }
  return(list(gate.pars=gate.pars,
              expert.pars=expert.pars,
              logL=logL,
              parM=parM,
              list_priors=mstep$list_priors,
              list_posteriors=mstep$list_posteriors,
              list_density=mstep$list_density))
}


data(iris)

tree <- c("0",
         "0.1", "0.2",
         "0.1.1", "0.1.2")
debugonce(hme)
tst <- hme("Sepal.Width ~ Petal.Width | -1 + Species + Petal.Length + Sepal.Length",
           tree=tree, data=iris)


cols <- c("blue", "orange", "green")
with(iris, plot(Petal.Width, Sepal.Width, col=cols[as.integer(iris$Species)]))
for (e in tst$expert.pars) {
  abline(e[1], e[2])
}




