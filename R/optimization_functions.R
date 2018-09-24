"standard hme approach where each node is updated independently"
m_step <- function(tree, hme_type, Y, X, Z, exp.pars, gat.pars, root_prior)
{
  gat.out <- lapply(gat.pars, function(x) list())
  exp.out <- lapply(exp.pars, function(x) list())
  
  list_priors <- napply(names(gat.pars), par_to_gate_paths, gat.pars, Z)
  list_density <- napply(names(exp.pars), par_to_expert_dens, "gaussian",
                         exp.pars, Y, X)
  list_posteriors <- napply(names(gat.pars), posterior_weights,
                            tree, list_priors, list_density)
  
  for (e in names(exp.pars)) {
    jpw <- joint_posterior_weight(e, tree, list_posteriors, root_prior)
    opt_blk <- optimize_block(e, tree, exp.pars, Y, X, jpw, "gaussian", hme_type)
    exp.out[[e]] <- opt_blk$par
  }
  
  for (g in names(gat.pars)) {
    "test for binomial vs multinomial here"
    nchds <- length(unlist(children(g, tree)))
    node_type <- "multinomial"
    if (nchds == 2)
      node_type <- "binomial"
    # joint posterior weights
    wt <- joint_posterior_weight(g, tree, list_posteriors, root_prior)
    # posterior branches
    target <- list_posteriors[[g]]
    wts=list(branch=target, joint_post=wt)
    opt_blk <- optimize_block(g, tree, gat.pars, NULL, Z, wts, node_type, hme_type)
    gat.out[[g]] <- opt_blk$par
  }
  
  list_priors <- napply(names(gat.pars), par_to_gate_paths, gat.pars, Z)
  list_density <- napply(names(exp.pars), par_to_expert_dens, "gaussian",
                         exp.pars, Y, X)
  list_posteriors <- napply(names(gat.pars), posterior_weights,
                            tree, list_priors, list_density)
  lik <- log_likelihood(tree, list_priors, list_posteriors, list_density,
                        root_prior)
  
  return(list(exp.pars=exp.out, gat.pars=gat.out, loglik=lik,
              list_priors=list_priors, list_posteriors=list_posteriors,
              list_density=list_density))
}


optimize_block <- function(node, treestr, par.list, Y=NULL, X, wts,
                           block_type=c("binomial", "multinomial", "gaussian"),
                           hme_type=c("hme", "hmre"), ...)
{
  block_type <- match.arg(block_type)
  hme_type <- match.arg(hme_type)
  theta <- par.list[[node]]
  # log-likelihood and gradient functionals
  D <- gradient(block_type) # needs to be function that returns a function
  Q <- Q(block_type)
  if (block_type == "binomial") {
    
    out <- irls_logit(X, wts$branch[, 1], wts$joint_post, maxit=3)
    
  } else if (block_type == "multinomial") {
    # H is a matrix of posterior branches
    out <- irls_multinomial(X, wts$branch, wts$joint_post, maxit=3, tol=1e-08)
    
  } else if (block_type == "gaussian") {
    exp.idx <- expert_index(hme_type, treestr)[node]
    #wt <- as.array(H[, exp.idx])
    out <- optim(theta, fn=Q, gr=D, Y=Y, X=X, wt=wts, method="Nelder-Mead", 
                 control=list(fnscale=-1), ...)
    
  }
  return(out)
}


irls_logit <- function(X, H, w=NULL, maxit=25, tol=1e-08)
{
  family <- binomial
  WNULL <- is.null(w)
  s <- t <- 0
  QR <- qr(X)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  for(iter in 1:maxit) {
    g <- family()$linkinv(t)
    gprime <- family()$mu.eta(t)
    z <- t + (H - g) / gprime
    if (WNULL) {
      A <- 1
    } else {
      A <- w
    }
    W <- as.vector((A / family()$variance(g)) * gprime^2)
    wmin <- min(W)
    if(wmin < sqrt(.Machine$double.eps))
      warning("Tiny weights encountered")
    s_old <- s
    C <- chol(crossprod(Q, W*Q))
    s <- forwardsolve(t(C), crossprod(Q,W*z))
    s <- backsolve(C,s)
    t <- Q %*% s
    if(sqrt(crossprod(s - s_old)) < tol)
      break
  }
  x <- backsolve(R, crossprod(Q,t))
  list(par=c(x), iterations=iter)
}


irls_multinomial <- function(X, H, w=NULL, maxit=25, tol=1e-08)
{
  nc <- ncol(H)
  lst <- vector("list", nc-1)
  for (cc in seq_len(nc-1)) {
    bM <- H[, c(cc, nc)]
    bM <- sweep(bM, 1, rowSums(bM), `/`)
    # if X is list X[[g]] here
    #bM <- H[, cc] / H[, nc]
    lst[[cc]] <- irls_logit(X, bM[,1], w, maxit=maxit, tol=tol)
  }
  grab <- function(x) c(x$par)
  list(par=lapply(lst, grab), iterations=lapply(lst, `[[`, "iterations")) 
}


# function factory
gradient <- function(grad_type)
{
  if (grad_type == "gaussian") {
    f_ <- function(parm, Y, X, wt)
    {
      beta <- parm[1:(length(parm) - 1)]
      variance <- exp(parm[length(parm)])
      eps <- Y - (X %*% beta)
      g <- wt * (-0.5 * ((1 / variance) - (eps / variance)**2))
      beta_contrib <- sweep(X, 1, (wt * eps) / variance, `*`)
      return(colSums(cbind(beta_contrib, g)))
    }
  } else {
    return(NULL)
  }
  #  else if (grad_type == "binomial") {
  #  f_ <- function(w, Z, branch, joint_post)
  #  {
  #    eta <- exp(Z %*% w)
  #    g <- eta / (1 + eta)
  #    contributions <- sweep(Z, 1, as.array(joint_post * branch * (1 - g)), `*`)
  #    return(colSums(contributions))
  #  }
  #} else if (grad_type == "multinomial") {
  #  stop("Multinomial gradient is not yet written boi.")
  return(f_)
}


# likelihood factory
Q <- function(like_type)
{
  if (like_type == "gaussian") {
    f_ <- function(parm, Y, X, wt)
    {
      beta <- parm[1:(length(parm) - 1)]
      variance <- exp(parm[length(parm)])
      epsq <- (Y - (X %*% beta))**2
      contributions <- -0.5 * log(2 * pi * variance) - (epsq / (2 * variance))
      return(sum(wt * contributions))
    }
  } else {
    return(NULL)
  }
  # if (like_type == "binomial") {
  #  f_ <- function(w, Z, branch, joint_post)
  #  {
  #    eta <- Z %*% w
  #    l1p_eta <- log1p(exp(eta))
  #    return(sum(joint_post * branch * (eta - l1p_eta)))
  #  }
  #} else if (like_type == "multinomial") {
  #  stop("Multinomial likelihood function not written yet boi")
  # }

  return(f_)
}



multinomial_info_matrix <- function(node, treestr, gate.pars, ln, lp, Z,
                                    rp)
{
  gate.par <- gate.pars[[node]]
  if (is.list(gate.par))
    m <- length(gate.pars) + 1
  else
    m <- 2
  H <- joint_posterior_weight(node=node, treestr=treestr, lp=lp, rp=rp)
  h <- lp[[node]]
  g <- ln[[node]]
  
  I <- matrix(0, nrow=ncol(Z), ncol=ncol(Z))
  for (i in seq_len(nrow(Z))) {
   for (j in seq_len(ncol(h) - 1)) {
     wt <- H[i] * (h[i, j] - g[i, j])
     I <- I + wt**2 * (Z[i, ] %o% Z[i, ])
   }
  }
  
  G <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  for (j in seq_len(ncol(h) - 1)) {
    G <- G + sweep(Z, 1, H * (h[,j] - g[,j]), `*`) 
  }
  
  HS <- matrix(0, nrow=ncol(Z), ncol=ncol(Z))
  for (i in seq_len(nrow(Z))) {
    OZ <- Z[i, ] %o% Z[i, ]
    tmp <- matrix(0, nrow=ncol(Z), ncol=ncol(Z))
    for (j in seq_len(ncol(g) - 1)) {
      for (k in seq_len(ncol(g) - 1)) {
        tmp <- H[i] * ((j == k) - g[i, j]) * g[i, k] * OZ
      }
    }
    HS <- HS + tmp
  }
  OPG <- crossprod(G)
  sandwich <- solve(HS) %*% OPG %*% solve(HS)
  return(list(OPG=OPG, I=I, H=HS, sandwich=sandwich))
}


#debugonce(multinomial_info_matrix)
#IJJ <- multinomial_info_matrix("0", tree, tst$gate.pars, tst$list_priors,
#                               tst$list_posteriors, tst$Z)





