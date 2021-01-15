"standard hme approach where each node is updated independently"
m_step <- function(tree, hme_type, expert_type, Y, X, Z, exp.pars, gat.pars, root_prior)
{
  gat.out <- lapply(gat.pars, function(x) list())
  exp.out <- lapply(exp.pars, function(x) list())
  
  list_priors <- napply(names(gat.pars), par_to_gate_paths, gat.pars, Z)
  list_density <- napply(names(exp.pars), par_to_expert_dens, expert_type,
                         exp.pars, Y, X)
  list_posteriors <- napply(names(gat.pars), posterior_weights,
                            tree, list_priors, list_density)
  
  for (e in names(exp.pars)) {
    jpw <- joint_posterior_weight(e, list_posteriors, root_prior)
    opt_blk <- optimize_block(e, tree, exp.pars, Y, X, jpw, expert_type, hme_type)
    exp.out[[e]] <- opt_blk$par
  }
  
  for (g in names(gat.pars)) {
    "test for binomial vs multinomial here"
    nchds <- length(unlist(children(g, tree)))
    node_type <- "multinomial"
    if (nchds == 2)
      node_type <- "binomial"
    # joint posterior weights
    wt <- joint_posterior_weight(g, list_posteriors, root_prior)
    # posterior branches
    target <- list_posteriors[[g]]
    wts=list(branch=target, joint_post=wt)
    opt_blk <- optimize_block(g, tree, gat.pars, NULL, Z, wts, node_type, hme_type)
    gat.out[[g]] <- opt_blk$par
  }
  
  list_priors <- napply(names(gat.pars), par_to_gate_paths, gat.pars, Z)
  list_density <- napply(names(exp.pars), par_to_expert_dens, expert_type,
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
                           block_type=c("binomial", "multinomial", "gaussian",
                                        "bernoulli"),
                           hme_type=c("hme", "hmre"), ...)
{
  block_type <- match.arg(block_type)
  hme_type <- match.arg(hme_type)
  theta <- par.list[[node]]
  # log-likelihood and gradient functionals
  D <- gradient(block_type) # needs to be function that returns a function
  Q <- Q(block_type)
  if (block_type == "binomial") {
    
    out <- irls_logit(X, wts$branch[, 1], wts$joint_post, maxit=10)
    
  } else if (block_type == "multinomial") {
    # H is a matrix of posterior branches
    out <- irls_multinomial(X, wts$branch, wts$joint_post, maxit=10, tol=1e-08)
    
  } else if (block_type == "gaussian") {
    exp.idx <- expert_index(hme_type, treestr)[node]
    #wt <- as.array(H[, exp.idx])
    out <- optim(theta, fn=Q, gr=D, Y=Y, X=X, wt=wts, method="Nelder-Mead", 
                 control=list(fnscale=-1), ...)
    
  } else if (block_type == "bernoulli") {
    out <- glm.fit(x=X, y=Y, weights=wts, family=binomial(link="logit"))
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
      beta <- parm[-length(parm)]
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



logistic_score <- function(node, list_post, list_priors, Z)
{
  # joint posterior
  H <- joint_posterior_weight(node=node, lp=list_post, rp=1)
  # prior weight for gating node
  g <- list_priors[[node]]
  
  splits <- ncol(g) - 1L
  omega_score <- vector("list", length(splits))
  for (s in seq_len(splits)) {
    # sweep out gating varibles Z in direction 1
    omega_score[[s]]  <- sweep(Z, 1, as.array(H * (1 - g[, s])), `*`)
  }
  # sweep out gating varibles Z in direction 1
  # omega_score <- sweep(Z, 1, as.array(H * (1 - g[, 1])), `*`)
  return(c(omega_score))
}


logistic_hessian <- function(node, list_post, list_priors, Z)
{
  
  # joint posterior
  H <- joint_posterior_weight(node=node, lp=list_post, rp=1)
  # prior weight for gating node
  g <- list_priors[[node]]
  
  nrz <- nrow(Z)
  ncz <- ncol(Z)
  ncg <- ncol(g)
  
  # Outer hessian
  ZZ <- apply(Z, 1, function(x) x %*% t(x))
  dim(ZZ) <- c(ncz, ncz, nrz)
  
  gamma_t <- function(z)
  {
    o <- z %*% t(z)
    diag(o) <- -z * (1 - z)
    return(o)
  }
  Gamma <- apply(g[, -ncg, drop=FALSE], 1, gamma_t)
  dim(Gamma) <- c(ncg - 1L, ncg - 1L, nrow(g))
  
  d <- (ncg - 1) * ncz
  hess <- array(dim=c(d, d, nrow(Z)))
  for (i in seq_len(nrz)) {
    hess[,,i] <- H[i] * Gamma[,,i] %x% ZZ[,,i]
  }
  return(hess)
}


gaussian_score <- function(node, expert.pars, list_post, list_priors, Y, X)
{
  # extract the variables for the expert
  pars <- expert.pars[[node]]
  np <- length(pars)
  betas <- pars[-np]
  variance <- exp(pars[np])
  
  # joint posterior probability and residuals
  H <- joint_posterior_weight(node, lp=list_post, rp=1)
  eps <- Y - expert_pred(pars, X, "gaussian")
  
  beta_score <- sweep(X, 1, as.array(H * (eps / variance)), `*`)
  var_score <- -0.5 * H * ((eps**2 / variance) - 1)
  gaussian_score <- cbind(beta_score, var_score)
  colnames(gaussian_score) <- c(colnames(beta_score), "variance")
  return(gaussian_score)
}


gaussian_hessian <- function(node, expert.pars, list_post, list_priors, Y, X)
{
  # extract the variables for the expert
  pars <- expert.pars[[node]]
  np <- length(pars)
  betas <- pars[-np]
  variance <- exp(pars[np])
  
  # joint posterior probability
  H <- joint_posterior_weight(node, lp=list_post, rp=1)
  yhat <- expert_pred(pars, X, "gaussian")
  eps <- Y - yhat
  
  out_array <- array(NA_real_, dim=c(length(pars), length(pars), length(Y)))
  for (ii in seq_along(Y)) {
    XX <- X[ii, ] %o% X[ii, ]
    Xe <- X[ii, ] * eps[ii]
    VV <- -0.5 * eps[ii]**2
    out_array[,,ii] <- (H[ii] / variance) * cbind(rbind(XX, Xe), c(t(Xe), VV))
  }
  return(out_array)
}


# Take a list of square matrices and make them block diagonal
block_diag <- function(lst)
{
  if (!inherits(lst, "list")) {
    stop("lst must be a list of matrices.")
  }
  clss <- sapply(lst, inherits, "matrix")
  if (!all(clss)) {
    stop("All elements of `lst` must be of class `matrix`.")
  }
  nc <- sapply(lst, ncol)
  nr <- sapply(lst, nrow)
  if (any(nc != nr)) {
    stop("Elements of `lst` are not all square matrices")
  }
  out <- matrix(0, nrow=sum(nr), ncol=sum(nc))
  idx <- 0
  for (ii in seq_along(nr)) {
    ed <- cumsum(nr)[ii]
    st <- idx + 1
    out[st:ed, st:ed] <- lst[[ii]]
    idx <- ed
  }
  return(out)
}


sandwich_vcov <- function(gte.nms, exp.nms, lp, ln, exp.pars, Y, X, Z, N)
{
  gate_scores <- napply(gte.nms, logistic_score, lp, ln, Z)
  
  expt_scores <- napply(exp.nms, gaussian_score, exp.pars, lp, ln, Y, X)
  
  # Create the full score vector of theta = (omega + beta)
  scores <- do.call(cbind, c(unlist(gate_scores, recursive=FALSE), expt_scores))
  nc <- ncol(scores)
  rm(gate_scores, expt_scores)
  
  # Sum outer product of each input pattern
  scores <- apply(scores, 1, function(x) x %*% t(x))
  dim(scores) <- c(nc, nc, N)
  scores <- rowSums(scores, dims=2)
  
  gate_hess <- napply(gte.nms, logistic_hessian, lp, ln, Z)
  gate_hess <- lapply(gate_hess, rowSums, dims=2)
  
  expt_hess <- napply(exp.nms, gaussian_hessian, exp.pars, lp, ln, Y, X)
  expt_hess <- lapply(expt_hess, rowSums, dims=2)
  
  H <- block_diag(c(gate_hess, expt_hess))
  
  vcv <- solve(H) %*% scores %*% solve(H)
  
  return(list(hessian=H, OPG=scores, sandwich=vcv))
}


