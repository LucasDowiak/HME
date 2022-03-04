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
  lik <- log_likelihood(tree, list_priors, list_density)
  
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
  # needs to be function that returns a function
  D <- gradient(block_type) 
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
  return(f_)
}


Omega_0 <- function(node, list_priors, GD)
{
  g <- list_priors[[node]]
  # For the last split we have set omega^{J} = 0 so no score is needed for it
  splits <- ncol(g) - 1L
  split_nms <- paste(node, seq_len(splits), sep="|")
  
  # Find experts that are descendants from the gating node
  all_experts <- names(GD)
  tree <- c(names(list_priors), names(GD))
  node_progeny <- unlist(progeny(node, tree))
  node_experts <- intersect(all_experts, node_progeny)
  
  omega_0 <- vector("list", length(split_nms))
  names(omega_0) <- split_nms
  for (i in seq_len(splits)) {
    # split between explicit and implicit experts
    expl_experts <- node_experts[grepl(paste(node, i, sep="."), node_experts)]
    impl_experts <- setdiff(node_experts, expl_experts)
    
    # value for explicit experts
    expl_omega_0 <- (1 - g[, i]) * Reduce(`+`, GD[expl_experts])
    
    # All other splits need to be considered
    impl_omega_0 <- vector("list", ncol(g))
    for (t in seq_len(ncol(g))) {
      if (i == t) {
        # Zero will not affect the sum later on
        impl_omega_0[[t]] <- rep(0, length(GD[[1]]))
      } else {
        impl_experts2 <- impl_experts[grepl(paste(node, t, sep="."), impl_experts)]
        impl_omega_0[[t]] <- g[, t] * Reduce(`+`, GD[impl_experts2])
      }
    }
    impl_omega_0 <- Reduce(`+`, impl_omega_0)
    omega_0[[i]] <- expl_omega_0 - impl_omega_0
  }
  return(omega_0)
}


# TODO: 2) Re-organize to be called from inside sandwich_cov
logistic_score <- function(node, list_priors, list_density, Z)
{
  # For the last split we have set omega^{J} = 0 so no score is needed for it
  splits <- ncol(list_priors[[node]]) - 1L
  split_nms <- paste(node, seq_len(splits), sep="|")
  
  # Find experts that are descendants from node
  all_experts <- names(list_density)
  
  # Each experts contribution to the likelihood
  GD <- expert_lik_contr(all_experts, list_density, list_priors)
  denom <- Reduce(`+`, GD)
  
  omega_0 <- Omega_0(node, list_priors, GD)
  out <- lapply(omega_0, function(x) sweep(Z, 1, x / denom, FUN = `*`))
  names(out) <- split_nms
  return(out)
}


# TODO: 2) Re-organize to be called from inside sandwich_cov
gaussian_score <- function(expert, expert.pars, list_priors, list_density, Y, X)
{
  # extract the variables for the expert and calculate error terms
  pars <- expert.pars[[expert]]
  np <- length(pars)
  betas <- pars[-np]
  variance <- exp(pars[np])
  eps <- Y - expert_pred(pars, X, "gaussian")
  
  # prior weights for descendants experts
  all_experts <- names(list_density)
  
  GD <- expert_lik_contr(all_experts, list_density, list_priors)
  expert_lik_prop <- GD[[expert]] / Reduce(`+`, GD)
  
  beta_score <- sweep(X, 1, as.array(expert_lik_prop * (eps / variance)), `*`)
  var_score <- 0.5 * expert_lik_prop * ( (eps**2 / variance) - 1 )
  
  gaussian_score <- cbind(beta_score, var_score)
  colnames(gaussian_score) <- c(colnames(beta_score), "variance")
  return(gaussian_score)
}


# TODO: 2) Re-organize to be called from inside sandwich_cov
cross_logistic_hessian <- function(node1, node2, list_priors, list_density, Z)
{
  # node splits
  g1 <- list_priors[[node1]]
  g2 <- list_priors[[node2]]
  # For the last split we have set omega^{J} = 0 so no hessian is needed for it
  splits1 <- ncol(g1) - 1L
  splits2 <- ncol(g2) - 1L
  
  # expert names
  all_experts <- names(list_density)
  GD <- expert_lik_contr(all_experts, list_density, list_priors)
  denom <- Reduce(`+`, GD)
  omega_01 <- Omega_0(node1, list_priors, GD)
  omega_02 <- Omega_0(node2, list_priors, GD)
  
  # Find experts that are descendants from nodes 1 and 2
  tree <- c(names(list_priors), names(GD))
  node1_progeny <- unlist(progeny(node1, tree))
  node1_experts <- intersect(all_experts, node1_progeny)
  node2_progeny <- unlist(progeny(node2, tree))
  node2_experts <- intersect(all_experts, node2_progeny)
  
  shared_experts <- intersect(node1_experts, node2_experts)
  has_shared_experts <- length(shared_experts) > 0
  
  # Each interaction between the splits in node1 and node2 has a ZxZxT array as an output
  cross_split_names <- apply(expand.grid(paste(node1, seq_len(splits1), sep="|"),
                                         paste(node2, seq_len(splits2), sep="|"),
                                         stringsAsFactors = FALSE),
                             1, paste, collapse="-")
  cross_split_names <- matrix(cross_split_names, nrow=splits1)
  unq_cross_split_names <- cross_split_names[lower.tri(cross_split_names, diag=TRUE)]
  out_arrays <- napply(unq_cross_split_names, function(x) array(0, dim=c(ncol(Z), ncol(Z))))
  
  for (csn in unq_cross_split_names) {
    nodesplits <- unlist(strsplit(csn, split="-"))
    ns1 <- unlist(strsplit(nodesplits[1], "\\|"))
    ns2 <- unlist(strsplit(nodesplits[2], "\\|"))
    split_nm_perd_1 <- paste(ns1[1], ns1[2], sep=".")
    split_nm_pipe_1 <- paste(ns1[1], ns1[2], sep="|")
    split_nm_perd_2 <- paste(ns2[1], ns2[2], sep=".")
    split_nm_pipe_2 <- paste(ns2[1], ns2[2], sep="|")
    
    for (t in seq_len(nrow(Z))) {
      ZZ <- Z[t, ] %o% Z[t, ]
      square_term <- denom[t]**2 * omega_01[[split_nm_pipe_1]][t] * omega_02[[split_nm_pipe_2]][t]
      partial_term <- 0
      
      if (has_shared_experts) {
        omega_1 <- omega_2 <- 0
      
        for (expert in shared_experts) {
          is_explicit_1i <- as.integer(grepl(split_nm_perd_1, expert))
          is_explicit_2j <- as.integer(grepl(split_nm_perd_2, expert))
          omega_1 <- omega_1 + prod(is_explicit_1i - g1[t, as.integer(ns1[2])],
                                    is_explicit_2j - g2[t, as.integer(ns2[2])],
                                    GD[[expert]][t])
          if (ns1[1] == ns2[1]) {
            same_split <- as.integer(ns1[2] == ns2[2])
            omega_2 <- omega_2 - g1[t, as.integer(ns1[2])] * (same_split - g2[t, as.integer(ns2[2])]) * GD[[expert]][t]
          }
        }
        partial_term <- omega_1 + omega_2
      }
      out_arrays[[csn]] <- out_arrays[[csn]] + (denom[t] * partial_term - square_term) * ZZ
    }
  }
  # Divide by sample size to take average
  # out_arrays <- lapply(out_arrays, function(x) (1 / nrow(Z)) * x)
  
  # fit together the output of out_arrays into the hessian
  # this will be a splits1 * splits2 * length(gate.nodes)
  colms <- vector("list", ncol(cross_split_names))
  
  for (cc in seq_len(ncol(cross_split_names))) {
    
    lst_cc <- vector("list", nrow(cross_split_names))
    
    for (rr in seq_len(nrow(cross_split_names))) {
      
      nm <- cross_split_names[rr, cc]
      
      if (nm %in% unq_cross_split_names) {
        
        lst_cc[[rr]] <- out_arrays[[nm]]
        
      } else {
        
        nmpairs <- unlist(strsplit(nm, "-"))
        
        newnm <- paste(c(nmpairs[2], nmpairs[1]), collapse="-")
        
        lst_cc[[rr]] <- out_arrays[[newnm]]
      }
    }
    colms[[cc]] <- do.call(rbind, lst_cc)
  }
  out <- do.call(cbind, colms)
  return(out)
}


# TODO: 2) Re-organize to be called from inside sandwich_cov
gaussian_hessian <- function(expert, expert.pars, list_priors, list_density, Y, X)
{
  # extract the variables for the expert and calculate error terms
  pars <- expert.pars[[expert]]
  np <- length(pars)
  betas <- pars[-np]
  variance <- exp(pars[np])
  eps <- Y - expert_pred(pars, X, "gaussian")
  std_err <- eps / variance
  
  all_experts <- names(list_density)
  # ---------------- Just pass this in form sandwich cov function
  GD <- expert_lik_contr(all_experts, list_density, list_priors)
  expert_lik_prop <- GD[[expert]] / Reduce(`+`, GD)
  
  out_array <- array(NA_real_, dim=c(length(pars), length(pars), length(Y)))
  for (t in seq_along(Y)) {
    Xt <- X[t, ]
    XX <- Xt %o% Xt
    XX <- expert_lik_prop[t] * ( std_err[t]**2 * (1 - expert_lik_prop[t]) - (1 / variance) ) * XX
    XV <- 0.5 * expert_lik_prop[t] * std_err[t] * ( (eps[t]**2 / variance - 1) * (1 - expert_lik_prop[t]) - 2 ) * Xt
    VV <- 0.25 * expert_lik_prop[t] * ( (eps[t]**2 / variance - 1)**2 * (1 - expert_lik_prop[t]) - (4 * eps[t]**2 / variance) )
    out_array[,,t] <- cbind(rbind(XX, XV), c(t(XV), VV))
  }
  out_array <- rowSums(out_array, dims=2L)
  return(out_array)
}


# TODO: 2) Re-organize to be called from inside sandwich_cov
cross_gaussian_hessian <- function(expert1, expert2, expert.pars, list_priors, list_density, Y, X)
{
  # extract the variables for the expert and calculate error terms
  pars1 <- expert.pars[[expert1]]
  np <- length(pars1)
  betas1 <- pars1[-np]
  variance1 <- exp(pars1[np])
  eps1 <- Y - expert_pred(pars1, X, "gaussian")
  
  pars2 <- expert.pars[[expert2]]
  betas2 <- pars2[-np]
  variance2 <- exp(pars2[np])
  eps2 <- Y - expert_pred(pars2, X, "gaussian")
  
  # prior weights for descendants experts
  all_experts <- names(list_density)

  # ---------------- Just pass this in form sandwich cov function
  GD <- expert_lik_contr(all_experts, list_density, list_priors)
  denom <- Reduce(`+`, GD)
  
  expert_lik_prop1 <- GD[[expert1]] / denom
  expert_lik_prop2 <- GD[[expert2]] / denom
  std_err1 <- eps1 / variance1
  std_err2 <- eps2 / variance2
  
  out_array <- array(NA_real_, dim=c(length(pars1), length(pars1), length(Y)))
  for (t in seq_along(Y)) {
    Xt <- X[t, ]
    XX <- Xt %o% Xt
    XX <- -expert_lik_prop1[t] * expert_lik_prop2[t] * std_err1[t] * std_err2[t] * XX
    XV <- -0.5 * expert_lik_prop1[t] * expert_lik_prop2[t] * std_err1[t] * (eps2[t]**2 / variance2 - 1) * Xt
    VV <- -0.25 * expert_lik_prop1[t] * expert_lik_prop2[t] * (eps1[t]**2 / variance1 - 1) * (eps2[t]**2 / variance2 - 1)
    out_array[,,t] <- cbind(rbind(XX, XV), c(t(XV), VV))
  }
  out_array <- rowSums(out_array, dim=2L)
  return(out_array)
}


# TODO: 2) Re-organize to be called from inside sandwich_cov
cross_gaussian_multinomial_hessian <- function(node, expert, expert.pars, list_priors,
                                               list_density, Y, X, Z)
{
  
  # node splits
  g <- list_priors[[node]]
  # For the last split we have set omega^{J} = 0 so no hessian is needed for it
  splits <- ncol(g) - 1L
  
  # extract the variables for the expert and calculate error terms
  pars <- expert.pars[[expert]]
  np <- length(pars)
  betas <- pars[-np]
  variance <- exp(pars[np])
  eps <- Y - expert_pred(pars, X, "gaussian")
  std_err <- eps / variance
  
  all_experts <- names(list_density)
  
  # Need to find out if m \in M(node)
  tree <- c(names(list_priors), all_experts)
  node_progeny <- unlist(progeny(node, tree))
  node_experts <- intersect(all_experts, node_progeny)
  is_child_expert = as.integer(expert %in% node_experts)
  
  GD <- expert_lik_contr(all_experts, list_density, list_priors)
  denom <- Reduce(`+`, GD)
  
  expert_lik_prop <- GD[[expert]] / denom
  omega_0 <- Omega_0(node=node, list_priors=list_priors, GD=GD)
  
  out_array <- lapply(1:splits, function(x) array(NA_real_, dim=c(length(pars), ncol(Z), length(Y))))
  names(out_array) <- paste(node, 1:splits, sep="|")
  
  for (t in seq_along(Y)) {
    
    XZ <- X[t, ] %o% Z[t, ]
    
    for (i in seq_len(splits)) {
      split_nm_pipe <- paste(node, i, sep="|")
      split_nm_perd <- paste(node, i, sep=".")
      # explicit or implicit relationship
      is_explicit <- as.integer(grepl(split_nm_perd, expert))
      
      XZ_tmp <- expert_lik_prop[t] * std_err[t] * ( (is_explicit - g[t, i]) * is_child_expert - omega_0[[split_nm_pipe]][t] / denom[t] ) * XZ
      
      ZV_tmp <- 0.5 * expert_lik_prop[t] * (eps[t]**2 / variance - 1) * ( (is_explicit - g[t, i]) * is_child_expert - omega_0[[split_nm_pipe]][t] / denom[t] ) * Z[t, ]
      
      out_array[[i]][,,t] <- rbind(XZ_tmp, ZV_tmp)
    }
  }
  out_array <- lapply(out_array, rowSums, dim=2)
  return(do.call(cbind, out_array))
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


# Calculates each experts contribution to the likelihood (e.g: f_{m ,t})
expert_lik_contr <- function(exp.nms, ld, ln)
{
  # prior weights for descendants experts
  Pi <- napply(exp.nms, prior_weights, ln=ln)
  # prior *  P^{k}
  return(napply(exp.nms, function(x) Pi[[x]] * ld[[x]]))
}


# TODO: Incorporate all new hessian values into the 
sandwich_vcov <- function(gte.nms, exp.nms, ld, ln, gate.pars, exp.pars, Y, X, Z, N)
{
  # Calculate each experts contribution to the likelihood value
  GD <- expert_lik_contr(exp.nms, ld, ln)
  
  #--------------------------------------------------------------------
  gate_scores <- napply(gte.nms, logistic_score, ln, ld, Z)
  
  expt_scores <- napply(exp.nms, gaussian_score, exp.pars, ln, ld, Y, X)
  
  # Create the full score vector of theta = (omega + beta)
  scores <- do.call(cbind, c(unlist(gate_scores, recursive=FALSE), expt_scores))
  nc <- ncol(scores)
  rm(gate_scores, expt_scores)
  
  # Sum outer product of each input pattern
  scores <- apply(scores, 1, function(x) x %*% t(x))
  dim(scores) <- c(nc, nc, N)
  scores <- rowSums(scores, dims=2)
  
  # Block diagonal each node's personally hessian: partial^2 l_t / partial omega^2
  gate_hess <- napply(gte.nms, function(x) cross_logistic_hessian(x, x, ln, ld, Z))
  expt_hess <- napply(exp.nms, gaussian_hessian, exp.pars, ln, ld, Y, X)
  H <- block_diag(c(gate_hess, expt_hess))
  
  # IX <- list mapping node names to parameter locations: f("0.1") -> [12, 13, 14, 15]
  make_ind_vec <- function(nme, lst)
  {
    out <- c()
    n <- length(lst[[nme]])
    for (i in seq_len(n)) {
      o <- lst[[nme]][[i]]
      o[] <- 1L
      names(o) <- rep(nme, length(o))
      out <- c(out, o)
    }
    return(out)
  }
  IX <- unlist(c(lapply(gte.nms,  make_ind_vec, gate.pars),
                 lapply(exp.nms, make_ind_vec, exp.pars)))
  IX <- cumsum(IX)
  IX <- split(IX, as.factor(names(IX)))
  
  # find unique cross combos
  cross_exp <- apply(combn(exp.nms, 2, simplify = T), 2, paste, collapse="-")
  for (combo in cross_exp) {
    combo <- unlist(strsplit(combo, split="-"))
    tmp_hess <- cross_gaussian_hessian(combo[1], combo[2], exp.pars, ln, ld, Y, X)
    H[IX[[combo[1]]], IX[[combo[2]]]] <- tmp_hess
    H[IX[[combo[2]]], IX[[combo[1]]]] <- tmp_hess
  }
  
  # find unique cross combos
  if (length(gte.nms) > 1) {
    cross_gates <- apply(combn(gte.nms, 2, simplify = T), 2, paste, collapse="-")
    for (combo in cross_gates) {
      combo <- unlist(strsplit(combo, split="-"))
      tmp_hess <- cross_logistic_hessian(combo[1], combo[2], ln, ld, Z)
      H[IX[[combo[1]]], IX[[combo[2]]]] <- tmp_hess
      H[IX[[combo[2]]], IX[[combo[1]]]] <- tmp_hess
    } 
  }
  
  # find unique cross combos
  cross_gate_experts <- apply(expand.grid(gte.nms, exp.nms), 1, paste, collapse="-")
  for (combo in cross_gate_experts) {
    combo <- unlist(strsplit(combo, split="-"))
    tmp_hess <- cross_gaussian_multinomial_hessian(combo[1], combo[2], exp.pars, ln, ld, Y, X, Z)
    H[IX[[combo[1]]], IX[[combo[2]]]] <- t(tmp_hess)
    H[IX[[combo[2]]], IX[[combo[1]]]] <- tmp_hess
  }
  
  # The sandwich
  vcv <- solve(H) %*% scores %*% solve(H)
  
  return(list(hessian=H, OPG=scores, sandwich=vcv))
}


calculate_sandwich_vcov <- function(mod)
{
  return(sandwich_vcov(mod$gate.nodes, mod$expert.nms, mod$list_density,
                       mod$list_priors, mod$gate.pars, mod$expert.pars, mod$Y,
                       mod$X, mod$Z, mod$N))
}


