# for individual gating nodes, Subtract par means from path i
#    Equation \label{eq:gate_marginal_effect} for a single node
gate_margin <- function(node, gate.pars, ln)
{
  pars <- gate.pars[[node]]
  if (!is(pars, "list"))
    pars <- list(pars)
  # Set parameters branch J to 0
  pars[[length(pars) + 1]] <- rep(0, length(pars[[1]]))
  
  # Add up all the weighted parameter vectors
  nn <- length(pars)
  margs <- vector("list", nn)
  # Weight each by their path from `node`
  for (ii in seq_len(nn)) {
    margs[[ii]] <- sapply(pars[[ii]], `*`, ln[[node]][, ii])
  }
  par.means <- Reduce(`+`, margs)
  
  # Sweep out the differences
  f_ <- function(x) -sweep(par.means, 2, pars[[x]], FUN=`-`)
  p_diff <- lapply(seq_len(nn), f_)
  
  return(p_diff)
}


# Marginal Gating Values routed to each expert: \partial pi / \parial Z
#    Equation \label{eq:marginal_effect} for a single node
tree_expert_margin <- function(expert, gate.pars, ln, weighted=TRUE)
{
  gpp <- gate_path_product("0", expert, ln)
  # all nodes from the root node to the expert
  a <- unlist(ancestors(expert))[-1]
  # find the nodes
  abls <- unlist(lapply(a, all_but_last_split))
  # and the direction from each node
  lspl <- unlist(lapply(a, last_split))
  
  # add up all the parameter diffs
  treasure <- vector("list", length(a))
  for (nn in seq_along(a)) {
    int <- as.integer(lspl[nn])
    treasure[[nn]] <- gate_margin(abls[[nn]], gate.pars, ln)[[int]]
  }
  
  delta_m <- Reduce(`+`, treasure)
  
  # normalized by the gate path product
  if (weighted)
    delta_m <- sweep(delta_m, 1, gpp, `*`)
  
  return(delta_m)
}




margin_matrix <- function(experts, gate.pars, ln)
{
  # Avergage Marginal Effect
  f_ <-function(expert)
  {
    colMeans(tree_expert_margin(expert, gate.pars, ln, TRUE))
  }
  return(napply(experts, f_))
}


"the marginal effects for each expert should sum to zero"
# colMeans(tree_expert_margin("0.1", tst2$gate.pars, tst2$list_priors))
# colMeans(tree_expert_margin("0.2", tst2$gate.pars, tst2$list_priors))
# colMeans(tree_expert_margin("0.3", tst2$gate.pars, tst2$list_priors))

# colMeans(tree_expert_margin("0.1.1", tst$gate.pars, tst$list_priors))
# colMeans(tree_expert_margin("0.1.2", tst$gate.pars, tst$list_priors))
# colMeans(tree_expert_margin("0.2", tst$gate.pars, tst$list_priors))


fitted_expert <- function(expert, X, expert.pars, expert_type)
{
  beta <- expert.pars[[expert]]
  return(expert_pred(beta, X, "gaussian"))
}



marginal_effects <- function(obj)
{
  expt.nms <- c(obj[["expert.pars.nms"]])
  gate.nms <- obj[["gate.pars.nms"]]
  expert.type <- obj[["expert.type"]]
  expert.pars <- obj[["expert.pars"]]
  gate.pars <- obj[["gate.pars"]]
  experts <- obj[["expert.nms"]]
  gates <- obj[["gate.nodes"]]
  ln <- obj[["list_priors"]]
  tree <- obj[["tree"]]
  VCV <- obj[["full.vcv"]][["sandwich"]]
  
  expt.only.nms <- setdiff(expt.nms, c(gate.nms, "Dispersion"))
  gate.only.nms <- setdiff(gate.nms, expt.nms)

  # Calucate the building blocks
  #  gpp: Full mixture weights for each observation
  #  delta_m: marginal effect of the gating network
  #  expert_hats: fitted values for each expert
  gpp <- napply(experts, function(g) gate_path_product("0", g, ln))
  delta_m <- napply(experts, tree_expert_margin, gate.pars, ln)
  delta_m <- lapply(delta_m, `colnames<-`, gate.nms)
  
  expert_hats <- napply(experts, fitted_expert, obj[["X"]], expert.pars, expert.type)

  # Gating network allocations
  tree_expert_margins <- do.call(rbind, lapply(delta_m, colMeans))
  colnames(tree_expert_margins) <- gate.nms
  
  # Marginal Effects w.r.t gating network. Cumulative product
  row_multiply <- function(M, m)
  {
    dimnames(M) <- list(NULL, gate.nms)
    out <- sweep(M, 1, m, FUN=`*`)
    out
  }
  gate_margins <- mapply(row_multiply, delta_m, expert_hats, SIMPLIFY=FALSE)

  # Marginal Effects w.r.t the experts
  weighted_parameters <- function(W, b, exp_type)
  {
    if (exp_type == "gaussian")
      b <- b[-length(b)]
    out <- c(W) %o% b
    dimnames(out) <- list(NULL, setdiff(expt.nms, "Dispersion"))
    return(out)
  }
  expert_margins <- mapply(weighted_parameters, gpp, expert.pars, expert.type,
                           SIMPLIFY=FALSE)
  
  # Full marginal effects for variables that appear in the network
  # and the experts
  combine_margins <- function(gate, expert)
  {
    nms <- intersect(gate.nms, expt.nms)
    return(gate[, nms, drop=FALSE] + expert[, nms, drop=FALSE])
  }
  full_margins <- mapply(combine_margins,
                         gate_margins,
                         expert_margins,
                         SIMPLIFY=FALSE)
  full_margins   <- Reduce(`+`, full_margins)
  expert_margins <- Reduce(`+`, expert_margins)
  gate_margins   <- Reduce(`+`, gate_margins)
  
  # start of inference for average marginal effect
  # vectorize it to all nodes

  ASYVAR <- vector("list", length(experts))
  
  for (ii in seq_along(experts)) {
    
    m <- experts[ii]
    
    asy.var <- 0
    
    for (k in experts) {
      
      Delta_1 <- Delta_m_partial_theta_k(m, k, gate.nms, expt.nms, gates,
                                         experts, ln, obj[['Z']], obj[['X']],
                                         gate.pars, expert.pars, delta_m,
                                         expert_hats)
      
      for (j in experts) {
        
        Delta_2 <- Delta_m_partial_theta_k(m, j, gate.nms, expt.nms, gates,
                                           experts, ln, obj[['Z']], obj[['X']],
                                           gate.pars, expert.pars, delta_m,
                                           expert_hats)
        
        asy.var <- asy.var + asymptotic_variance(Delta_1, VCV, Delta_2, gates, experts,
                                                 gate.nms, expt.nms)
      }
    }
    
    ASYVAR[[ii]] <- asy.var
  }
  
  # Notes for Docs:
  #    *_margins are non-overlapping sets
  #    gate_expert_margins - sums to zero vector
  #                        - measures the "pull" of expert i on each variable j;
  #                          ignore intercept term; (or better yet what does it mean?)
  return(list(gate_margins        = colMeans(gate_margins[, gate.only.nms, drop=FALSE]),
              expert_margins      = colMeans(expert_margins[, expt.only.nms, drop=FALSE]),
              full_margins        = colMeans(full_margins),
              gate_expert_margins = t(tree_expert_margins),
              asyvar              = ASYVAR
              ))
}


asymptotic_variance <- function(D1, VCV, D2, gates, experts, gate.nms, expt.nms)
{
  nvrb <- ncol(D1[["0"]])
  
  vrbnms <- colnames(D1[["0"]])
  
  idx_to_gradient <- function(node, D, i)
  {
    if (node %in% gates) {
      nms <- gate.nms
    } else {
      nms <- setdiff(expt.nms, "Dispersion")
    }
    nr <- length(nms)
    G <- matrix(0, nr, nvrb)
    if (length(D[[node]]) > 1) {
      g <- D[[node]][i, ]
      for (v in seq_along(nms)) {
        i <- which(vrbnms==nms[v])
        G[v, i] <- g[i]
      }
    }
    if (node %in% experts) {
      G <- rbind(G, 0)
    }
    return(G)
  }

  out <- 0
  for (r in seq_len(nrow(D1[["0"]]))) {
    d1 <- do.call(rbind, lapply(names(D1), idx_to_gradient, D=D1, i=r))
    d2 <- do.call(rbind, lapply(names(D1), idx_to_gradient, D=D2, i=r))
    out <- out + t(d1) %*% VCV %*% d2
  }
  return(out)
}




# Also known as brackets for a node
# Eq 50 (G = )
var_wtd_gate_pars <- function(node, gate.pars, ln, Z)
{
  gs <- ln[[node]]
  
  # Generalized to J splits
  pars <- gate.pars[[node]]
  if (!is(pars, "list"))
    pars <- list(pars)
  # Set parameters branch J to 0
  pars[[length(pars) + 1]] <- rep(0, length(pars[[1]]))
  
  
  nc <- ncol(gs)
  nr <- nrow(gs)
  np <- length(pars[[1]])
  diffpars <- matrix(0, nrow=nr, ncol=np)
  
  out <- vector("list", length(pars))
  
  for (i in seq_along(pars)) {
    
    for (j in seq_along(pars)) {
      
      if (i == j) {
        
        g <- gs[, i] * (1 - gs[, i])
        samepars <- matrix(unlist(lapply(pars[[j]], function(x) x * g)), ncol=np)
        
      } else {
        
        g <- gs[, i] * gs[, j]
        diffpars <- diffpars + matrix(unlist(lapply(pars[[j]], function(x) x * g)), ncol=np)
        
      }
    }
    out[[i]] <- samepars - diffpars
  }
  
  return(out)
}



Delta_m_partial_theta_k <- function(m, k, gate_vrb, expt_vrb, gates, experts, ln, Z, X,
                                    gate.pars, expert.pars, delta_m, expert_hats)
{
  
  m_nodes <- unlist(ancestors(m))
  k_nodes <- unlist(ancestors(k))
  
  shared_nodes <- intersect(m_nodes, k_nodes)
  shared_gate_nodes <- intersect(shared_nodes, gates)
  
  # Q: How do we handle the intercept term between the two structures
  #    Currently handing it as a shared variable, which is probs wrong
  gate_grad <- napply(gates, function(x) 0)
  
  for (nn in shared_gate_nodes) {
    gate_grad[[nn]] <- find_marginal_gate_gradient(m, gate.pars, expert.pars,
                                                   expt.nms, ln, Z, expert_hats,
                                                   gate_vrb, expt_vrb)
  }
  
  expt_grad <- napply(experts, function(x) 0)
  
  if (m == k) {
    expt_grad[[m]] <- find_marginal_expert_gradient(m, gate.pars, expert.pars,
                                                    delta_m, ln, X, expert_hats,
                                                    gate_vrb, expt_vrb)
  }
  
  out <- c(gate_grad, expt_grad)
}

find_marginal_gate_gradient <- function(expert, gate.pars, expert.pars, expt.nms,
                                        ln, Z, expert_hats, gate_vrb, expt_vrb)
{
  
  expt_vrb <- setdiff(expt_vrb, "Dispersion")
  shared_vrb <- intersect(gate_vrb, expt_vrb)
  all_vrb <- union(gate_vrb, expt_vrb)
  
  # Start out with a matrix of zeros. Num of columns equals total number of parameters
  out <- matrix(0, nrow=nrow(Z), ncol=length(all_vrb), dimnames=list(NULL, all_vrb))
  
  # All variables in gating network have the following gradient
  grd <- delta_m_partial_omega(expert, gate.pars, ln, Z)
  grd <- lapply(grd, function(x) sweep(x, 1, expert_hats[[expert]], `*`))[[1]]

  # Variables that appear in both structures need an additive gradient term
  beta <- expert.pars[[expert]]
  beta <- beta[-length(beta)] # Remove dispersion parameter
  names(beta) <- expt_vrb 
  gpp_partial <- gpp_m_partial_omega(expert, gate.pars, ln, Z)[[1]]
  gpp_beta <- sweep(gpp_partial[, shared_vrb], 2, beta[shared_vrb], `*`)
  
  # Sum by column names
  out[, gate_vrb] <- grd[, gate_vrb]
  out[, shared_vrb]  <- out[, shared_vrb] + gpp_beta[, shared_vrb]
  
  return(out)
}


find_marginal_expert_gradient <- function(expert, gate.pars, expert.pars, delta_m,
                                          ln, X, expert_hats, gate_vrb, expt_vrb)
{
  expt_vrb <- setdiff(expt_vrb, "Dispersion")
  shared_vrb <- intersect(gate_vrb, expt_vrb)
  all_vrb <- union(gate_vrb, expt_vrb)
  
  # Start out with a matrix of zeros. Num of columns equals total number of HME parameters
  out <- matrix(0, nrow=nrow(X), ncol=length(all_vrb), dimnames=list(NULL, all_vrb))
  
  # All variables in gating network have the following gradient
  gpp_m <- gate_path_product("0", expert, ln)
  ones <- matrix(1, nrow=nrow(X), ncol=ncol(X), dimnames=list(NULL, expt_vrb))
  grd <- sweep(ones, 1, gpp_m, `*`)
  
  # Variables that appear in both structures need an additive gradient term
  delta_m_X <- delta_m[[expert]][, shared_vrb] * X[, shared_vrb]
  
  # Sum by column names
  out[, expt_vrb] <- grd[, expt_vrb]
  out[, shared_vrb]  <- out[, shared_vrb] + delta_m_X[, shared_vrb]
  
  return(out)
}


obj <- hme_2w

debugonce(Delta_m_partial_theta_k)

aa <- Delta_m_partial_theta_k("0.1", "0.2", obj[["gate.pars.nms"]], obj[["expert.pars.nms"]],
                              obj[["gate.nodes"]], obj[["expert.nms"]], obj[["list_priors"]],
                              obj[["Z"]], obj[["X"]], obj[["gate.pars"]])

debugonce(marginal_effects)
aa <- marginal_effects(hme_2w)



delta_m_partial_omega <- function(expert, gate.pars, ln, Z)
{
  
  # find all nodes from root to expert
  nodes <- unlist(ancestors(expert))
  paths <- as.integer(unlist(last_split(nodes)))
  nodes <- nodes[-length(nodes)]
  W <- tree_expert_margin(expert, gate.pars, ln, weighted=FALSE)
  
  gpp <- gate_path_product("0", expert, ln)
  
  delta_m_partial <- function(node, node_split)
  {

    G <- var_wtd_gate_pars(node, gate.pars=gate.pars, ln=ln, Z=Z)
    npaths <- length(gate.pars[node])
    out <- vector("list", length(npaths))
    
    for (p in seq_len(npaths)) {
      
      path_bool <- p == node_split
      
      sgn <- 2 * as.integer(path_bool) - 1
      
      g <- ln[[node]][, node_split]
      
      const <- as.integer(path_bool) - g
      
      multi_z <- (sgn * (1 - g) * W[[node_split]] - G[[node_split]]) * Z
      
      M <- sweep(multi_z, 1, const, `+`)
      
      out[[p]] <- sweep(M, 1, gpp, `*`)
    }
    return(out)
  }
  
  tst <- mapply(delta_m_partial, nodes, paths)
  
  return(tst)
}


gpp_m_partial_omega <- function(expert, gate.pars, ln, Z)
{
  
  # find all nodes from root to expert
  nodes <- unlist(ancestors(expert))
  paths <- as.integer(unlist(last_split(nodes)))
  nodes <- nodes[-length(nodes)]
  
  gpp <- gate_path_product("0", expert, ln)
  
  gpp_partial <- function(node, node_split)
  {
    npaths <- length(gate.pars[node])
    out <- vector("list", length(npaths))
    
    for (p in seq_len(npaths)) {
      
      path_bool <- p == node_split
      
      g <- ln[[node]][, node_split]
      
      const <- gpp * (as.integer(path_bool) - g)
      
      out[[p]] <- sweep(Z, 1, const, `*`)
    }
    return(out)
  }

  tst <- mapply(gpp_partial, nodes, paths)
  
  return(tst)
}




