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
  N <- obj[["N"]]
  
  # Venn-Diagram of expert and gate variable names
  expt.nms <- setdiff(expt.nms, "Dispersion")
  expt.only.nms <- setdiff(expt.nms, gate.nms)
  gate.only.nms <- setdiff(gate.nms, expt.nms)
  shared.nms <- intersect(gate.nms, expt.nms)
  all.nms <- union(gate.nms, expt.nms)
  
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
    dimnames(out) <- list(NULL, expt.nms)
    return(out)
  }
  expert_margins <- mapply(weighted_parameters, gpp, expert.pars, expert.type,
                           SIMPLIFY=FALSE)
  
  # Full marginal effects for variables that appear in the network
  # and the experts
  combine_margins <- function(gate, expert)
  {
    return(gate[, shared.nms, drop=FALSE] + expert[, shared.nms, drop=FALSE])
  }
  both_margins <- mapply(combine_margins,
                         gate_margins,
                         expert_margins,
                         SIMPLIFY=FALSE)
  
  both_margins   <- Reduce(`+`, both_margins)
  expert_margins <- Reduce(`+`, expert_margins)
  gate_margins   <- Reduce(`+`, gate_margins)
  
  margins <- matrix(0, N, length(all.nms), dimnames = list(NULL, all.nms))
  
  if (length(gate.only.nms) > 0)
    margins[, gate.only.nms] <- gate_margins[, gate.only.nms]
  
  if (length(expt.only.nms) > 0)
    margins[, expt.only.nms] <- expert_margins[, expt.only.nms]
  
  if (length(shared.nms) > 0)
    margins[, shared.nms] <- both_margins[, shared.nms]
  
  
  # inference for average marginal effect of the entire model
  shared.asy.var <- gateonly.asy.var <- exptonly.asy.var <- 0
  #asy.var <- as.list(rep(0, 3))
  for (k in experts) {
      
    Delta <- Delta_m_partial_theta_k(k, gate.nms, expt.nms,
                                     gates, experts, ln, obj[['Z']], obj[['X']],
                                     gate.pars, expert.pars, delta_m, expert_hats)
    
    k_asy_var <- asymptotic_variance(Delta, VCV, gates, experts, gate.nms, expt.nms)
    
    shared.asy.var <- shared.asy.var + k_asy_var[['shared']]
    gateonly.asy.var <- gateonly.asy.var + k_asy_var[['gateonly']]
    exptonly.asy.var <- exptonly.asy.var + k_asy_var[['exptonly']]
  }
  
  cat("\n--------------------------\nGate Margins:\n")
  gm <- colMeans(gate_margins)
  se <- sqrt(diag(gateonly.asy.var))[names(gm)]
  zstat <- gm / se
  tbl <- cbind(gm, se, zstat, 2 * pnorm(-abs(zstat)))
  colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
  printCoefmat(tbl, digits = 3)
  
  cat("\n\n\n--------------------------\nExpert Margins:\n")
  em <- colMeans(expert_margins)
  se <- sqrt(diag(exptonly.asy.var))[names(em)]
  zstat <- em / se
  tbl <- cbind(em, se, zstat, 2 * pnorm(-abs(zstat)))
  colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
  printCoefmat(tbl, digits = 3)
  
  cat("\n\n\n--------------------------\nFull Margins:\n")
  fm <- colMeans(margins)
  se <- sqrt(diag(shared.asy.var))[names(fm)]
  zstat <- fm / se
  tbl <- cbind(fm, se, zstat, 2 * pnorm(-abs(zstat)))
  colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
  printCoefmat(tbl, digits = 3)
  
  cat("\n\n\n--------------------------\nGating Network:\n")
  print(round(t(tree_expert_margins), 4))
  
  # Notes for Docs:
  #    *_margins are non-overlapping sets
  #    tree_expert_margins - sums to zero vector
  #                        - measures the "pull" of expert i on each variable j;
  #                          ignore intercept term; (or better yet what does it mean?)
  return(
    invisible(
      list(gate_margins          = gm,
           expert_margins        = em,
           margins               = fm,
           gate_expert_margins   = t(tree_expert_margins),
           shared_asy_var        = shared.asy.var,
           gateonly_asy_var      = gateonly.asy.var,
           exptonly_asy_var      = exptonly.asy.var
      )
    )
  )
}


asymptotic_variance <- function(D1, VCV, gates, experts, gate.nms, expt.nms)
{
  nodes <- names(D1)
  allnms <- union(gate.nms, expt.nms)
  nvrb <- length(allnms)
  tree <- c(gates, experts)
  
  col_means <- function(x)
  {
    if (is(x, "matrix")) {
      return(colMeans(x))
    } else {
      lapply(x, col_means)
    }
  }
  
  find_by_name <- function(x, name)
  {
    if (name %in% names(x)) {
      return(x[[name]])
    } else {
      lapply(x, find_by_name, name)
    }
  }
  
  replace_with_zeros <- function(x)
  {
    replace_ <- function(z) {
      y <- rep(0, length(z))
      names(y) <- names(z)
      return(y)
    }
    
    if (is(x, "numeric")) {
      return(replace_(x))
    } else {
      lapply(x, replace_with_zeros)
    }
  }
  
  idx_to_gradient <- function(node, M)
  {
    if (node %in% gates) {
      nms <- gate.nms
    } else {
      nms <- expt.nms
    }
    nr <- length(nms)
    G <- matrix(0, nr, nvrb, dimnames = list(NULL, allnms))
    if (node %in% experts) {
      # dispersion parameter doesn't effect the marginal effects for a normal distribution
      # add a zero row
      G <- rbind(G, 0)
    }
    
    if (node %in% names(M)) {
      D <- M[[node]]
      if (!is(D, "list"))
        D <- list(D)
      nsplits <- length(D)
      
      out <- vector("list", nsplits)
      
      for (s in seq_len(nsplits)) {
        g <- D[[s]]
        for (v in seq_along(nms)) {
          i <- which(allnms==nms[v])
          G[v, i] <- g[nms[v]]
        }
        out[[s]] <- G
      }
    } else {
      out <- list(G)
    }
    return(out)
  }
  
  D1 <- lapply(D1, col_means)
  
  # apply it
  Dshared <- find_by_name(D1, 'shared')
  
  nodenms <- intersect(nodes, experts)
  Dgateonly <- find_by_name(D1, 'only')
  Dgateonly[nodenms] <- lapply(Dgateonly[nodenms], replace_with_zeros)
  
  nodenms <- intersect(nodes, gates)
  Dexptonly <- find_by_name(D1, 'only')
  Dexptonly[nodenms] <- lapply(Dexptonly[nodenms], replace_with_zeros)
  
  
  gradients <- list(napply(tree, idx_to_gradient, Dshared),
                    napply(tree, idx_to_gradient, Dgateonly),
                    napply(tree, idx_to_gradient, Dexptonly))
  names(gradients) <- c("d_shared", "d_gatesonly", "d_exptonly")
  
  outnms <- c("shared", "gateonly", "exptonly")
  out <- vector("list", length(outnms))
  names(out) <- outnms
  
  for (nn in seq_along(out)) {
    GRD <- do.call(rbind, unlist(gradients[[nn]], recursive=FALSE))
    out[[nn]] <- t(GRD) %*% VCV %*% GRD
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



Delta_m_partial_theta_k <- function(m, gate_vrb, expt_vrb, gates, experts, ln, Z, X,
                                    gate.pars, expert.pars, delta_m, expert_hats)
{
  
  m_anc <- unlist(ancestors(m))
  tree <- c(gates, experts)
  
  # Q: How do we handle the intercept term between the two structures
  #    Currently handing it as a shared variable, which is probs wrong
  gate_grad <- find_marginal_gate_gradient(m, gate.pars, expert.pars,
                                           expt.nms, ln, Z, expert_hats,
                                           gate_vrb, expt_vrb, tree)
  if (length(gate_grad) != length(m_anc) - 1)
    gate_grad <- list(gate_grad)
  
  names(gate_grad) <- intersect(m_anc, gates)
  
  
  expt_grad <- find_marginal_expert_gradient(m, gate.pars, expert.pars,
                                             delta_m, ln, X, expert_hats,
                                             gate_vrb, expt_vrb)
  expt_grad <- list(expt_grad)
  names(expt_grad) <- m
  
  out <- c(gate_grad, expt_grad)
}



find_marginal_gate_gradient <- function(expert, gate.pars, expert.pars, expt.nms,
                                        ln, Z, expert_hats, gate_vrb, expt_vrb,
                                        tree)
{
  
  shared_vrb <- intersect(gate_vrb, expt_vrb)
  all_vrb <- union(gate_vrb, expt_vrb)
  
  # All variables in gating network have the following gradient
  grd <- delta_m_partial_omega(expert, gate.pars, ln, Z, tree)
  grd <- lapply(grd, function(x) sweep(x, 1, expert_hats[[expert]], `*`)) # multiply by f (fitted mean)

  # Marginal effects of the regression; aka Beta
  beta <- expert.pars[[expert]]
  beta <- beta[-length(beta)] # Remove dispersion parameter
  names(beta) <- expt_vrb
  
  # gate_path_partial Jocobian wrt to Omega
  gpp_partial <- gpp_m_partial_omega(expert, gate.pars, ln, Z, tree)
  
  # contribution to the Jocobian from the Beta (parameters in both parts of the HME)
  gpp_beta <- lapply(gpp_partial, function(x) sweep(x[, shared_vrb], 2, beta[shared_vrb], `*`)) # Multiply by Beta
  
  # length(grd) is equal to the number of nodes along the path
  out <- vector("list", length(grd))
  
  for (nn in seq_along(out)) {
    
    lst <- list()
    
    G <- matrix(0, nrow=nrow(Z), ncol=length(all_vrb), dimnames=list(NULL, all_vrb))
    
    G[, gate_vrb] <- grd[[nn]][, gate_vrb]
    lst[['only']] <- G
    
    G[, shared_vrb]  <- G[, shared_vrb] + gpp_beta[[nn]][, shared_vrb]
    lst[['shared']] <- G
    
    out[[nn]] <- lst
  }

  
  return(out)
}


find_marginal_expert_gradient <- function(expert, gate.pars, expert.pars, delta_m,
                                          ln, X, expert_hats, gate_vrb, expt_vrb)
{
  shared_vrb <- intersect(gate_vrb, expt_vrb)
  all_vrb <- union(gate_vrb, expt_vrb)
  
  # Start out with a matrix of zeros. Num of columns equals total number of HME parameters
  G <- matrix(0, nrow=nrow(X), ncol=length(expt_vrb), dimnames=list(NULL, expt_vrb))

  # All variables in gating network have the following gradient
  gpp_m <- gate_path_product("0", expert, ln)
  ones <- matrix(1, nrow=nrow(X), ncol=ncol(X), dimnames=list(NULL, expt_vrb))
  grd <- sweep(ones, 1, gpp_m, `*`)
  
  # Variables that appear in both structures need an additive gradient term
  delta_m_X <- delta_m[[expert]][, shared_vrb] * X[, shared_vrb]
  
  # Sum by column names
  out <- list()
  
  G[, expt_vrb] <- grd[, expt_vrb]
  out[['only']] <- G
  
  G[, shared_vrb]  <- G[, shared_vrb] + delta_m_X[, shared_vrb]
  out[['shared']] <- G
  
  return(out)
}



# For all gating nodes on the path from the root node to expert m,
# return the Jocobian of the marginal effects wrt to Omega
delta_m_partial_omega <- function(expert, gate.pars, ln, Z, tree)
{
  
  # find all gate nodes from root to the expert
  nodes <- unlist(ancestors(expert))
  splits <- as.integer(unlist(last_split(nodes)))
  gnodes <- nodes[-length(nodes)]
  W <- tree_expert_margin(expert, gate.pars, ln, weighted=FALSE)
  
  gpp <- gate_path_product("0", expert, ln)
  
  delta_m_partial <- function(node, node_split)
  {

    G <- var_wtd_gate_pars(node, gate.pars=gate.pars, ln=ln, Z=Z)
    # npaths <- length(gate.pars[node])
    # Record all but the last split since J = 0
    paths <- as.integer(unlist(last_split(unlist(children(node, tree)))))
    paths <- paths[-length(paths)]
    
    out <- vector("list", length(paths))
    
    for (p in paths) {
      
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
  
  tst <- mapply(delta_m_partial, gnodes, splits)
  #tst <- lapply(gnodes, delta_m_partial)
  
  return(tst)
}



# For each expert m,
# return the Jocobian of the marginal effects wrt to Omega
gpp_m_partial_omega <- function(expert, gate.pars, ln, Z, tree)
{
  
  # find all nodes from root to expert
  nodes <- unlist(ancestors(expert))
  splits <- as.integer(unlist(last_split(nodes)))
  gnodes <- nodes[-length(nodes)]
  
  gpp <- gate_path_product("0", expert, ln)
  
  gpp_partial <- function(node, node_split)
  {
    # Record all but the last split since J = 0
    paths <- as.integer(unlist(last_split(unlist(children(node, tree)))))
    paths <- paths[-length(paths)]
    
    #npaths <- length(gate.pars[node])
    out <- vector("list", length(paths))
    
    for (p in paths) {
      
      path_bool <- p == node_split
      
      g <- ln[[node]][, node_split]
      
      const <- gpp * (as.integer(path_bool) - g)
      
      out[[p]] <- sweep(Z, 1, const, `*`)
    }
    return(out)
  }

  tst <- mapply(gpp_partial, gnodes, splits)
  
  return(tst)
}




