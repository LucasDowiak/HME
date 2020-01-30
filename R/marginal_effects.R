# for individual gating nodes, Subtract par means from path i
#    Parenthesis in equation 45
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
# equation 46
tree_expert_margin <- function(expert, gate.pars, ln)
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
  X <- Reduce(`+`, treasure)
  
  # normalized by the gate path product
  delta_m <- sweep(X, 1, gpp, `*`)
  return(delta_m)
}


margin_matrix <- function(experts, gate.pars, ln)
{
  # Avergage Marginal Effect
  f_ <-function(expert)
  {
    colMeans(tree_expert_margin(expert, gate.pars, ln))
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
  expt.nms <- obj[["expert.pars.nms"]]
  gate.nms <- obj[["gate.pars.nms"]]
  expert.type <- obj[["expert.type"]]
  expert.pars <- obj[["expert.pars"]]
  gate.pars <- obj[["gate.pars"]]
  experts <- obj[["expert.nms"]]
  gates <- obj[["gate.nodes"]]
  ln <- obj[["list_priors"]]
  tree <- obj[["tree"]]
  
  expt.only.nms <- setdiff(expt.nms, c(gate.nms, "Dispersion"))
  gate.only.nms <- setdiff(gate.nms, expt.nms)

  # Calucate the building blocks
  #  gpp: Full mixture weights for each observation
  #  tree_margins: marginal effect of the gating network
  #  expert_hats: fitted values for each expert
  gpp <- napply(experts, function(g) gate_path_product("0", g, ln))
  delta_m <- napply(experts, tree_expert_margin, gate.pars, ln) #TODO: add gpp
  expert_hats <- napply(experts, fitted_expert, obj[["X"]], expert.pars, expert.type)

  # Gating network allocations
  tree_expert_margins <- do.call(rbind, lapply(delta_m, colMeans))
  colnames(tree_expert_margins) <- gate.nms
  
  # Marginal Effects w.r.t gating network
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
  full_margins <-   Reduce(`+`, full_margins)
  expert_margins <- Reduce(`+`, expert_margins)
  gate_margins <-   Reduce(`+`, gate_margins)
  
  # vectorize it to all nodes
  dpo <-   napply(gates, delta_m_partial_omega,
                  obj[['tree']],
                  delta_m,
                  gate.pars,
                  ln,
                  obj[["Z"]])
  
  gpo <-   napply(gates, gpp_partial_omega,
                  obj[['tree']],
                  gate.pars,
                  ln,
                  obj[["Z"]])
  
  # implememt equation 49

  # Delta Method w.r.t gating network
  delta_method_for_omega <- function(node, tree)
  {
    
  }
  
  # Notes for Docs:
  #    *_margins are non-overlapping sets
  #    gate_expert_margins - sums to zero vector
  #                        - measures the "pull" of expert i on each variable j;
  #                          ignore intercept term; (or better yet what does it mean?)
  return(list(gate_margins        = colMeans(gate_margins[, gate.only.nms]),
              expert_margins      = colMeans(expert_margins[, expt.only.nms]),
              full_margins        = colMeans(full_margins),
              gate_expert_margins = t(tree_expert_margins),
              tst = tst))
}






# Also known as brackets for a node
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
  
  return(out[-length(out)])
}




delta_m_partial_omega <- function(node, tree, delta_m, gate.pars, ln, Z)
{
  
  offspring <- unlist(progeny(node, tree))
  experts <- offspring[unlist(is_terminal(offspring, tree))]
  
  # 
  delta_m_partial <- function(expert)
  {
    gpp <- gate_path_product("0", expert, ln)
    g <- ln[[node]][, 1]
    const <- c(gpp * (1 - g))
    bracket <- var_wtd_gate_pars(node, gate.pars, ln, Z)
    multi_z <- delta_m[[expert]] * (1 - g) + sweep(bracket[[1]], 1, gpp, `*`)
    return(const + (multi_z * Z))
  }
  
  out <- napply(experts, delta_m_partial)
  
  return(out)
}


gpp_partial_omega <- function(node, tree, gate.pars, ln, Z)
{
  offspring <- unlist(progeny(node, tree))
  experts <- offspring[unlist(is_terminal(offspring, tree))]
  
  # 
  gpp_partial <- function(expert)
  {
    gpp <- gate_path_product("0", expert, ln)
    g <- ln[[node]][, 1]
    const <- gpp * (1 - g)
    return(sweep(Z, 1, const, `*`))
  }
  
  out <- napply(experts, gpp_partial)
  
  return(out)
}

debugonce(marginal_effects)
aa <- marginal_effects(new_hme_1)
