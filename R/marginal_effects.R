gate_margin <- function(node, gate.pars, ln)
{
  pars <- gate.pars[[node]]
  if (!is(pars, "list"))
    pars <- list(pars)
  pars[[length(pars) + 1]] <- rep(0, length(pars[[1]]))
  
  nn <- length(pars)
  margs <- vector("list", nn)
  for (ii in seq_len(nn)) {
    margs[[ii]] <- sapply(pars[[ii]], `*`, ln[[node]][, ii])
  }
  par.means <- Reduce(`+`, margs)
  f_ <- function(x) -sweep(par.means, 2, pars[[x]], FUN=`-`)
  delta <- lapply(seq_len(nn), f_)
  return(delta)
}


# Marginal Gating Values routed to each expert: \partial pi / \parial Z
tree_expert_margin <- function(expert, gate.pars, ln)
{
  gpp <- gate_path_product("0", expert, ln)
  a <- unlist(ancestors(expert))[-1]
  abls <- unlist(lapply(a, all_but_last_split))
  lspl <- unlist(lapply(a, last_split))
  treasure <- vector("list", length(a))
  for (nn in seq_along(a)) {
    int <- as.integer(lspl[nn])
    treasure[[nn]] <- gate_margin(abls[[nn]], gate.pars, ln)[[int]]
  }
  X <- Reduce(`+`, treasure)
  out <- sweep(X, 1, gpp, `*`)
  return(out)
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
  ln <- obj[["list_priors"]]
  
  expt.only.nms <- setdiff(expt.nms, c(gate.nms, "Dispersion"))
  gate.only.nms <- setdiff(gate.nms, expt.nms)

  # Calucate the building blocks
  #  gpp: Full mixture weights for each observation
  #  tree_margins: marginal effect of the gating network
  #  expert_hats: fitted values for each expert
  gpp <- napply(experts, function(g) gate_path_product("0", g, ln))
  tree_margins <- napply(experts, tree_expert_margin, gate.pars, ln) #TODO: add gpp
  expert_hats <- napply(experts, fitted_expert, obj[["X"]], expert.pars, expert.type)

  # Gating network allocations
  tree_expert_margin <- do.call(rbind, lapply(tree_margins, colMeans))
  colnames(tree_expert_margin) <- gate.nms
  
  # Marginal Effects w.r.t gating network
  row_multiply <- function(M, m)
  {
    dimnames(M) <- list(NULL, gate.nms)
    out <- sweep(M, 1, m, FUN=`*`)
    out
  }
  gate_margins <- mapply(row_multiply, tree_margins, expert_hats, SIMPLIFY=FALSE)

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
  
  # Full marginal effects for variables that appear that appear in the network
  # and the experts
  combine_margins <- function(gate, expert)
  {
    nms <- setdiff(intersect(gate.nms, expt.nms), "(Intercept)")
    return(gate[, nms, drop=FALSE] + expert[, nms, drop=FALSE])
  }
  full_margins <- mapply(combine_margins, gate_margins, expert_margins,
                         SIMPLIFY=FALSE)
  full_margins <- Reduce(`+`, full_margins)
  expert_margins <- Reduce(`+`, expert_margins)
  gate_margins <- Reduce(`+`, gate_margins)
  
  return(list(gate_margins=colMeans(gate_margins[, gate.only.nms]),
              expert_margins=colMeans(expert_margins[, expt.only.nms]),
              full_margins=colMeans(full_margins),
              gate_expert_margins=tree_expert_margin))
}






