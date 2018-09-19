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
  # get list of expert names
  f_ <-function(expert)
  {
    colMeans(tree_expert_margin(expert, gate.pars, ln))
  }
  tst <- napply(experts, f_) # function(x) colMeans(tree_expert_margin(x, gate.pars, ln)))
  return(tst)
}


"the marginal effects for each expert should sum to zero"
# colMeans(tree_expert_margin("0.1", tst2$gate.pars, tst2$list_priors))
# colMeans(tree_expert_margin("0.2", tst2$gate.pars, tst2$list_priors))
# colMeans(tree_expert_margin("0.3", tst2$gate.pars, tst2$list_priors))

# colMeans(tree_expert_margin("0.1.1", tst$gate.pars, tst$list_priors))
# colMeans(tree_expert_margin("0.1.2", tst$gate.pars, tst$list_priors))
# colMeans(tree_expert_margin("0.2", tst$gate.pars, tst$list_priors))
