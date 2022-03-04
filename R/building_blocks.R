#---------------------------------------------------------------------
# Set of functions to transverse the HME tree

# given a node address, provide names of that nodes children
# if `d` is terminal node, return NULL
children <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (is(d, "NULL"))
      return(NULL)
    o <- nodes[grep(paste0("^", x, ".[0-9]$"), nodes)]
    if (length(o) == 0) {
      return(NULL)
    } else {
      return(o)
    }
  }
  lapply(d, f_)
}

# given a node address, provide name of that node's parent
# if `d` is the root node, return NULL
parent <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (x == "0") {
      return(NULL)
    } else {
      return(regmatches(x, regexpr(".*(?=\\.)", x, perl=TRUE)))
    }
  } 
  lapply(d, f_)
}


# given a node address, provide ancestry of that node back to the root
ancestors <- function(d)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  fetch_history <- function(y)
  {
    f_ <- function(x, z) 
    {
      paste(c(x, z), collapse=".") 
    }
    if (y == "0") {
      return("0")
    } else {
      return(Reduce(f_, unlist(strsplit(y, "\\.")), accumulate=T))
    }
  }
  lapply(d, fetch_history)
}


# given a node address, list all children, grand-children, ect of that node
progeny <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(z)
  {
    node_well <- c()
    look_forward <- function(y)
    {
      o <- unlist(children(y, nodes))
      if (!is(o, "NULL")) {
        node_well <<- c(o, node_well)
        look_forward(o)
      }
    }
    look_forward(z)
    if (is(node_well, "NULL")) {
      return(NULL)
    } else {
      return(sort(node_well))
    }
  }
  lapply(d, f_)
}


# given a node address, find all nodes that share the same parent
siblings <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    p <- unlist(parent(x, nodes))
    return(unlist(children(p, nodes)))
  }
  lapply(d, f_)
}


# given a node address, find all nodes that share the same depth
generation <- function(x, nodes)
{
  nodes[which(nchar(nodes) == nchar(x))]
}


# is `d` a terminal node of the tree (is it childless?)
is_terminal <- function(d, nodes)
{
  if(!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    bool <- FALSE
    if (is(unlist(children(x, nodes)), "NULL"))
      bool <- TRUE
    return(bool)
  }
  lapply(d, f_)
}



last_split <- function(d)
{
  if(!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (x == "0" || is(x, "NULL")) {
      return(NULL)
    } else {
      o <- unlist(strsplit(x, "\\."))
      return(o[length(o)])
    }
  }
  lapply(d, f_)
}


all_but_last_split <- function(d)
{
  if(!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (x == "0" || is(x, "NULL")) {
      return(NULL)
    } else {
      o <- unlist(strsplit(x, "\\."))
      return(paste0(o[-length(o)], collapse="."))
    }
  }
  lapply(d, f_)
}


# how deep is the tree
tree_depth <-function(nodes)
{
  if(!is(nodes, "list"))
    nodes <- as.list(nodes)
  
  f_ <- function(x) length(unlist(gregexpr("\\.", x)))
  max(unlist(lapply(nodes, f_)))
}


expert_index <- function(hme_type, nodes)
{
  idx <- unlist(is_terminal(nodes, nodes))
  if (hme_type == "hme") {
    idx[idx] <- seq_len(sum(idx))
  }
  if (hme_type == "hmre") {
    lsp <- as.integer(unlist(last_split(nodes)))
    idx <- c(0, as.integer(idx[-1]) * lsp)
  }
  names(idx) <- nodes
  return(idx)
}


napply <- function(x, f_, ...)
{
  out <- lapply(x, f_, ...)
  names(out) <- x
  return(out)
}




# ---------------------------------------------------------------------------
# Set of functions to obtain prior and posterior weights

par_to_exp_gate_paths <- function(node, gate_par_list, X)
{
  pars <- gate_par_list[[node]]
  if (!is.list(pars)) {
    pars <- list(pars)
  }
  f_ <- function(p)
  {
    return(exp(X %*% p))
  }
  numerator <- lapply(pars, f_)
  G <- matrix(unlist(numerator), ncol=length(pars))
  G <- cbind(G, 1)
  return(G)
}

par_to_gate_paths <- function(node, gate_par_list, X)
{
  G <- par_to_exp_gate_paths(node, gate_par_list, X)
  return(sweep(G, 1, rowSums(G), `/`))
}


par_to_expert_dens <- function(node, expert_type, expert_par_list, Y, X, ...)
{
  if (expert_type == "gaussian") {
    parm <- expert_par_list[[node]]
    beta <- parm[1:(length(parm) - 1)]
    variance <- exp(parm[length(parm)])
    mu <- X %*% beta
    return(dnorm(Y, mean=mu, sd=sqrt(variance), ...))
  } else if (expert_type == "bernoulli") {
    # if Y is a bernoulli variable {1, 0}
    g <- par_to_gate_paths(node, expert_par_list, X)
    # pmax(g[,1] * Y, g[,2] * (1 - Y))
    g[,1] * Y + g[,2] * (1 - Y)
  }
}


"Input: 
       branches - list of node names (i.e. '0.1.1')
       ln       - list of gating network split probabilities
 Output: 
       list of path values between the given node ('0.1.1') and its
       parent ('0.1')"
gate_path_values <- function(branches, ln)
{
  if (!is(branches, "list"))
    branches <- as.list(branches)
  
  f_ <- function(b)
  {
    if (b == "0")
      return(NULL)
    branch <- unlist(last_split(b))
    node <- unlist(all_but_last_split(b))
    ln[[node]][, as.integer(branch), drop=FALSE]
  }
  lapply(branches, f_)
}


inter_node_paths <- function(geezer, youngin, ln)
{
  "product path from node1 down to node2"
  a1 <- unlist(ancestors(geezer))
  a2 <- unlist(ancestors(youngin))
  if (!geezer %in% a2) {
    stext <- "Node1 [%s] must be an ancestor of node2 [%s]."
    stop(sprintf(stext, node1, node2))
  }
  return(gate_path_values(setdiff(a2, a1), ln))
}


"Input:
       geezer  - the name of a gating node
       youngin - the name of a gating node or expert that is an
                 ancestor of `geezer`
       ln      - list of gating network split probabilities
 Output:
       cumulative product of the paths from the `geezer` node to the
       `youngin` node"
gate_path_product <- function(geezer, youngin, ln)
{
  inp <- inter_node_paths(geezer, youngin, ln)
  Reduce(`*`, inp)
}



"Input:
       x - the integer names of an expert
 Output:
       list of the current density estimates of the experts"
get_expert_densities <- function(x, le)
{
  f_ <- function(y)
  {
    le[[y]]
  }
  lapply(x, f_)
}


"Calculate the prior weight for a gating node

 Input:
       node - any given non-root node
 Output:
       cumulative product of the paths from the root node to the
       input node (ie the prior weights for the subtree starting
       from`node`)"
prior_weights <- function(node, ln)
{
  "prior weights go from root node down"
  gate_path_product("0", node, ln)
}


"Calculate the posterior weights for a gating node
 
 Input:
       node    - any gating node
       treestr - list of all gating and expert nodes in the HME
       ln      - list of gating network split probabilities
       le      - list of expert densities
 Output:
       the posterior weights of the subtree starting from `node`"
posterior_weights <- function(node, treestr, ln, le)
{
  "posterior weights go from experts up"
  p <- unlist(progeny(node, treestr))
  terminals <- p[unlist(is_terminal(p, treestr))]
  childs <- unlist(children(node, treestr))
  
  Gs <- napply(terminals, function(x) gate_path_product(node, x, ln))
  # first match terminal paths with their densities
  f_ <- function(x)
  {
    idx <- which(unlist(terminals) == x)
    Reduce(`+`, Gs[idx]) * le[[x]]
  }
  H <- napply(names(le), f_)
  # sum up and standardize for the appropriate number of branches
  g_ <- function(x)
  {
    Reduce(`+`, H[grepl(x, names(H))])
  }
  HH <- matrix(unlist(lapply(childs, g_)), ncol=length(childs))
  sweep(HH, 1, rowSums(HH), `/`)
}


"Calculate the joint posterior weight for a gating node (weight for the glm)

 Output:
       the cumulative product of posterior weights from the root to `node`"
joint_posterior_weight <- function(node, lp, rp)
{
  # not feasible for the root node
  lrp <- length(rp)
  nlp <- nrow(lp[[1]])
  stopifnot(lrp == 1 || lrp == nlp)
  if (node == "0") {
    if (lrp == 1)
      return(rep(1, nlp))
    else
      return(rp)
  }
  return(rp * gate_path_product("0", node, ln=lp))
}



Qlog_likelihood <- function(treestr, ln, lp, ld, rp)
{
  expert.nodes <- treestr[unlist(is_terminal(treestr, treestr))]
  full_log_path <- function(node)
  {
    inp <- inter_node_paths("0", node, ln)
    logs <- lapply(c(inp, list(ld[[node]])), log)
    sumlogs <- Reduce(`+`, logs)
    post <- joint_posterior_weight(node, treestr, lp, rp)
    return(post * sumlogs)
  }
  lst <- lapply(expert.nodes, full_log_path)
  sum(unlist(lst))
}

log_likelihood <- function(treestr, ln, ld)
{
  expert.nodes <- treestr[unlist(is_terminal(treestr, treestr))]
  S <- simplify2array(expert_lik_contr(expert.nodes, ld, ln))
  sum(S)
}


init_gate_node_pars <- function(node, tree, n)
{
  nchilds <- length(unlist(children(node, tree)))
  limit <- 0.1
  if (nchilds == 2) {
    return(runif(n, -limit, limit))
  } else {
    return(lapply(seq_len(nchilds - 1), function(x) runif(n, -limit, limit)))
  }
}

