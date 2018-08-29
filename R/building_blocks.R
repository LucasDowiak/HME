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
      return(sort(c(z, node_well)))
    }
  }
  lapply(d, f_)
}


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


generation <- function(x, nodes)
{
  d <- tree_depth(nodes)
  if (x > d) {
    mssg <- "The generation [%d] must not be greater than the depth [%d]."
    stop(sprintf(mssg, x, d))
  }
  
  nodes[which(nchar(nodes) == 2 * x + 1)]
}


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
  #Map(function(y, z) if (y) z, b, d)d[unlist(b)]
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


tree_depth <-function(nodes)
{
  if(!is(nodes, "list"))
    nodes <- as.list(nodes)
  
  f_ <- function(x) length(unlist(gregexpr("\\.", x)))
  max(unlist(lapply(nodes, f_)))
}

" generalize to HMRE "
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

#expert_index("hme", sort(test1))
#expert_index("hmre", sort(test1))
#expert_index("hme", sort(test))
#expert_index("hmre", sort(test))


napply <- function(x, f_, ...)
{
  out <- lapply(x, f_, ...)
  names(out) <- x
  return(out)
}


which_experts_split <- function(d, nodes)
{
  experts <- unlist(is_terminal(nodes, nodes))
  childs <- unlist(children(d, nodes))
  f_ <- function(x)
  {
    p <- nodes[grepl(sprintf("^%s\\.", x), nodes)]
    if (length(p) == 0)
      return(x)
    else
      return(intersect(p, nodes[experts]))
  }
  lapply(childs, f_)
}


# ---------------------------------------------------------------------------
# Set of functions to obtain prior and posterior weights

test <- c("0",
          "0.1", "0.2",
          "0.1.1", "0.1.2", "0.2.1", "0.2.2",
          "0.1.1.1", "0.1.1.2")

test1 <- c("0",
           "0.1", "0.2", "0.3",
           "0.1.1", "0.1.2", "0.1.3", "0.2.1", "0.2.2", "0.2.3",
           "0.1.3.1", "0.1.3.2", "0.1.3.3", "0.2.3.1", "0.2.3.2", "0.2.3.3")

# Take in matrix
# return list
# match_terminals_dens <- function(hme_type=c("hme", "hmre"), treestr)
# {
#   
#   hme_type <- match.arg(hme_type)
#   idx <- expert_index(hme_type, treestr)
#   if (hme_type == "hme") {
#     # run expert_index 
#     # 
#   }
#   if (hme_type == "hmre")
#   {
#     last_split()
#   }
# }
# 
# 
# random_split <- function(nr, nc)
# {
#   M <- matrix(nrow=nr, ncol=nc)
#   M[, 1] <- runif(nr)
#   i <- 2
#   while (i < nc) {
#     M[, i] <- runif(nr, max = 1 - rowSums(M, na.rm=T))
#     i <- i + 1
#   }
#   M[, nc] <- 1 - rowSums(M, na.rm=T)
#   return(M)
# }
# 
# NE <- 3
# NR <- 10
# 
# expnames <- test1[unlist(is_terminal(test1, test1))]
# gatenames <- setdiff(test1, expnames)
# 
# list_nodes <- lapply(seq_len(length(gatenames)), function(x) random_split(NR, NE))
# names(list_nodes) <- gatenames
# 
# list_experts <- lapply(seq_len(NE), function(x) matrix(dnorm(rnorm(NR))))
# names(list_experts) <- as.character(seq_len(NE))


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


par_to_gate_paths <- function(node, node_type, gate_par_list, X)
{
  if (node_type == "binomial") {
    pars <- gate_par_list[[node]]
    eta <- X %*% pars
    leftpath <- exp(eta) / (1 + exp(eta))
    return(matrix(c(leftpath, 1 - leftpath), ncol=2))
  }
}


par_to_expert_dens <- function(node, expert_type, expert_par_list, Y, X, ...)
{
  if (expert_type == "gaussian") {
    parm <- expert_par_list[[node]]
    beta <- parm[1:(length(parm) - 1)]
    variance <- exp(parm[length(parm)])
    mu <- X %*% beta
    return(dnorm(Y, mean=mu, sd=sqrt(variance), ...))
  }
}


"Input:
       geezer  - the name of a gating node
       youngin - the name of a gating node or expert that is an
                 ancestor of `geezer`
       treestr - list of all gating and expert nodes in the HME
       ln      - list of gating network split probabilities
 Output:
       cumulative product of the paths from the `geezer` node to the
       `youngin` node"
gate_path_product <- function(geezer, youngin, treestr, ln)
{
  "product path from node1 down to node2"
  a1 <- unlist(ancestors(geezer))
  a2 <- unlist(ancestors(youngin))
  if (!geezer %in% a2) {
    stext <- "Node1 [%s] must be an ancestor of node2 [%s]."
    stop(sprintf(stext, node1, node2))
  }
  Reduce(`*`, gate_path_values(setdiff(a2, a1), ln))
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
prior_weights <- function(node, treestr, ln)
{
  "prior weights go from root node down"
  gate_path_product("0", node, treestr, ln)
}


"Calculate the posterior weights for a gating node
 
 Input:
       node    - any given non-terminal node
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
  
  Gs <- napply(terminals, function(x) gate_path_product(node, x, treestr, ln))
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
joint_posterior_weight <- function(node, treestr, lp)
{
  # not feasible for the root node
  if (node == "0") {
    return(rep(1, nrow(lp[[1]])))
  }
  return(gate_path_product("0", node, treestr, ln=lp))
}
