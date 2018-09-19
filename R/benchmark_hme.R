# ----------------------------------------------------------------------------
setwd("hme/R/")
source("building_blocks.R")
source("marginal_effects.R")
source("optimization_functions.R")
library(Formula)


hme <- function(formula, tree, hme_type=c("hme", "hmre"),
                expert_type=c("gaussian"), data,
                maxiter=100, tolerance=1e-4, trace=0)
{
  cl <- match.call()
  require(Formula)
  stopifnot(!missing(data))
  
  tree <- sort(tree)
  hme_type <- match.arg(hme_type)
  
  mf <- match.call(expand.dots = FALSE)
  #m <- match(c("formula", "data"), names(mf), 0)
  form <- Formula::as.Formula(formula)
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  Z <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  
  tree <- sort(tree)
  expert.nodes <- tree[unlist(is_terminal(tree, tree))]
  gate.nodes <- setdiff(tree, expert.nodes)
  
  # This needs to be generalized for HMRE
  expert.index <- expert_index(hme_type, tree)
  
  gate.pars <- napply(gate.nodes, init_gate_node_pars, tree, ncol(Z))
  expert.pars <- napply(expert.nodes, function(x) c(runif(ncol(X) , -2, 2),
                                                    runif(1, 1, 5)))
  NN <- length(c(unlist(gate.pars), unlist(expert.pars)))
  
  list_nodes <- napply(gate.nodes, par_to_gate_paths, gate.pars, Z)
  list_experts <- napply(expert.nodes, par_to_expert_dens, "gaussian",
                          expert.pars, Y, X)

  
  logL <- matrix(NA_real_, nrow=maxiter + 1)
  parM <- matrix(NA_real_, ncol=length(unlist(c(gate.pars, expert.pars))),
                 nrow=maxiter + 1)
  parM[1, ] <- unlist(c(gate.pars, expert.pars))
  newLL <- -Inf
  for (ii in seq_len(maxiter)) {
    oldLL <- newLL
    mstep <- m_step(tree, hme_type, Y=Y, X=X, Z=Z, exp.pars=expert.pars, gat.pars=gate.pars)
    expert.pars <- mstep[["exp.pars"]]
    gate.pars <- mstep[["gat.pars"]]
    newLL <- mstep[["loglik"]]
    logL[ii + 1, ] <- newLL
    parM[ii + 1, ] <- unlist(c(gate.pars, expert.pars))
    if (trace > 0) {
      cat('\r', sprintf("Step: %d - Log-Likelihood: %f", ii, newLL))
    } else {
      cat('\n', sprintf("Step: %d - Log-Likelihood: %f", ii, newLL))
    }
    if (abs(oldLL - newLL) < tolerance)
      break
  }
  cat("\n")
  
  gate.margins <- margin_matrix(expert.nodes, gate.pars, mstep$list_priors)
  gate.info.matrix <- napply(gate.nodes, multinomial_info_matrix, tree, gate.pars,
                             mstep$list_priors, mstep$list_posteriors, Z)
  
  structure(list(tree=tree,
                 gate.nodes=gate.nodes,
                 gate.pars=gate.pars,
                 gate.pars.nms=colnames(Z),
                 expert.nms=expert.nodes,
                 expert.pars=expert.pars,
                 expert.pars.nms=c(colnames(X), "Dispersion"),
                 dependent.var.nme=all.vars(form)[1],
                 logL=logL,
                 parM=parM,
                 list_priors=mstep$list_priors,
                 list_posteriors=mstep$list_posteriors,
                 list_density=mstep$list_density,
                 Y=Y,
                 X=X,
                 Z=Z,
                 no.of.pars=NN,
                 gate.margins=gate.margins,
                 gate.info.matrix=gate.info.matrix),
            class="hme")
}


data(iris)

tree <- c("0",
          "0.1", "0.2",
          "0.1.1", "0.1.2")
tree2 <- c("0",
           "0.1", "0.2", "0.3")
debugonce(hme)

# -1 + Species + Petal.Length + Sepal.Length
"Sepal.Width ~ Petal.Width | Petal.Width + Petal.Length + Sepal.Length"
tst <- hme("Sepal.Width ~ Petal.Width | Petal.Width + Petal.Length + Sepal.Length",
           tree=tree, data=iris, maxiter=200, tolerance = 1e-6, trace=1)
tst2 <- hme("Sepal.Width ~ Petal.Width | -1 + Species + Petal.Width + Petal.Length + Sepal.Length",
           tree=tree2, data=iris, maxiter=250, tolerance=1e-5, trace=1)


cols <- c("blue", "orange", "green")
with(iris, plot(Petal.Width, Sepal.Width, col=cols[as.integer(iris$Species)]))
for (e in tst$expert.pars) {
  abline(e[1], e[2])
}


summary.hme <- function(obj)
{  
  # cat(sprintf("Log-Likelihood: %f", log_Like(obj)))
  for (g in names(obj[["gate.pars"]])) {
    
    pars <- obj[["gate.pars"]][[g]]
    if (!is(pars, "list"))
      pars <- list(pars)
    info <- obj[["gate.info.matrix"]]
    if (!is(info, "list"))
      info <- list(info)
    
    nn <- length(pars)
    cat(sprintf("\n--------------------------\ngate-%s\n", g))
    for (i in seq_len(nn)) {
      p <- pars[[i]]
      se <- sqrt(diag(info[[i]][["sandwich"]]))
      zstat <- p / se
      tbl <- cbind(p, se, zstat, 2 * pnorm(-abs(zstat)))
      colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
      printCoefmat(tbl, digits = 3)
      
      # zscore <- 
      #colnames(p) <- list("Estimate Std.", "Std. Error") # "t-value", "Pr(>|t|)")
      #print(p)
      cat("\n")
    }
  }
  
  cat("\n\n----------------------------------------------------\n\n")
  
  for (e in names(obj[["expert.pars"]])) {
    p <- obj[["expert.pars"]][[e]]

    p <- as.matrix(round(c(p, exp(p[length(p)])), 4))
    dimnames(p) <- list(c(obj[["expert.pars.nms"]], ""), "")
    cat(sprintf("\n--------------------------\nexpert-%s\n", e))
    print(p)
  }
  
  cat("\n\n----------------------------------------------------\n\n")
  
  for (e in names(obj[["expert.pars"]])) {
    p <- obj[["gate.margins"]][[e]]
    
    p <- as.matrix(round(p, 4))
    dimnames(p) <- list(c(obj[["gate.pars.nms"]]), "")
    cat(sprintf("\n--------------------------\nexpert-%s marginal effects\n", e))
    print(p)
  }
}

debugonce(summary.hme)
summary(tst)




