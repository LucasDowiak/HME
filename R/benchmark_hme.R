# ----------------------------------------------------------------------------
setwd("~/hme/R/")
source("building_blocks.R")
source("marginal_effects.R")
source("optimization_functions.R")
source("summary_method.R")
source("plot_method.R")
source("predict_method.R")
source("tree_construction.R")
library(Formula)


hme <- function(tree, formula, hme_type=c("hme", "hmre"),
                expert_type=c("gaussian", "bernoulli"), data, root_prior=1,
                init_gate_pars=NULL, init_expert_pars=NULL,
                maxiter=100, tolerance=1e-4, trace=0)
{
  cl <- match.call()
  require(Formula)
  stopifnot(!missing(data))
  
  tree <- sort(tree)
  hme_type <- match.arg(hme_type)
  expert_type <- match.arg(expert_type)
  
  mf <- match.call(expand.dots = FALSE)
  call_ <- as.list(mf)
  m <- match(c("formula", "data"), names(mf), 0)
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
  
  if (is(init_gate_pars, "NULL")) {
    gate.pars <- napply(gate.nodes, init_gate_node_pars, tree, ncol(Z))
  } else {
    gate.pars <- init_gate_pars
  }
  
  if (is(init_expert_pars, "NULL")) {
    if (expert_type=="gaussian") {
      expert.pars <- bootstrap_glm(length(expert.nodes),
                                   formula=terms(form, data=mf, rhs=1),
                                   data=data,
                                   family=gaussian,
                                   model=FALSE)
      expert.pars <- napply(expert.pars, function(x) c(x, exp(var(Y))))
      names(expert.pars) <- expert.nodes
      #expert.pars <- napply(expert.nodes, function(x) c(runif(ncol(X) , -2, 2),
      #                                                  runif(1, 1, 5)))
    } else {
      expert.pars <- napply(expert.nodes, function(x) runif(ncol(X) , -2, 2))
    }
    
    
  } else {
    expert.pars <- init_expert_pars
  }

  NN <- length(c(unlist(gate.pars), unlist(expert.pars)))
  
  list_nodes <- napply(gate.nodes, par_to_gate_paths, gate.pars, Z)
  list_experts <- napply(expert.nodes, par_to_expert_dens, expert_type,
                          expert.pars, Y, X)

  
  logL <- matrix(NA_real_, nrow=maxiter + 1)
  parM <- matrix(NA_real_, ncol=length(unlist(c(gate.pars, expert.pars))),
                 nrow=maxiter + 1)
  parM[1, ] <- unlist(c(gate.pars, expert.pars))
  logL[1, ] <- newLL <- -Inf
  
  for (ii in seq_len(maxiter)) {
    oldLL <- newLL
    mstep <- m_step(tree, hme_type, expert_type, Y=Y, X=X, Z=Z,
                    exp.pars=expert.pars, gat.pars=gate.pars,
                    root_prior=root_prior)
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
  logL <- logL[!is.na(logL), , drop=FALSE]
  parM <- parM[apply(!is.na(parM), 1, FUN=all), ]
  
  gate.info.matrix <- napply(gate.nodes, multinomial_info_matrix, tree, gate.pars,
                             mstep$list_priors, mstep$list_posteriors, Z,
                             root_prior)
  
  expert.info.matrix <- napply(expert.nodes, expert_info_matrix, expert_type,
                               expert.pars, mstep$list_posteriors, X, Y)
  
  structure(list(tree=tree,
                 hme.type=hme_type,
                 gate.nodes=gate.nodes,
                 gate.pars=gate.pars,
                 gate.pars.nms=colnames(Z),
                 expert.type=expert_type,
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
                 N=length(Y),
                 no.of.pars=NN,
                 gate.info.matrix=gate.info.matrix,
                 expert.info.matrix=expert.info.matrix,
                 call_=call_),
            class="hme")
}


if (FALSE) {
# -1 + Species + Petal.Length + Sepal.Length
"Sepal.Width ~ Petal.Width | Petal.Width + Petal.Length + Sepal.Length + Species"
tst <- hme(tree,
           "Sepal.Width ~ Petal.Width + Petal.Length + Sepal.Length | Petal.Width + Petal.Length + Sepal.Length",
           data=iris, maxiter=200, tolerance = 1e-6, trace=1)


newtree <- grow_the_tree(tst)
newcall <- c(tst$call_[-c(1:2)], newtree['tree'])
newcall$maxiter <- 500
tst2 <- do.call(hme, newcall)

do.call(bootstrap_glm, newcall[names(newcall) %in% names(formals(glm))])

cols <- c("blue", "orange", "green")
with(iris, plot(Petal.Width, Sepal.Width, col=cols[as.integer(iris$Species)]))
for (e in tst$expert.pars) {
  abline(e[1], e[2])
}
points(iris$Sepal.Length, predict(tst, iris))
(lm1 <- lm("Sepal.Width ~ -1 + Species + Species:Petal.Width", data=iris))


dtf <- fread("../data/heloc_dataset_v1.csv")
dtf[, RiskPerformance2 := ifelse(RiskPerformance=="Good", 1, 0)]
formm <- "RiskPerformance2 ~ ExternalRiskEstimate + MSinceOldestTradeOpen + MSinceMostRecentTradeOpen + AverageMInFile + NumSatisfactoryTrades + NumTrades60Ever2DerogPubRec + NumTrades90Ever2DerogPubRec + PercentTradesNeverDelq + MSinceMostRecentDelq"
glm1 <- bootstrap_glm(formm, family=binomial(link="logit"),
                      data=dtf[, c(2:10, ncol(dtf)), with=FALSE])




init_exp <- list(glm1$coefficients + rnorm(np, sd=.25))
debugonce(hme)
tst <- hme(c("0", "0.1", "0.2"),
           "RiskPerformance2 ~ ExternalRiskEstimate | .",
           data=dtf[, c(2:10, ncol(dtf)), with=FALSE],
           expert_type="bernoulli",
           maxiter=200, tolerance = 1e-6, trace=1)

}
