logLik.hme <- function(obj)
{
  GD1 <- expert_lik_contr(obj$expert.nms, obj$list_density, obj$list_priors)
  sum(simplify2array(GD1))
}


criterion <- function(obj, type=c("aic", "bic", "mse"))
{
  type <- match.arg(type)
  K <- obj[["no.of.pars"]]
  N <- obj[["N"]]
  L <- logLik(obj)
  
  
  if (type %in% c("aic", "bic")) {
    if (type == "aic") {
      penalty <- 2
    } else if (type == "bic") {
      penalty <- log(N)
    }
    out <- (penalty * K - 2 * L) / N
  } else if (type == "mse") {
    mses <- obj$MSE[!is.na(obj$MSE)]
    out <- mses[length(mses)]
  }
  return(out)
}


var_test_class <- function()
{
  return(structure(list(VARstat      =NA_real_,
                        alpha        =NA_real_,
                        null         ="H0: N * varHat --> M_{p + q}(. ; lambda**2)",
                        pvalue       =NA_real_,
                        reject_null  =NA,
                        pars         =list(lambda =NA_real_,
                                           varHat =NA_real_,
                                           N      =NA_integer_)),
                   class="VAR-test"))
}


voung_lr_test_class <- function()
{
  return(structure(list(LRstat      =NA_real_,
                        alpha       =NA_real_,
                        null        ="H0: [1 / sqrt(varHat * N)] * LRn --> N(0, 1)",
                        pvalue      =NA_real_,
                        reject_null =NA,
                        pars        =list(LRn=NA_real_)),
                   class="LR-test"))
}


voung_selection <- function(hme1, hme2, var_test_alpha=0.05, model_test_alpha=0.05,
                            adjustedLR=TRUE)
{
  nm1 <- deparse(substitute(hme1))
  nm2 <- deparse(substitute(hme2))
  stopifnot(all(sapply(list(hme1, hme2), inherits, "hme")))
  if (is.null(hme1$full.vcv)) {
    stop(sprintf("The sandwich variance estimator has not been calulcated for %s", nm1))
  }
  if (is.null(hme2$full.vcv)) {
    stop(sprintf("The sandwich variance estimator has not been calulcated for %s", nm2))
  }
  out <- list(variance_test=NULL, voung_lr_test=NULL)
  
  # Likelihood contributions
  GD1 <- expert_lik_contr(hme1$expert.nms, hme1$list_density, hme1$list_priors)
  lik1 <- rowSums(simplify2array(GD1))
  GD2 <- expert_lik_contr(hme2$expert.nms, hme2$list_density, hme2$list_priors)
  lik2 <- rowSums(simplify2array(GD2))
  
  LR <- log(lik1 / lik2)
  LRn <- sum(LR)
  if (adjustedLR) {
    p <- hme1$no.of.pars
    q <- hme2$no.of.pars
    # Uses the Schwarz (1978) correction factor
    adjusted <- (p / 2) * log(hme1$N) - (q / 2) * log(hme1$N)
    LRn <- LRn - adjusted
  }
  
  # 1 Variance Test
  variance_test <- var_test_class()
  variance_test$alpha = model_test_alpha
  
  # First, test that f(.|.,theta1) = g(.|.,theta2) which is equivalent to
  # testing if varHat is significantly different than zero
  
  # Calculate the eigenvalues of matrix W defined in Eq. 3.6: lambda^{2}_{star}
  # See Eq. 4.6
  # hessian = A = H
  # opg     = B = G
  #           ff              
  # W = [
  #       -G1 * H1^{-1}       -OPG12 * H2^{-1} 
  # 
  #       -OPG21 * H1^{-1}    -G2 * H2^{-1}
  # ]
  OPG12 <- mapply(function(x) hme1$scores[x, ] %*% t(hme2$scores[x, ]), seq_len(hme1$N))
  dim(OPG12) <- c(ncol(hme1$scores), ncol(hme2$scores), hme1$N)
  OPG12 <- rowSums(OPG12, dims=2)
  
  W <- rbind(cbind(-hme1$full.vcv$OPG %*% solve(hme1$full.vcv$hessian),
                   -OPG12 %*% solve(hme2$full.vcv$hessian)
                   ),
             cbind(-t(OPG12) %*% solve(hme1$full.vcv$hessian),
                   -hme2$full.vcv$OPG %*% solve(hme2$full.vcv$hessian)
                   )
             )
  # Calculate sample variance Eq 4.2 from Voung
  variance_test$pars$lambda <- eigen(W)$values
  variance_test$pars$N <- hme1$N
  variance_test$pars$varHat <- mean(LR**2) - (mean(LR))**2
  variance_test$VARstat <- variance_test$pars$varHat * variance_test$pars$N
  
  # If varHat = 0 [or equiv f() = g()], then n * varHat is asymptotically
  # distributed as sum of weighted chi-sq ~ M(n * varHat | lambda^{2})
  # find the test_statistic
  
  # If varHat = 0 and information matrix equivalencies of corollary 4.4
  # asy distr of n * varHat is (centered) chi-squared with dof = p + q - 2 * rank(OPG12)
  # dof = length(lambda) - 2 * Matrix::rankMatrix(OPG12)
  
  variance_test$pvalue <- mgcv::psum.chisq(variance_test$VARstat,
                                           lb=variance_test$pars$lambda**2)
  variance_test$reject_null <- variance_test$pvalue < variance_test$alpha
  out$variance_test <- variance_test
  
  if (!variance_test$reject_null) {
    cat(sprintf("Variance Test: Failed to reject the null of zero variance - %s.\n", variance_test$null))
    cat(sprintf("The available data can not discriminate between model %s and model %s\n\n", nm1, nm2))
    return(out)
    
  } else {
    cat(sprintf("Variance Test: Reject null of zero variance - %s\n", variance_test$null))
    direct_text <- "%s is favored over %s."
    voung_lr_test <- voung_lr_test_class()
    voung_lr_test$pars$LRn = LRn
    
    # use Eq. 6.4 and the directional tests of Eqs. 5.7 and 5.8
    voung_lr_test$LRstat <- LRn * (1 / (sqrt(variance_test$pars$varHat * hme1$N)))
    voung_lr_test$alpha <- model_test_alpha
    voung_lr_test$pvalue <- pnorm(voung_lr_test$LRstat)
    reject_null <- voung_lr_test$pvalue < voung_lr_test$alpha / 2 | voung_lr_test$pvalue > 1 - voung_lr_test$alpha / 2
    voung_lr_test$reject_null <- reject_null
    voung_lr_test$pars$varHat <- variance_test$varHat
    
    if (sign(voung_lr_test$LRstat) > 0) {
      cat(sprintf(direct_text, nm1, nm2))
    } else {
      cat(sprintf(direct_text, nm2, nm1))
    }
    out$voung_lr_test <- voung_lr_test
    return(out)
  }
  
  # If strictly non-nested, then f() != g() at all
  
  # Strictly Non-nested
  # Overlapping
  # Nested
}



grab_vitals <- function(obj, expert)
{
  # gate product to expert to pass to root_posterior
  root_prior <- gate_path_product("0", expert, obj$list_priors)
  # standard errors of expert parameters
  stderrs <- sqrt(diag(obj$expert.info.matrix[[expert]][["sandwich"]]))
  # two new expert parameters (from fitted values) (gaussian deviations with std errors)
  pars <- obj[["expert.pars"]][[expert]]
  newpars <- list(rnorm(length(pars), mean=pars, sd=stderrs),
                  rnorm(length(pars), mean=pars, sd=stderrs))
  names(newpars) <- c("0.1", "0.2")
  return(list(root_prior=root_prior, init_expert_pars=newpars))
}


grow_the_tree <- function(obj)
{
  oldtree <- obj[["tree"]]
  call_ <- obj[["call_"]]
  call_[["tree"]] <- c("0", "0.1", "0.2")
  expnms <- obj[["expert.nms"]]
  buds <- vector("list", length(expnms))
  names(buds) <- expnms
  for (ee in expnms) {
    newcall <- c(call_[-1], grab_vitals(obj, ee))
    buds[[ee]] <- do.call(hme, newcall)
  }
  ll <- which.max(sapply(buds, logLik))
  
  new_gate_pars <- buds[[ll]][["gate.pars"]]
  names(new_gate_pars) <- expnms[ll]
  new_gate_pars <- c(obj[["gate.pars"]], new_gate_pars)[sort(c(obj[["gate.nodes"]], expnms[ll]))]
  newtree <- c(obj[["tree"]], paste(names(ll), 1:2, sep="."))
  return(list(tree=newtree, init_gate_pars=new_gate_pars))
}

# Estimate an GLM model, find the fitted values and standard errors, and then
# produce n random draws of the parameter vector 
bootstrap_glm <- function(n=2, ...)
{
  object <- glm(...)
  df.r <- object$df.residual
  dispersion <- if (object$family$family %in% c("poisson", "binomial")) 
    1
  else if (df.r > 0) {
    est.disp <- TRUE
    if (any(object$weights == 0)) 
      warning("observations with zero weight not used for calculating dispersion")
    sum((object$weights * object$residuals^2)[object$weights > 0]) / df.r
  }
  else {
    est.disp <- TRUE
    NaN
  }
  p <- object$rank
  p1 <- 1L:p
  Qr <- stats:::qr.lm(object)
  covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  covmat <- dispersion * covmat.unscaled
  # Randomly draw n set of initial coef
  replicate(n,
            rnorm(p, mean=object$coefficients, sd = sqrt(diag(covmat))),
            simplify=FALSE)
}




