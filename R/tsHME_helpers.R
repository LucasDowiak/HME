node_names <- function(gn)
{
  squish <- function(parent, children)
  {
    out <- c()
    for (p in parent) {
      out <- c(out, paste(p, seq_len(children), sep="."))
    }
    out
  }
  unlist(Reduce(squish, gn, init="0", accumulate=T))
}


multinomial_g <- function(X, par, idx, gn, ss, depth, initsplit)
{
  gpn <- ncol(X)
  matGb <- matGe <- NULL
  initGb <- lapply(list(par[idx[[1]]]), gating_weights, XX=X, splits=initsplit)
  if (depth > 1) {
      bibool <- sapply(idx[-1], function(x) length(x) == gpn)
      lstpar <- lapply(idx, function(x) par[x])
      matGb <- lapply(lstpar[-1][bibool], gating_weights, XX=X, splits=2)
      matGe <- lapply(lstpar[-1][!bibool], gating_weights, XX=X, splits=ss)
  }
  out <- c(initGb, matGb, matGe)
  names(out) <- if (depth > 1) {
    node_names(gn)
  } else {
    "0"
  }
  out
}


gating_weights <- function(beta, XX, splits)
{
  gpn <- ncol(XX)
  eps <- .Machine$double.eps^(2/3)
  break_beta <- function(x, n)
  {
    ii <- seq_len(length(x))
    lapply(split(ii, ceiling(ii / n)), function(i) x[i])
  }
  G <- matrix(NA_real_, nrow(XX), splits-1)
  lstbeta <- break_beta(beta, gpn)
  for (b in seq_along(lstbeta)) {
    G[,b] <- exp(XX %*% lstbeta[[b]])
  }
  G <- cbind(G, 1)
  G <- G / rowSums(G)
  G[G < eps] <- eps; G[G > 1-eps] <- 1-eps
  return(G)
}


lstcumsum <- function(x, y) x[length(x)] + cumsum(y)


get_n_par <- function(x)
{
  if (x[[1]] == 0) {
      x[[1]] <- integer(0)
    }
    length(x[[1]]) + x[[2]] + x[[3]] + 1
}


expert_index <- function(ar, p)
{
  np <- sapply(ar, get_n_par)
  idx <- mapply(rep_len, 1, length.out=np, SIMPLIFY=FALSE)
  Reduce(lstcumsum, idx, init=p, accumulate=TRUE)[-1]
}

# linear model
expert_index_lm <- function(nexp, eparn, gpn)
{
  idx <- mapply(rep_len, 1, length.out=rep(eparn, nexp), SIMPLIFY=FALSE)
  Reduce(lstcumsum, idx, init=gpn, accumulate=TRUE)[-1]
}
# linear model



gate_index <- function(gn, D, ss, nc)
{
  nones <- if (D==1) {
    (ss - 1) * nc
  } else {
    (c(gn, ss) - 1) * nc
  }
  nlists <- cumprod(c(1, gn))
  idx <- c()
  for (j in seq_along(nones)) {
    idx <- c(idx, replicate(nlists[j], rep_len(1, nones[j]), simplify = FALSE))
  }
  Reduce(lstcumsum, idx, init=0, accumulate=TRUE)[-1]
}


expert_names <- function(x)
{
  as.character(seq_len(x))
}


normalizer <- function(W, D)
{
  out <- W * D
  norm <- rowSums(out)
  list(H=sweep(out, 1, norm, `/`), norm=norm)
}



e_step <- function(WW, DD, gp, nexp_call, nexp)
{
  eps <- .Machine$double.eps**(2/3)
  if (is.null(nexp_call)) {                                     # Modified for lm version
    grpidx <- split(seq_len(nexp), ceiling(seq_len(nexp) / 2))  # allows for each terminal node to be an expert
    DD <- lapply(grpidx, function(i) DD[, i])
  }
  normalize_level <- function(lst)
  {
    norm <- sapply(lst, `[[`, "norm")
    childs <- colnames(norm)
    parents <- unique(sub("\\.[^\\.]*$", "", childs))
    normalz <- function(p)
    {
      pat <- paste0("^", p)
      normalizer(WW[[p]], norm[, grepl(pat, childs)])
    }
    out <- lapply(parents, normalz)
    names(out) <- parents
    out
  }
  chain_priors <- function(lst)
  {
    p <- names(lst)
    pat <- sprintf("^%s\\.[1-9]$", gsub("\\.", "\\\\\\.", p))
    out <- c()
    for (np in seq_along(pat)) {
      bottom <- WW[grepl(pat[np], names(WW))]
      nn <- length(bottom)
      nc <- ncol(lst[[np]])
      chain <- lapply(seq_len(nn), function(x) {
        sweep(bottom[[x]], 1, lst[[np]][,x], `*`)
      })
      names(chain) <- names(bottom)
      out <- c(out, chain)
    }
    out
  }
  node.per.lev <- cumprod(c(1, gp))
  hout <- vector("list", length(gp) + 1)
  if (is.null(nexp_call)) {
    hout[[length(gp) + 1]] <- mapply(normalizer, tail(WW, node.per.lev[length(node.per.lev)]),
      DD, SIMPLIFY=FALSE)
  } else {
    hout[[length(gp) + 1]] <- lapply(tail(WW, node.per.lev[length(node.per.lev)]),
                                   normalizer, D=DD)  
  }
  for (lev in rev(seq_along(gp))) {
    hout[[lev]] <- normalize_level(hout[[lev + 1]])
  }
  for (lev in seq_along(gp)) {
    if (lev == 1) {
      gout <- chain_priors(WW["0"])
    } else {
      gout <- chain_priors(gout)
    }
  }
  if (length(WW)==1) {
    hout <- hout[length(hout)]
    gout <- WW
  }
  H <- lapply(unlist(hout, recursive=FALSE), `[[`, "H")
  if (is.null(nexp_call)) {                             # Modified for the lm version
    GD <- mapply(`*`, gout, DD, SIMPLIFY=FALSE)         # allows each terminal node to be an expert
    GD <- do.call(cbind, GD)
  } else {
    GD <- lapply(gout, function(g) g * DD)
    GD <- Reduce(`+`, GD, init=0)
  }
  return(list(GD=GD, G=gout, posterior=H, prior=WW))
}


prior_gate_idx <- function(d)
{
  if (d <= 1)
    return(0)
  g <- 2**d - 1
  c(0, rep(seq(1, (d-1)/2), each=2))
}



m_step <- function(parm, Y, XG, AR, estep, one, ng, nexp, maxit)
{
  estimate_gates <- function(node, g)
  {
    if (g > 2) {
      multinom <- irls_multinomial
      ss <- g
    } else {
      multinom <- irls_logit
      ss <- NULL
    }
    if (node != "0") {
      ns <- strsplit(node, "\\.")[[1]]
      cc <- as.numeric(ns[length(ns)])
      p <- paste(ns[-length(ns)], collapse=".")
      return(multinom(XG, estep$posterior[[node]], estep$posterior[[p]][,cc],
                      maxit=maxit))
    } else if (node == "0") {
      return(multinom(XG, estep$posterior[[node]], maxit=maxit))
    }
  }
  gates <- Map(estimate_gates, names(estep$posterior),
               g=rep(c(ng, nexp), cumprod(c(1, ng))))
  if (ng==1) {
    gates <- Map(estimate_gates, names(estep$posterior), g=nexp)
  } else if (ng > 1) {
    gates <- Map(estimate_gates, names(estep$posterior),
                 g=rep(c(ng, nexp), cumprod(c(1, ng))))
  }
  expert_wts <- sweep(estep$GD, 1, rowSums(estep$GD), `/`)
  experts <- vector("list", one)
  cnams <- colnames(expert_wts)
  for (e in unique(cnams)) {
      cidx <- which(cnams == e)
      experts[[as.integer(e)]] <- cll_ar(Y, AR[[cidx]], expert_wts[, cidx])
  }
  names(experts) <- as.character(rep(1:one))
  return(list(gates=gates, experts=experts))
}


m_step_lm <- function(parm, Y, XG, XE, estep, one, ng, ss, maxit)
{
  estimate_gates <- function(node, g)
  {
    if (g > 2) {
      multinom <- irls_multinomial
      ss <- g
    } else {
      multinom <- irls_logit
      ss <- NULL
    }
    if (node != "0") {
      ns <- strsplit(node, "\\.")[[1]]
      cc <- as.numeric(ns[length(ns)])
      p <- paste(ns[-length(ns)], collapse=".")
      return(multinom(XG, estep$posterior[[node]], estep$posterior[[p]][,cc],
                      maxit=maxit))
    } else if (node == "0") {
      return(multinom(XG, estep$posterior[[node]], maxit=maxit))
    }
  }
  gates <- Map(estimate_gates, names(estep$posterior),
               g=rep(c(ng, ss), cumprod(c(1, ng))))
  #if (ng==1) {
  #  gates <- Map(estimate_gates, names(estep$posterior), g=nexp)
  #} else if (ng > 1) {
  #  gates <- Map(estimate_gates, names(estep$posterior),
  #               g=rep(c(ng, nexp), cumprod(c(1, ng))))
  #}
  expert_wts <- sweep(estep$GD, 1, rowSums(estep$GD), `/`)
  experts <- vector("list", one)
  cnams <- colnames(expert_wts)
  for (e in unique(cnams)) {
      cidx <- which(cnams == e)
      experts[[as.integer(e)]] <- init_lm(Y, XE, expert_wts[, cidx])
  }
  names(experts) <- as.character(rep(1:one))
  return(list(gates=gates, experts=experts))
}


irls_logit <- function(X, H, w=NULL, maxit=25, tol=1e-08)
{
  H <- H[,1]
  family <- binomial
  WNULL <- is.null(w)
  ss <- t <- 0
  QR <- qr(X)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  for(iter in 1:maxit) {
      g <- family()$linkinv(t)
      gprime <- family()$mu.eta(t)
      z <- t + (H - g) / gprime
      if (WNULL) {
        A <- 1
      } else {
        A <- w
      }
      W <- as.vector((A / family()$variance(g)) * gprime^2)
      wmin <- min(W)
      if(wmin < sqrt(.Machine$double.eps))
        warning("Tiny weights encountered")
      s_old <- ss
      C <- chol(crossprod(Q, W*Q))
      ss <- forwardsolve(t(C), crossprod(Q,W*z))
      ss <- backsolve(C,ss)
      t <- Q %*% ss
      if(sqrt(crossprod(ss - s_old)) < tol)
        break
    }
  x <- backsolve(R, crossprod(Q,t))
  list(parm=x, iterations=iter)
}


irls_multinomial <- function(X, H, w=NULL, maxit=25, tol=1e-08)
{
  # X can be a single matrix (same model but different parameters) or
  nc <- ncol(H)
  lst <- vector("list", nc-1)
  for (cc in seq_len(nc-1)) {
    bM <- H[,c(cc, nc)]
    bM <- sweep(bM, 1, rowSums(bM), `/`)
    # if X is list X[[g]] here
    lst[[cc]] <- irls_logit(X, bM, w, maxit=maxit, tol=tol)
  }
  list(parm=lapply(lst, `[[`, "parm"), iterations=lapply(lst, `[[`, "iterations")) 
}


Xpdf_ar1 <- function(parm, y, x, wt, hess=F)
{
  stopifnot(length(y) == nrow(x),
            length(y) == length(wt))
  pdf_fun <- function(parm, y, x, wt) {
    np <- length(parm)
    TT <- as.integer(length(y))
    sigma2 <- parm[np]; phi <- parm[2]; cc <- parm[1]
    SS <- (y - tcrossprod(x, matrix(parm[1:2], nrow=1)))^2 / (2*sigma2)
    #SS <- (1 / (2*sigma2)) *
    #  c((1 - phi**2) * (y[1] - (cc / (1 - phi)))**2, eps[-1]**2)
    DET <- (1/2) * (log(1 - phi^2) - TT * log(sigma2))
    return(-sum(wt * ((-TT/2) * log(2*pi) + DET + SS))) # Log-Likelihood
  }
  S <- 1e-5
  return(optim(parm, pdf_fun, method="L-BFGS-B", y=y, x=x, hessian=hess, wt=wt,
               lower=c(-Inf, -1 + S, S), upper=c(Inf, 1 - S, Inf)))
}


get_density <- function(lsfitmod)
{
  eps <- .Machine$double.eps^(2/3)
  p <- lsfitmod$parm
  DD <- dnorm(lsfitmod$wlsfit$residuals, sd=sqrt(p[length(p)]))
  DD[DD < eps] <- eps
  return(DD)
}


Gaussian <- function(par, y, X, ar=NULL)
{
  eps^(2/3)
  n <- length(par)
  nr <- nrow(y)
  nc <- ncol(y)
  if (ar[[4]]) {
    Y <- matrix(NA_real_, nr, nc-1)
    for (cc in 1:(nc-1)) {
      Y[,cc] <- y[,cc] - y[,cc+1]
    }
    y <- Y
    if (n == 1) {
      DD <- dnorm(y[,1], sd=sqrt(par))
    } else {
      DD <- dnorm(y[,1] - crossprod(t(X), par[-n]), sd=sqrt(par[n]))
    }
  } else {
    DD <- dnorm(y[,1] - crossprod(t(X), par[-n]), sd=sqrt(par[n]))
  }
  DD[DD < eps] <- eps
  return(DD)
}


get_coef_names <- function(ar)
{
  arnames <- if (all(ar[[1]] > 0)) {
    paste0("ar", ar[[1]])
  }
  c(if (ar[[2]]) "const", if (ar[[3]]) "trend",  arnames)
}


X_from_Y <- function(y, ar)
{
  nr <- nrow(y)
  nc <- ncol(y)
  if (ar[[4]]) {
    Y <- matrix(NA_real_, nr, nc-1)
    for (cc in 1:(nc-1)) {
      Y[,cc] <- y[,cc] - y[,cc+1]
    }
    y <- Y
  }
  x <- if (max(ar[[1]]) > 0) {
    x <- y[, ar[[1]] + 1, drop=FALSE]
    colnames(x) <- paste0("ar", ar[[1]])
    x
  }
  con <- if (ar[[2]]) rep(1, nr)
  trend <- if (ar[[3]]) seq_len(nr)/nr
  return(cbind(con, trend, x))
}


cll_ar <- function(y, ar, wt=NULL)
{
  if (is.null(wt))
    wt <- rep(1, nrow(y))
  XX <- X_from_Y(y, ar)
  yy <- if (!ar[[4]]) {
    y[,1]
  } else {
    y[,1] - y[,2]
  }
  if (is.null(XX)) {
    sigma2 <- if (is.null(wt)) {
      sum(yy^2)
    } else {
      V1 <- sum(wt)
      V2 <- sum(wt^2)
      as.numeric(crossprod(wt, yy^2) / (V1 - (V2/V1)))
    }
    return(list(parm=sigma2, wlsfit=list(residuals=yy)))
  } else {
    # bias correction
    wls <- lsfit(XX, yy, wt, intercept=FALSE)
    sigma2 <- if (is.null(wt)) {
      sum(wls$residuals^2) / (nrow(XX) - ncol(XX))
    } else {
      V1 <- sum(wt)
      V2 <- sum(wt^2)
      as.numeric(crossprod(wt, wls$residuals^2)) / (V1 - (V2/V1))
    }
    return(list(parm=c(wls$coef, sigma2), wlsfit=wls[-1]))
  }
}


init_lm <- function(y, x, wt=NULL)
{
  if (is.null(wt))
    wt <- rep(1, nrow(y))
  wls <- lsfit(x, y, wt, intercept=FALSE)
  sigma2 <- if (is.null(wt)) {
      sum(wls$residuals^2) / (nrow(XX) - ncol(XX))
    } else {
      V1 <- sum(wt)
      V2 <- sum(wt^2)
      as.numeric(crossprod(wt, wls$residuals^2)) / (V1 - (V2/V1))
    }
  return(list(parm=c(wls$coef, sigma2), wlsfit=wls[-1]))
}


arCheck <- function(ar)
{
  p <- max(which(c(1, -ar) != 0)) - 1
  if (!p) 
    return(TRUE)
  all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
}


starting_weights <- function(y, ar)
{
  nparts <- length(ar)
  nr <- nrow(y)
  seqs <- seq_len(nr)
  regimes <- split(seqs, ceiling(seqs / (nr/nparts)))
  wts <- replicate(nparts, rep_len(0, nr), simplify=FALSE)
  for (ii in seq_along(wts)) {
    wts[[ii]][regimes[[ii]]] <- 1
  }
  return(wts)
}


calc_bias <- function(lst, PLOT=TRUE)
{
  idx <- lst[[1]]$idx
  SSs <- unique(unlist(lapply(lst, `[[`, "ss")))
  PAR <- c()
  for (ss in SSs) {
    lstss <- Filter(function(x) x$ss==ss && length(x) > 1, lst)
    lstCoef <- vector("list", length(lstss))
    for (jj in seq_along(lstss)) {
      par <- lstss[[jj]]$par
      lstCoef[[jj]] <- rbind(c(ss, par[idx[[1]]]),
                             c(ss, par[idx[[2]]]))
    }
    PAR <- c(PAR, list(ss=do.call(rbind, lstCoef)))
  }
  names(PAR) <- SSs
  BETA <- c(0.8,0.2,0.5,0.1)
  find_bias <- function(P)
  {
    cls <- kmeans(P[,-1], 2)
    # generalize this 
    h <- which.max(cls$centers[,1])
    l <- which.min(cls$centers[,1])
    tbl <- table(cls$cluster)
    if (!all(tbl==tbl[1]))
      warnings("Sample size %d did not have equal cluster size.", P[1,1])
    return(c(cls$centers[h,], cls$centers[l,]))
  }
  matBias <- sweep(sapply(PAR, find_bias), 1, BETA, `-`)
  if (PLOT) {
    ss <- as.integer(colnames(matBias))
    tits <- c(expression(displaystyle(phi[1])), expression(displaystyle(sigma[1]^{2})),
              expression(displaystyle(phi[2])), expression(displaystyle(sigma[2]^{2})))
    par(mfrow=c(2,2))
    for (pp in 1:nrow(matBias)) {
      b <- matBias[pp, ]
      rb <- max(abs(range(b)))
      pd <- 0.10 * rb
      plot(ss, b, type="l", ylim=c(c(-pd,pd) + c(-rb,rb)), main=tits[pp],
           ylab="bias", xlab="sample size")
      abline(h=0, col="red"); grid()
    }
  }
  return(matBias)
}


get_llh <- function(GD)
{
  if (is.matrix(GD)) {
    return(sum(log(rowSums(estep$GD))))
  } else if (is.list(GD)) {
    return(sum(log(unlist(GD))))
  } else {
    stop("Case not met.")
  }

}

log_Like <- function(obj)
{
  sum(log(rowSums(obj$weights$GD)))
}


logLik.tsHME <- function(obj)
{
  obj$loglik[length(obj$loglik)]
}


coef.tsHME <- function(obj)
{
  obj$par
}


coef.lmHME <- function(obj)
{
  obj$par
}


summary.tsHME <- function(obj)
{
  theta <- obj$par
  ngp <- obj$num.gate.par
  nep <- obj$num.expert.par
  lstgatpar <- lapply(obj$par.gate.index, function(x) theta[x])
  lstexppar <- lapply(obj$par.expert.index, function(x) theta[x])
  sum_gate_par <- function(par, ng)
  {
    lstpar <- split(seq_along(par), ceiling(seq_along(par) / ng))
    for (i in seq_along(lstpar)) {
      p <- as.matrix(par[lstpar[[i]]])
      dimnames(p) <- list(paste(sprintf("w%d_", i), 0:(ng-1), sep=""), "")
      print(p)
      cat("\n")
    }
  }
  cat(sprintf("Log-Likelihood: %f", log_Like(obj)))
  for (i in seq_len(obj$num.gates)) {
    p <- round(lstgatpar[[i]], 4)
    # rownames(p) <- paste("w", 0:(ngp-1), sep=""); colnames(p) <- ""
    cat(sprintf("\n--------------------------\ngate-%d\n", i))
    sum_gate_par(p, ngp)
  }
  cat("\n\n--------------------------\n--------------------------\n\n")
  for (j in seq_len(obj$num.experts)) {
    p <- as.matrix(round(lstexppar[[j]], 6))
    arspec <- obj$AR[[j]]
    rnames <- c(get_coef_names(obj$AR[[j]]), "var")
    rownames(p) <- rnames; colnames(p) <- ""
    cat(sprintf("\n--------------------------\nexpert-%d\n", j))
    print(p)
  }
}


plot.tsHME <- function(obj, type=c("state_probs", "network", "prior",
  "posterior", "par.conv.gate", "par.conv.expert", "wts.scatter"))
{
  type <- match.arg(type)
  if (type == "network") {
    cols <- c("steelblue", "tomato2", "orange", "olivedrab4")
    GG <- multinomial_g(obj$model.matrix.gate, obj$par, obj$par.gate.index,
      obj$num.gate.par, obj$num.experts, obj$depth, obj$first_split)
    boxshape <- c(obj$depth, max(obj$level.node))
    #lascol <- (1:boxshape[1]) * (boxshape[2]/boxshape[1])
    par(mfrow=boxshape)
    jj <- 1
    for (d in 1:boxshape[1]) {
      pltidx <- obj$level.node[d]
      for (cc in 1:boxshape[2]) {
        if (cc <= pltidx) {
          plot(GG[[jj]][, 1], type="l", ylim=c(0,1), lwd=2, col=cols[1])
          for (ll in 2:ncol(GG[[jj]])) {
            lines(GG[[jj]][,ll], col=cols[ll], ylim=c(0,1), lwd=2)
          }
          grid()
          jj <- jj + 1
        } else {
          frame()
        }
      } 
    }
  } else if (type == "state_probs")
  {
    ne <- obj$num.experts
    G <- if (is.null(obj$callNexp)) {
      do.call(cbind, obj$weights$G)
    } else {
      Reduce(`+`, obj$weights$G, init=0)
    }
    post <- sweep(obj$weights$GD, 1, rowSums(obj$weights$GD), `/`)
    par(bty="l", mfrow=c(ne, 1), mar=rep(0,4), oma=c(4,4,3,1))
    for (ii in seq_len(ne)) {
      if (ii == ne) {
        plot(post[,ii], type="l", las=1, cex.axis=0.75, ylim=c(0,1))
        lines(G[,ii], col="red"); grid()
      } else {
        plot(post[,ii], type="l", las=1, cex.axis=0.75, ylim=c(0,1), xaxt="n")
        lines(G[,ii], col="red"); grid()
      }
    }
  }
  else if (type == "par.conv.expert")
  {
    ne <- obj$num.experts
    nep <- obj$num.expert.par
    colidx <- seq_len(ne * nep) +
      (obj$num.gates * obj$num.gate.par)
    M <- obj$par.evolution[, colidx]
    par(bty="l", mfcol=c(nep, ne), mar=rep(0,4), oma=c(4,4,3,1))
    for (ii in seq_along(colidx)) {
      if (ii %% nep == 0) {
        plot(M[,ii], type="l", las=1, cex.axis=0.75)
      } else {
        plot(M[,ii], type="l", las=1, cex.axis=0.75, xaxt="n")
      }
    }
    mtext(paste0("b", (nep-1):0), 2, outer=TRUE, line=2, at=seq_len(nep)/(nep+1))
    mtext("Iteration", 1, outer=TRUE, line=2.5)
    mtext("Expert Network Parameters", 3, out=TRUE, line=1.25)
  }
  else if (type == "par.conv.gate")
  {
    ng <- obj$num.gates
    ngp <- obj$num.gate.par
    colidx <- seq_len(ng * ngp)
    M <- obj$par.evolution[, colidx]
    par(bty="l", mfcol=c(ngp, ng), mar=rep(0,4), oma=c(4,4,3,1))
    for (ii in seq_along(colidx)) {
      if (ii %% ngp == 0) {
        plot(M[,ii], type="l", las=1, cex.axis=0.75)
      } else {
        plot(M[,ii], type="l", las=1, cex.axis=0.75, xaxt="n") 
      }
    }
    mtext(paste0("b", (ngp-1):0), 2, outer=TRUE, line=2, at=seq_len(ngp)/(ngp+1))
    mtext("Iteration", 1, outer=TRUE, line=2.5)
    mtext("Gate Network Parameters", 3, out=TRUE, line=1.25)
  }
  if (obj$depth == 2)
  {
    if (type == "wts.scatter")
    {
      # bring in wts
      matpst <- obj$weights$posterior
      matpri <- obj$weights$prior
      ep <- 1/16
      pad <- matrix(rep(c(-1,1,-1,1), each=3), nrow=3)
      PM <- rbind(c(c(3, 5, 5, 7)),
                  c(c(1, 3, 1, 3)),
                  c(c(5, 7, 1, 3))) / 8
      PM <- PM + (pad * ep)
      split.screen(PM)
      for (i in 1:3) {
        screen(i)
        par(mar=rep(3*ep, 4), bty="l", cex.axis=0.8)
        plot(matpst[,i], matpri[,i], ylim=c(0,1), xlim=c(0,1), col="red",
             las=1, pch=2, cex=0.8)
        axis(1, at=1, labels="")
        axis(2, at=1)
        abline(0,1)
      }
      close.screen(all.screens = TRUE)
    } else if (type %in% c("prior", "posterior"))
    {
      M <- obj$weights[[type]]
      ep <- 1/16
      pad <- matrix(rep(c(-1,1,-1,1), each=3), nrow=3)
      PM <- rbind(c(c(3, 5, 5, 7)),
                  c(c(1, 3, 1, 3)),
                  c(c(5, 7, 1, 3))) / 8
      PM <- PM + (pad * ep)
      split.screen(PM)
      for (i in 1:3) {
        histO <- hist(M[,i], plot=F)
        screen(i)
        par(mar=rep(ep, 4))
        hist(M[,i], main="", xlab=NULL, ylab=NULL, freq=FALSE, axes=F)
        axis(1, labels=rep("",2), at=c(0L, 1L))
        axis(2,labels=FALSE, tick=FALSE, at=c(0,round(max(histO$density), 1)))
      }
      close.screen(all.screens = TRUE)
    }
  } else {
    message(sprintf("Plot type %s not available depth %d.", type, obj$depth))
  }
}


plot.lmHME <- function(obj, type=c("state_probs", "network", "prior",
  "posterior", "par.conv.gate", "par.conv.expert", "wts.scatter"))
{
  type <- match.arg(type)
  if (type == "network") {
    cols <- c("steelblue", "tomato2", "orange", "olivedrab4")
    GG <- multinomial_g(obj$model.matrix.gate, obj$par, obj$par.gate.index,
      obj$num.gate.par, obj$num.experts, obj$depth, obj$first_split)
    boxshape <- c(obj$depth, max(obj$level.node))
    #lascol <- (1:boxshape[1]) * (boxshape[2]/boxshape[1])
    par(mfrow=boxshape)
    jj <- 1
    for (d in 1:boxshape[1]) {
      pltidx <- obj$level.node[d]
      for (cc in 1:boxshape[2]) {
        if (cc <= pltidx) {
          plot(GG[[jj]][, 1], type="l", ylim=c(0,1), lwd=2, col=cols[1])
          for (ll in 2:ncol(GG[[jj]])) {
            lines(GG[[jj]][,ll], col=cols[ll], ylim=c(0,1), lwd=2)
          }
          grid()
          jj <- jj + 1
        } else {
          frame()
        }
      } 
    }
  } else if (type == "state_probs")
  {
    ne <- obj$num.experts
    G <- if (is.null(obj$callNexp)) {
      do.call(cbind, obj$weights$G)
    } else {
      Reduce(`+`, obj$weights$G, init=0)
    }
    post <- sweep(obj$weights$GD, 1, rowSums(obj$weights$GD), `/`)
    par(bty="l", mfrow=c(ne, 1), mar=rep(0,4), oma=c(4,4,3,1))
    for (ii in seq_len(ne)) {
      if (ii == ne) {
        plot(post[,ii], type="l", las=1, cex.axis=0.75, ylim=c(0,1))
        lines(G[,ii], col="red"); grid()
      } else {
        plot(post[,ii], type="l", las=1, cex.axis=0.75, ylim=c(0,1), xaxt="n")
        lines(G[,ii], col="red"); grid()
      }
    }
  }
  else if (type == "par.conv.expert")
  {
    ne <- obj$num.experts
    nep <- obj$num.expert.par
    idx <- unlist(obj$par.expert.index)
    maxllh <- which.max(obj$loglik)
    M <- obj$par.evolution[, idx]
    par(bty="l", mfcol=c(nep, ne), mar=rep(0,4), oma=c(4,4,3,1))
    for (ii in seq_along(idx)) {
      if (ii %% nep == 0) {
        plot(M[,ii], type="l", las=1, cex.axis=0.75)
      } else {
        plot(M[,ii], type="l", las=1, cex.axis=0.75, xaxt="n")
      }
      abline(v=maxllh, col="red")
    }
    mtext(paste0("b", (nep-1):0), 2, outer=TRUE, line=2, at=seq_len(nep)/(nep+1))
    mtext("Iteration", 1, outer=TRUE, line=2.5)
    mtext("Expert Network Parameters", 3, out=TRUE, line=1.25)
  }
  else if (type == "par.conv.gate")
  {
    #ng <- obj$num.gates
    ngp <- obj$num.gate.par
    idx <- obj$par.gate.index
    M <- obj$par.evolution
    maxllh <- which.max(obj$loglik)
    plot_gate <- function(i, ngp, title)
    {
      nc <- length(i) / ngp
      par(bty="l", mfcol=c(ngp, nc), mar=rep(0,4), oma=c(4,4,3,1))
      for (ii in i) {
        if (ii %% ngp == 1) {
          plot(M[,ii], type="l", las=1, cex.axis=0.75)
        } else {
          plot(M[,ii], type="l", las=1, cex.axis=0.75, xaxt="n") 
        }
        abline(v=maxllh, col="red")
      }
      mtext(paste0("b", (ngp-1):0), 2, outer=TRUE, line=2, at=seq_len(ngp)/(ngp+1))
      mtext("Iteration", 1, outer=TRUE, line=2.5)
      mtext(paste("Gate Network Parameters", title, sep=": "), 3, out=TRUE, line=1.25)
    }
    mapply(plot_gate, i=idx, ngp=list(ngp), title=names(obj$gates))
  }
  if (obj$depth == 2)
  {
    if (type == "wts.scatter")
    {
      # bring in wts
      matpst <- obj$weights$posterior
      matpri <- obj$weights$prior
      ep <- 1/16
      pad <- matrix(rep(c(-1,1,-1,1), each=3), nrow=3)
      PM <- rbind(c(c(3, 5, 5, 7)),
                  c(c(1, 3, 1, 3)),
                  c(c(5, 7, 1, 3))) / 8
      PM <- PM + (pad * ep)
      split.screen(PM)
      for (i in 1:3) {
        screen(i)
        par(mar=rep(3*ep, 4), bty="l", cex.axis=0.8)
        plot(matpst[,i], matpri[,i], ylim=c(0,1), xlim=c(0,1), col="red",
             las=1, pch=2, cex=0.8)
        axis(1, at=1, labels="")
        axis(2, at=1)
        abline(0,1)
      }
      close.screen(all.screens = TRUE)
    } else if (type %in% c("prior", "posterior"))
    {
      M <- obj$weights[[type]]
      ep <- 1/16
      pad <- matrix(rep(c(-1,1,-1,1), each=3), nrow=3)
      PM <- rbind(c(c(3, 5, 5, 7)),
                  c(c(1, 3, 1, 3)),
                  c(c(5, 7, 1, 3))) / 8
      PM <- PM + (pad * ep)
      split.screen(PM)
      for (i in 1:3) {
        histO <- hist(M[,i], plot=F)
        screen(i)
        par(mar=rep(ep, 4))
        hist(M[,i], main="", xlab=NULL, ylab=NULL, freq=FALSE, axes=F)
        axis(1, labels=rep("",2), at=c(0L, 1L))
        axis(2,labels=FALSE, tick=FALSE, at=c(0,round(max(histO$density), 1)))
      }
      close.screen(all.screens = TRUE)
    }
  } else {
    message(sprintf("Plot type %s not available depth %d.", type, obj$depth))
  }
}


internal_DM <- function(obj, lastN=25)
{
  niter <- nrow(obj$par.evolution)
  theta_ <- obj$par.evolution[niter,]
  e_step <- switch(as.character(obj$depth), "1"=e_step1, "2"=e_step2, "3"=e_step3)
  sub_one <- function(i, par)
  {
    out <- theta_
    out[i] <- par[i]
    return(out)
  }
  map <- function(par)
  {
    Norm <- function(x) Gaussian(par[x], obj$Y, obj$model.matrix.expert)
    Gate <- function(x) gating_weights(obj$model.matrix.gate, par[x])
    matGs <- matrix(unlist(lapply(obj$par.gate.index, Gate)), ncol=obj$num.gates)
    colnames(matGs) <- if (obj$depth > 1) {
      c("0", unlist(lapply(1:(obj$depth - 1), node_names)))
    } else {
      "0"
    }
    matFs <- sapply(obj$par.expert.index, Norm)
    colnames(matFs) <- node_names(obj$depth)
    estep <- e_step(matGs, matFs)
    return(unlist(m_step(par, obj$Y, obj$model.matrix.gate, obj$model.matrix.expert, estep, obj$depth)))
  }
  R <- array(NA_real_, dim=c(rep(length(theta_), 2), lastN))
  for (k in seq_len(lastN)) {
    kk <- niter - (lastN + 1) + k
    theta <- obj$par.evolution[kk,]
    for (i in seq_along(theta)) {
      for (j in seq_along(theta)) {
        R[i,j,k] <- (map(sub_one(i, theta_))[j] - theta_[j]) / (theta[i] - theta_[i]) 
      }
    }
  }
  return(R)
}


score_function <- function(obj)
{

  one <- obj$num.experts
  ong <- obj$num.gates
  g <- obj$weights$prior
  h <- obj$weights$posterior
  tx <- c(0, rep(seq(1, (ong - 1) / 2), each=2))
  eps <- h - g
  gscores <- vector("list", one + ong)
  for (ii in seq_len(ong)) {
    if (ii == 1) {
      gscores[[ii]] <- crossprod(obj$model.matrix.gate, eps[,1,drop=F])
    } else if (ii %% 2 == 0) {
      gscores[[ii]] <- crossprod(obj$model.matrix.gate, diag(h[,tx[ii]]) %*% eps[,ii,drop=F])
    } else {
      gscores[[ii]] <- crossprod(obj$model.matrix.gate, diag(1-h[,tx[ii]]) %*% eps[,ii,drop=F])
    }
  }
  for (ii in seq_len(one)) {
    exprt <- obj$experts[[ii]]; sigma <- exprt$parm[length(exprt$parm)]
    gscores[[ii + ong]] <- crossprod(obj$model.matrix.expert,
                                     -(exprt$wlsfit$wt/sigma) * exprt$wlsfit$residuals)
  }
  return(unlist(gscores))
}


hess_function <- function(obj)
{
  one <- obj$num.experts
  ong <- obj$num.gates
  g <- obj$weights$prior
  h <- obj$weights$posterior
  tx <- c(NA, rep((ong - 1) / 2, each=2))
  eps <- h - g
  hess <- vector("list", one + ong)
  for (ii in seq_len(ong)) {
    if (ii == 1) {
      # Need the W matrix from pg 186 of Venables and Ripley
      W <- 1 / binomial()$variance(g[,ii])
    } else if (ii %% 2 == 0) {
      W <- h[,tx[ii]] / binomial()$variance(g[,ii])
    } else {
      W <- (1 - h[,tx[ii]]) / binomial()$variance(g[,ii])
    }
    WX <- sweep(obj$model.matrix.gate, 1, W, `*`)
    hess[[ii]] <- chol2inv(crossprod(obj$model.matrix.gate, WX))
  }
  for (ii in seq_len(one)) {
    expert <- obj$experts[[ii]]
    par <- expert$parm
    lp <- length(par)
    sigma <- par[lp]
    V <- matrix(0, nrow=lp, ncol=lp)
    WX <- sweep(obj$model.matrix.expert, 1, expert$wlsfit$wt, `*`)
    V[seq_len(lp-1), seq_len(lp-1)] <- sigma * chol2inv(crossprod(obj$model.matrix.expert, WX))
    V[lp, lp] <- (2 * sigma^2) / (obj$N * sum(expert$wlsfit$wt))
    hess[[ii + ong]] <- V
  }
  return(hess)
}


cluster_par <- function(lstobj, nclus)
{
  LL <- kmeans(sapply(lstobj, logLik), nclus)
  goods <- lstobj[LL$cluster == which(LL$centers==max(LL$centers))]
  pars <- t(sapply(goods, function(x) x$par[unlist(x$par.expert.index)]))
  stack_par <- function(aa)
  {
    t(sapply(aa$par.expert.index, function(x) aa$par[x]))
  }
}


sample_se <- function(obj, alpha=0.1, no_replicates=99, burnin=50, 
                      type=c("bootstrap", "montecarlo"))
{
  resids <- sapply(obj$experts, function(x) x$wlsfit$residuals)
  H <- sweep(obj$weights$GD, 1, rowSums(obj$weights$GD), `/`)
  p <- max(sapply(obj$AR, function(x) max(x[[1]]))) # largest AR order
  resample <- function(r, wt)
  {
    nsamp <- nrow(r) + burnin + p
    M <- matrix(NA_real_, nrow=nsamp, ncol=ncol(r))
    for (cc in seq_len(ncol(r))) {
      M[,cc] <- sample(r[,cc], nsamp, replace=TRUE, prob=wt[,cc])
    }
    return(M)
  }
  # use posterior weights for re-sampling errors
  lstSamples <- replicate(no_replicates, resample(resids, H),
                          simplify=FALSE)
  # use prior weights for sampling state
  G <- if (is.null(obj$callNexp)) {
    do.call(cbind, obj$weights$G)
  } else {
    Reduce(`+`, obj$weights$G, init=0)
  }
  # Coefficients
  BETA <- lapply(obj$par.expert.index, function(x) obj$par[x])
  lstSIM <- vector("list", no_replicates)
  for (bb in seq_len(no_replicates)) {
    # sample states
    stateid <- c(sample.int(obj$num.experts, replace=TRUE, burnin, prob=G[1,]),
                 apply(G, 1, function(x) sample.int(obj$num.experts, 1, prob=x)))
    # forward filter
    simdat <- vector("numeric", length(stateid) + p)
    simdat[1:p] <- {ii <- sample.int(length(obj$Y)-p, 1); obj$Y[c(ii:(ii+p-1))]}
    conditional_mean <- function(state, tt)
    {
      nlag <- max(obj$AR[[state]][[1]])
      beta <- BETA[[state]]
      Xt <- simdat[(tt-1):(tt-nlag)]
      if (obj$AR[[state]][[2]])
        Xt <- c(1, Xt)
      as.numeric(crossprod(beta[-length(beta)], Xt))
    }
    for (jj in (p+1):length(simdat)) {
      ste <- stateid[jj-p]
      mu <- conditional_mean(ste, jj)
      simdat[jj] <- mu + lstSamples[[bb]][jj, ste]
    }
    simdat <- simdat[-(1:burnin)]
    simtst <- tsHME("y ~ y", DATA=data.frame(y=simdat),
      DEPTH=obj$depth, NEXPERTS=obj$num.experts, TOL=1e-3, AR=obj$AR,
      MAXEPOCHS=750, MAXIT=obj$maxit, PRIOR=obj$weights$prior)
    lstSIM[[bb]] <- list(parm=coef(simtst), loglik=logLik(simtst),
      converge=simtst$converge)
  }
  lstSIM
}
