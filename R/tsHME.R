#setwd("~/Dropbox/HME/hme_code/")
library(data.table)
#source("tsHME_helpers.R")


#AR
#list
# 1 - p for AR(p)

tsHME <- function(FORM, DATA, DEPTH=2, NEXPERTS=NULL, FIRST_SPLIT=2,
                  AR=list(list(1,TRUE,FALSE,FALSE), list(1,FALSE,FALSE,FALSE)),
                  INIT=NULL, INIT_METHOD=c("time.spec","model.spec"),
                  MAXEPOCHS=350, TOL=1e-3, PRIOR=NULL, MAXIT=5, TRACE=0,
                  TOL_FUNC=function(x) sqrt(crossprod(x)))
{
  form <- as.formula(FORM)
  if (DEPTH==1) {
    gnode <- 1
  } else {
    gnode <- rep(2, DEPTH-1)
    gnode[1] <- FIRST_SPLIT
  }
  gnodeN <- cumsum(cumprod(gnode))[length(gnode)] + 1
  
  if (is.null(NEXPERTS)) {
    enodeN <- NEXPERTS <- 2**DEPTH
    ss <- 2
  } else {
    enodeN <- ss <- NEXPERTS
  }
  stopifnot(length(AR) == NEXPERTS)
  mf <- model.frame(form, data=DATA)
  maxar <- max(sapply(AR, function(x) max(x[[1]]))) + any(sapply(AR, `[[`, 4))
  if (any(sapply(AR, `[[`, 4)))
    maxar <- max(1, maxar)
  Y <- embed(model.response(mf, "numeric"), 1 + maxar)
  # Extra parameters depend on assumed parametric density
  # In this case, it's the dispersion parameter
  epN <- sum(sapply(AR, get_n_par))
  N <- nrow(Y)
  if (grepl("~\\.", gsub(" ", "", FORM))) {
    XG <- matrix(c(rep(1,N), T=(1:N)/N), ncol=2)
    colnames(XG) <- c("(INTERCEPT)", "T")
    gparN <- 2
  } else {
    XG <- model.matrix(terms(form, data=mf), mf)
    XG <- XG[-(1:maxar),]
    gparN <- ncol(XG)
  }
  #!!!!!!
  gpN <- if (DEPTH==1) {
    ncol(XG) * (FIRST_SPLIT-1)
  } else {
    sum(c(1, cumprod(gnode)) * (c(gnode, ss) - 1) * ncol(XG))
  }
  #!!!!!
  if (is.null(INIT)) {
    INIT_METHOD <- match.arg(INIT_METHOD)
    if (INIT_METHOD == "model.spec") {
      arcoef <- lapply(AR, function(ar) cll_ar(Y, ar))
    } else if (INIT_METHOD == "time.spec") {
      wts <- starting_weights(Y, AR)
      arcoef <- mapply(cll_ar, y=list(Y), ar=AR, wt=wts, SIMPLIFY=FALSE)
    }
    PAR <- c(runif(gpN, min=-10, max=10), unname(unlist(lapply(arcoef, `[[`, "parm"))))
  } else {
    PAR <- INIT
    stopifnot(length(PAR) == epN + gpN)
  }
  totparN <- gpN  + epN
  gpIdx <- gate_index(gnode, DEPTH, ss, ncol(XG))
  epIdx <- expert_index(AR, gpN)
  
  # Create matrix to track evolution of parameter values
  logL <- matrix(NA_real_, nrow=MAXEPOCHS+1)
  parM <- matrix(NA_real_, ncol=length(PAR), nrow=MAXEPOCHS+1)
  parM[1,] <- PAR
  
  # gates
  matGs <- if (is.null(PRIOR)) {
    multinomial_g(XG, PAR, gpIdx, gnode, ss, DEPTH, FIRST_SPLIT)
  } else {
    nr <- nrow(PRIOR[[1]])
    if (nr > N) {
      for (ll in seq_along(PRIOR)) {
        PRIOR[[ll]] <- PRIOR[[ll]][-seq_len(nr-N),]
      }
    } else if (nr < N) {
      for (ll in seq_along(PRIOR)) {
        rr <- PRIOR[[ll]][1,]
        PRIOR[[ll]] <- rbind(rbind(rr, N-nr), PRIOR[[ll]])
      }
    }
    PRIOR
  }
  # experts
  matFs <- if (is.null(INIT)) {
    sapply(arcoef, get_density)
  } else {
    XX <- mapply(X_from_Y, y=list(Y), AR, SIMPLIFY=FALSE)
    pp <- mapply(function(x,p) p[x], epIdx, list(PAR), SIMPLIFY=FALSE)
    mapply(Gaussian, par=pp, y=list(Y), X=XX, ar=AR)
  }
  colnames(matFs) <- as.character(seq_len(enodeN))
  estep <- e_step(matGs, matFs, gnode)
  logL[1,] <- sum(log(rowSums(estep$GD)))
  maxp <- 500
  for (ii in 2:MAXEPOCHS) {
    PAR_old <- PAR
    mstep <- m_step(PAR, Y, XG, AR, estep, enodeN, gnode, NEXPERTS, MAXIT)
    PAR <- unname(
      c(unlist(lapply(mstep$gate, `[[`, "parm")), 
        unlist(lapply(mstep$experts, `[[`, "parm")))
    )
    if (any(PAR < -maxp) || any(PAR > maxp)) {
      MAXIT <- 5
    }
    matFs <- sapply(mstep$experts, get_density)
    colnames(matFs) <- as.character(seq_len(enodeN))
    matGs <- multinomial_g(XG, PAR, gpIdx, gnode, ss, DEPTH, FIRST_SPLIT)
    estep <- e_step(matGs, matFs, gnode)
    LL <- sum(log(rowSums(estep$GD)))
    parM[ii,] <- PAR
    logL[ii,] <- LL
    metric_norm <- TOL_FUNC(PAR - PAR_old)
    cat('\r', sprintf("Square Norm: %f, Log-Like: %f, Epoch: %d", metric_norm,
                      LL, ii))
    if (TRACE > 0) {
      cat('\r', sprintf("Square Norm: %f, Log-Like: %f, Epoch: %d", metric_norm,
                        LL, ii))
    } else {
      cat('\n', sprintf("Square Norm: %f, Log-Like: %f, Epoch: %d", metric_norm,
                        LL, ii))
    }
    if (metric_norm < TOL)
      break
  }
  out <- list(par=PAR,
              depth=DEPTH,
              epochs=ii,
              converge=as.integer(ii < MAXEPOCHS),
              gates=mstep$gates,
              experts=mstep$experts,
              weights=estep,
              first_split=FIRST_SPLIT,
              loglik=logL[!is.na(logL)],
              par.evolution=parM[rowSums(is.na(parM))==0,],
              par.gate.index=gpIdx,
              par.expert.index=epIdx,
              par.norm=metric_norm,
              model.matrix.gate=XG,
              AR=AR, #model.matrix.expert=XE,
              Y=Y,
              N=length(Y),
              gnode=gnode,
              gates.per.level=cumprod(c(1, gnode)),
              level.node=cumprod(c(gnode, ss)),
              num.experts=enodeN,
              num.expert.par=length(AR),
              num.gate.par=gparN,
              num.gates= if (DEPTH>1) sum(cumprod(c(1, gnode))) else 1, #!!!!
              matFs=matFs,
              maxit=MAXIT,
              callNexp=NEXPERTS,
              formula=as.character(form))
  class(out) <- "tsHME"
  return(out)
}



multi_call_search <- function(..., startmaxit=3, maxloops=5, maxepochs=50,
                              maxspeed=25, minspeed=3)
{
  go_get <- function(lst, inp_lst) {
    jj <- which.max(sapply(lst, function(x) max(x$loglik)))
    ii <- max(which.max(lst[[jj]]$loglik) - 1, 1)
    inp_lst$INIT <- lst[[jj]]$par.evolution[ii,]
    inp_lst$MAXEPOCHS <- 2
    inp_lst$MAXIT <- lst[[jj]]$maxit
    do.call(tsHME, inp_lst)
  }
  
  input_list <- as.list(substitute(list(...)))[-1L]
  input_list$MAXIT <- startmaxit
  input_list$MAXEPOCHS <- 100
  tmp <- replicate(10, do.call(tsHME, input_list), simplify=FALSE)
  tmp1 <- go_get(tmp, input_list)
  
  input_list$MAXIT <- startmaxit
  input_list$MAXEPOCHS <- maxepochs
  input_list$INIT <- tmp1$par
  mod_hist <- vector("list", maxloops)
  ll <- 1
  while (ll <= maxloops) {
    tmp2 <- do.call(tsHME, input_list)
    mod_hist[[ll]] <- tmp2
    ii <- which.max(tmp2$loglik)
    input_list$INIT <- tmp2$par.evolution[ii,]
    if (ii == input_list$MAXEPOCHS) {
      input_list$MAXIT <- min(input_list$MAXIT + 3, maxspeed)
    } else {
      input_list$MAXIT <- max(input_list$MAXIT - 5, minspeed)
    }
    ll <- ll + 1
  }
  mod_hist <- Filter(Negate(is.null), mod_hist)
  out <- go_get(mod_hist, input_list)
  return(out)
}


# Plot gold as random walk vs trendstationary
tmp2 <- multi_call_search("AU ~ .",TRACE=1,
                          DATA=dtf2[year(Date) >= 2007 & year(Date) <= 2012],
                          DEPTH=2, NEXPERTS=2, TOL=1e-2, FIRST_SPLIT=3, MAXEPOCHS = 300,
                          AR=list(list(1,TRUE,TRUE,FALSE),
                                  list(1,TRUE,FALSE,TRUE)),
                          startmaxit=5, maxloops=5, maxepochs=100)


obj <- tmp1
predict.tsHME <- function(obj, newdata=dtf)
{
  stopifnot(inherits(dtf, "data.frame"))
}


