setwd("~/Dropbox/HME/")
library(Formula)
source("hme_code/tsHME_helpers.R")

dtflm <- foreign::read.dta("data/20070061_data/Makowsky Stratmann - Pol Econ at Any Speed.dta")

#AR
#list
# 1 - p for AR(p)
form.cit <- "lnamount ~ -1 + outtown + outstate + orloss01 + lnmphover + lnpvaluepc + black + hispanic + female + lnage + femalelnage + statepol + otstatepol + osstatepol + splnpvaluepc + sporloss01 | -1 + cdl2 + outtown + outstate + orloss01 + lnmphover + lnpvaluepc + black + hispanic + female + lnage + femalelnage + statepol + otstatepol + osstatepol + splnpvaluepc + sporloss01"
form <- "lnamount ~ female | age + lnmphover + black"
form.dat <- "income ~ nexp | -1 + bigen + schl"

debugonce(lmHME)
hme.dat <- lmHME(form, dtflm, TRACE=1, NEXPERTS=3, DEPTH=2, MAXEPOCHS=75)

lmHME <- function(FORM, DATA, DEPTH=2, NEXPERTS=NULL, FIRST_SPLIT=2,
                  INIT=NULL, INIT_METHOD=c("random.sample"),
                  MAXEPOCHS=350, TOL=1e-3, PRIOR=NULL, MAXIT=5, TRACE=0,
                  TOL_FUNC=function(x) sqrt(crossprod(x)))
{
  
  cl <- match.call()
  stopifnot(!missing(DATA))
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("form", "data"), names(mf), 0)
  form <- Formula::as.Formula(FORM)
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = DATA)
  Y <- model.response(mf, "numeric")
  XE <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  XG <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  eparN <- ncol(XE) + 1 # Depends on dispersion
  gparN <- ncol(XG)
  nexp_call <- NEXPERTS
  if (is.null(NEXPERTS)) {
    enodeN <- NEXPERTS <- 2**DEPTH
    ss <- 2
  } else {
    enodeN <- ss <- NEXPERTS
  }
  epN <- NEXPERTS * eparN
  
  if (DEPTH==1) {
    gnode <- 1
  } else {
    gnode <- rep(2, DEPTH-1)
    gnode[1] <- FIRST_SPLIT
  }
  gnodeN <- cumsum(cumprod(gnode))[length(gnode)] + 1
  #!!!!!!
  gpN <- if (DEPTH==1) {
    ncol(XG) * (FIRST_SPLIT-1)
  } else {
    sum(c(1, cumprod(gnode)) * (c(gnode, ss) - 1) * ncol(XG))
  }
  #!!!!!
  totparN <- gpN  + epN
  if (is.null(INIT)) {
    INIT_METHOD <- match.arg(INIT_METHOD)
    if (INIT_METHOD == "random.sample") {
      smpidx <- sample.int(NEXPERTS, nrow(XE), replace=T)
      smpwts <- lapply(seq_len(NEXPERTS), function(ii) as.integer(smpidx==ii))
      smpmods <- mapply(init_lm, y=list(Y), x=list(XE), wt=smpwts, SIMPLIFY=FALSE)
    } else {
      stop("Only one INIT_METHOD available. Whah Whah.")
    }
    stpars <- lapply(smpmods, function(m) m$parm)
    PAR <- c(runif(gpN, min=-6, max=6), unname(unlist(stpars)))
  } else {
    PAR <- INIT
    stopifnot(length(PAR) == epN + gpN)
  }
  gpIdx <- gate_index(gnode, DEPTH, ss, ncol(XG))
  epIdx <- expert_index_lm(NEXPERTS, eparN, gpN)
  
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
    sapply(smpmods, get_density)
  } else {
    XX <- mapply(X_from_Y, y=list(Y), AR, SIMPLIFY=FALSE)
    pp <- mapply(function(x,p) p[x], epIdx, list(PAR), SIMPLIFY=FALSE)
    mapply(Gaussian, par=pp, y=list(Y), X=XX, ar=AR)
  }
  
  
  # Pre-Step
  colnames(matFs) <- as.character(seq_len(enodeN))
  estep <- e_step(matGs, matFs, gnode, nexp_call, NEXPERTS)
  logL[1,] <- sum(log(rowSums(estep$GD)))
  maxp <- 500
  
  
  for (ii in 2:MAXEPOCHS) {
    PAR_old <- PAR
    mstep <- m_step_lm(PAR, Y, XG, XE, estep, enodeN, gnode, ss, MAXIT)
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
    estep <- e_step(matGs, matFs, gnode, nexp_call, NEXPERTS)
    LL <- sum(log(rowSums(estep$GD)))
    parM[ii,] <- PAR
    logL[ii,] <- LL
    metric_norm <- TOL_FUNC(PAR - PAR_old)
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
  
  # Need to take find which parameter vector yields max like and re-run gates
  # and experts
  #p <- parM[which.max(logL), ]
  #mstep <- m_step_lm(PAR, Y, XG, XE, estep, enodeN, gnode, ss, MAXIT)
  #PAR <- unname(
  #  c(unlist(lapply(mstep$gate, `[[`, "parm")), 
  #    unlist(lapply(mstep$experts, `[[`, "parm")))
  #)
  #if (any(PAR < -maxp) || any(PAR > maxp)) {
  #  MAXIT <- 5
  #}
  #matFs <- sapply(mstep$experts, get_density)
  #colnames(matFs) <- as.character(seq_len(enodeN))
  #matGs <- multinomial_g(XG, PAR, gpIdx, gnode, ss, DEPTH, FIRST_SPLIT)
  #estep <- e_step(matGs, matFs, gnode, nexp_call, NEXPERTS)
  
  
  out <- list(par=p,                                          # parameter vector
              depth=DEPTH,                                    # depth of the network
              epochs=ii,                                      # no. of E-M cycles
              converge=as.integer(ii < MAXEPOCHS),            # did algo converge before maxepochs
              gates=mstep$gates,                              # glm output for gating networks
              experts=mstep$experts,                          # expert output
              weights=estep,                                  # list of various weights
              first_split=FIRST_SPLIT,                        # no. of splits from the root node
              loglik=logL[!is.na(logL)],                      # last log-likelihood from E-M process
              par.evolution=parM[rowSums(is.na(parM))==0,],   # matrix of parameter sequences
              par.gate.index=gpIdx,                           # index identifying gating parameters
              par.expert.index=epIdx,                         # index identifying expert parameter
              par.norm=metric_norm,                           # final value of metric_norm             
              model.matrix.gate=XG,                           # gating model matrix
              model.matrix.expert=XE,                         # expert model matrix
              Y=Y,                                            # vector of expert dependent variable
              N=length(Y),                                    # no. of observations
              gnode=gnode,                                    # splits per node for each level
              gates.per.level=cumprod(c(1, gnode)),           # no. of gates for each level
              level.node=cumprod(c(1, gnode, ss)),            # no. of gates and experts for each level
              num.experts=enodeN,                             # no. of unique experts
              num.expert.par=(ncol(XE) + 1),                  # no. of expert parameters. +1 depends on dispersion parameter(s)
              num.gate.par=gparN,                             # no. of unique gating parameters; not equal to total no. of parameters
              num.gates= if (DEPTH>1) sum(cumprod(c(1, gnode))) else 1, #!!!! total no. of gating nodes
              matFs=matFs,                                    # matrix where columns are expert density values
              prior.expert.mass=colSums(Reduce(`+`, estep$G)) / length(Y),  # percent of prior weights for each expert
              maxit=MAXIT,                                    # max no. of E-M cycles before stopping
              callNexp=nexp_call,                             # call to number of experts
              formula=FORM,                                   # original formula call
              call=cl)
  class(out) <- "lmHME"
  return(out)
}


summary.lmHME <- function(obj)
{
  theta <- obj$par
  ngp <- obj$num.gate.par
  nep <- obj$num.expert.par
  lstgatpar <- lapply(obj$par.gate.index, function(x) theta[x])
  lstexppar <- lapply(obj$par.expert.index, function(x) theta[x])
  sum_gate_par <- function(par, ng, rnames)
  {
    lstpar <- split(seq_along(par), ceiling(seq_along(par) / ng))
    for (i in seq_along(lstpar)) {
      p <- as.matrix(par[lstpar[[i]]])
      dimnames(p) <- list(rnames, "")
      print(p)
      cat("\n")
    }
  }
  cat(sprintf("Log-Likelihood: %f", log_Like(obj)))
  for (i in seq_len(obj$num.gates)) {
    p <- round(lstgatpar[[i]], 4)
    rnames <- get_coef_names_lm(obj$form, "gates")
    # rownames(p) <- paste(sprintf("w%d_", i), rnames, sep="_"); colnames(p) <- ""
    cat(sprintf("\n--------------------------\ngate-%d\n", i))
    sum_gate_par(p, ngp, rnames) # - only implement when we need more than two splits
  }
  cat("\n\n--------------------------\n--------------------------\n\n")
  for (j in seq_len(obj$num.experts)) {
    p <- as.matrix(round(lstexppar[[j]], 6))
    rnames <- get_coef_names_lm(obj$form, "experts")
    rownames(p) <- c(rnames, "var"); colnames(p) <- ""
    cat(sprintf("\n--------------------------\nexpert-%d\n", j))
    print(p)
  }
}
summary(hme.dat)



get_coef_names_lm <- function(form, part=c("experts", "gates"))
{
  part <- match.arg(part)
  patt <- if (part=="experts") {
    "(?<=~)(.*?)(?=\\|)"
  } else {
    "(?<=\\|).*"
  }
  eqtn <- regexpr(patt, form, perl=TRUE)
  eqtn <- trimws(strsplit(regmatches(form, eqtn), "\\+")[[1]])
  eqtn[-which(eqtn=="-1")]
}



debugonce(plot.lmHME)
plot(hme.dat, "par.conv.gate")
plot(hme.dat, "par.conv.expert")
