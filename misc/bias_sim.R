setwd("~/Dropbox/HME/simulations/")
source("../hme_code/tsHME.R")
source("../hme_code/tsHME_helpers.R")



sim_2r_ar_smooth <- function(TT, xrange=c(-10, 10), phi=1, PLOT=TRUE,
                            ar1=list(c(0.6, 0.2), c(0.8,0.1)))
{
  x <- seq(xrange[1], xrange[2], length.out=TT)
  y <- (tanh(phi * x) + 1) / 2
  if (PLOT) {
    plot(1:TT, 1-y, type="l", yaxt="n", xaxt="n", ylab="P(s=i)", xlab="", bty="l",
         main="Two-Expert Membership Probability")
    axis(1, at=TT, labels="T"); axis(2, at=1, las=2, labels="1")
    grid()
    legend(x=5, y=0.5, c("Expert 1"), lty=1, col="black")
  }
  G <- matrix(c(1-y, y), ncol=2)
  randstate <- c(1, apply(G, 1, function(x) sample.int(ncol(G), 1, prob=x)))
  z <- vector("numeric", TT + 1)
  sigs <- lapply(ar1, function(x) x[length(x)])
  errors <- sapply(sigs, function(x) rnorm(TT, sd=sqrt(x)))
  for (tt in 2:length(z)) {
    if (randstate[tt]==1) {
      z[tt] <- ar1[[1]][1] * z[tt-1] + errors[tt-1, 1] 
    } else {
      z[tt] <- ar1[[2]][1] * z[tt-1] + errors[tt-1, 2] 
    }
  }
  Z <- z[-1]
  S <- matrix(0, nrow=nrow(G), ncol=ncol(G))
  for (rr in seq_len(TT)) {
    S[rr, randstate[rr + 1]] <- 1
  }
  LL <- mapply(function(x) dnorm(errors[,x], sd=sqrt(sigs[[x]]), log=T),
               seq_len(length(ar1)))
  attr(Z, "log-lik") <- sum(LL * S)
  return(list(y=Z, G=G))
}


sim_3r_ar_smooth <- function(TT, xrange=c(-10, 10), phi1=1.5, phi2=1.5, PLOT=TRUE,
                             ar1=list(c(0.6, 0.2), c(-0.4,0.15), c(0.8,0.01)))
{
  x <- seq(xrange[1], xrange[2], length.out=TT)
  y <- (tanh(phi1 * (x - (-20/5))) + 1) / 2
  yy <- (tanh(phi2 * (x - (20/5))) + 1) / 2
  y1 <- 1-y
  y2 <- 1 - (1-y+yy)
  y3 <- 1-(y1+y2)
  if (PLOT) {
    plot(1:TT, y1, type="l", col="blue", yaxt="n", xaxt="n", ylab="P(s=i)", xlab="", bty="l",
         main="Three-Expert Membership Probability")
    axis(1, at=TT, labels="T"); axis(2, at=1, las=2, labels="1")
    lines(y2, col="red")
    lines(y3, col="orange")
    grid()
    legend(x=375, y=0.75, c("Expert 1", "Expert 2", "Expert 3"),
           lty=1, col=c("blue", "red", "orange"))
  }
  G <- matrix(c(y1, y2, y3), ncol=3)
  randstate <- c(1, apply(G, 1, function(x) sample.int(ncol(G), 1, prob=x)))
  z <- vector("numeric", TT + 1)
  sigs <- lapply(ar1, function(x) x[length(x)])
  errors <- sapply(sigs, function(x) rnorm(TT, sd=sqrt(x)))
  for (tt in 2:length(z)) {
    if (randstate[tt]==1) {
      z[tt] <- ar1[[1]][1] * z[tt-1] + errors[tt-1, 1] # rnorm(1, sd=sqrt(ar1[[1]][2]))
    } else if (randstate[tt]==2) {
      z[tt] <- ar1[[2]][1] * z[tt-1] + errors[tt-1, 2] # rnorm(1, sd=sqrt(ar1[[2]][2]))
    } else if (randstate[tt]==3) {
      z[tt] <- ar1[[3]][1] * z[tt-1] + errors[tt-1, 3] # rnorm(1, sd=sqrt(ar1[[3]][2]))
    }
  }
  Z <- z[-1]
  S <- matrix(0, nrow=nrow(G), ncol=ncol(G))
  for (rr in seq_len(TT)) {
    S[rr, randstate[rr + 1]] <- 1
  }
  LL <- mapply(function(x) dnorm(errors[,x], sd=sqrt(sigs[[x]]), log=T),
               seq_len(length(ar1)))
  attr(Z, "log-lik") <- sum(LL * S)
  return(Z)
}



calc_bias <- function(lst, BETA, PLOT=TRUE)
{
  stopifnot(!missing(BETA))
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
  find_bias <- function(P)
  {
    cls <- kmeans(P[,-1], 2)
    h <- which.max(cls$centers[,1])
    l <- which.min(cls$centers[,1])
    tbl <- table(cls$cluster)
    if (any(tbl!=tbl[1]))
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
           ylab="bias", xlab="series length", cex.main=2)
      lines(ss, predict(loess("y~t", data.frame(y=b, t=ss), span=1)), col="green")
      #lines(ss, predict(smooth.spline(x=b,y=ss,spar=0.5)))
      abline(h=0, col="red"); grid()
    }
  }
  return(matBias)
}


#AR1 <- list(c(0.8, 0.15), c(-0.7, 0.1) c(0.5,0.07))
################################################################################
DD <- 2
first_split <- 2
nexperts <- 2
sample.size <- round(seq(100, 2000, by=(100/3)))
NN <- 250
sim_type <- sim_2r_ar_smooth
ARsim <- list(c(0.8, 0.1), c(0.5, 0.03))
ARspec <- list(list(1,FALSE,FALSE,FALSE),
               list(1,FALSE,FALSE,FALSE))
ffile <-"bias_R2D2S2_1300_1500.RDS" #  "bias_R3D2S2.RDS"
################################################################################
# Bias Simulation

lst.d <- vector("list", length(sample.size) * NN)
attr(lst.d, "beta") <- unlist(ARsim)
ii <- 1
for (ss in sample.size) { 
  jj <- 0
  while (jj < NN) {
    ysim <- sim_type(100, ar1=ARsim, PLOT=FALSE)
    mod.d <- try(tsHME("y ~ .", DATA=data.frame(y=as.vector(ysim$y)),
                        DEPTH=DD, NEXPERTS=nexperts, FIRST_SPLIT=first_split, PRIOR=ysim['G'],
                        AR=ARspec, MAXIT=10, MAXEPOCHS=200, TOL=1e-3)) #, maxloops=6))
    if (inherits(mod.d, "try-error"))
      next
    pp <- mod.d$par[unlist(mod.d$par.expert.index)]
    if (any(pp < 0) || any(pp > 0.99))
      next
    lst.d[[ii]] <- list(ss=ss, d=mod.d$depth, par=mod.d$par, idx=mod.d$par.expert.index,
                        LL=mod.d$loglik[length(mod.d$loglik)])
    jj <- jj + 1
    ii <- ii + 1
  }
  print(sprintf("Sample size %d finished at %s", ss, Sys.time()))
  saveRDS(lst.d, file=ffile)
}

if (FALSE) {
lst1 <- readRDS("bias_R2D2S2_100_550.RDS")
lst2 <- readRDS("bias_R2D2S2_575_1025.RDS")
lst3 <- readRDS("bias_R2D2S2_1050_1500.RDS")
lst <- c(lst1, lst2, Filter(function(x) length(x) > 0, lst3))

lst <- readRDS("bias_R2D1.RDS")
debugonce(calc_bias)
calc_bias(lst, c(0.8, 0.1, 0.5, 0.03))

}