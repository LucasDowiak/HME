setwd("~/Dropbox/HME/simulations/")
source("../hme_code/tsHME.R")
source("../hme_code/tsHME_helpers.R")



sim_2r_ar_sharp <- function(TT, xrange=c(-10, 10), phi=1, PLOT=TRUE,
                            ar1=list(c(0.6, 0.2), c(0.8,0.1)))
{
  x <- seq(xrange[1], xrange[2], length.out=TT)
  y <- (tanh(phi * x) + 1) / 2
  if (PLOT) {
    plot(1:TT, y, type="l")
    grid()
  }
  G <- matrix(c(1-y, y), ncol=2)
  randstate <- c(1, apply(G, 1, function(x) sample.int(ncol(G), 1, prob=x)))
  z <- vector("numeric", TT + 1)
  z[1] <- arima.sim(1, model=list(ar=ar1[[1]][1]), sd=sqrt(ar1[[1]][2]))
  for (tt in 2:length(z)) {
    if (randstate[tt]==1) {
      z[tt] <- ar1[[1]][1] * z[tt-1] + rnorm(1, sd=sqrt(ar1[[1]][2]))
    } else {
      z[tt] <- ar1[[2]][1] * z[tt-1] + rnorm(1, sd=sqrt(ar1[[2]][2]))
    }
  }
  return(z[-1])
}


sim_3r_ar_smooth <- function(TT, xrange=c(-10, 10), phi1=1.5, phi2=1.5, PLOT=TRUE)
{
  
  mod1_ar1 <- 0.7
  mod1_ar2 <- 0
  mod1_ar3 <- 0
  var1 <- 0.03
  
  mod2_ar1 <- -0.5
  mod2_ar2 <- 0.3
  mod2_ar3 <- 0
  var2 <- 0.05
  
  mod3_ar1 <- 0.05
  mod3_ar2 <- 0.1
  mod3_ar3 <- 0.6
  var3 <- 0.1
  
  x <- seq(xrange[1], xrange[2], length.out=TT)
  y <- (tanh(phi1 * (x - (-20/5))) + 1) / 2
  yy <- (tanh(phi2 * (x - (20/5))) + 1) / 2
  y1 <- 1-y
  y2 <- 1 - (1-y+yy)
  y3 <- 1-(y1+y2)
  
  G <- matrix(c(y1, y2, y3), ncol=3)
  randstate <- c(1, apply(G, 1, function(x) sample.int(ncol(G), 1, prob=x)))
  z <- vector("numeric", TT + 1)
  sigs <- list(var1, var2, var3)
  errors <- sapply(sigs, function(x) rnorm(TT, sd=sqrt(x)))
  for (tt in 2:length(z)) {
    if (randstate[tt]==1) {
      z[tt] <- mod1_ar1 * z[tt-1] + errors[tt-1, 1]
    } else if (randstate[tt]==2) {
      z[tt] <- mod2_ar1 * z[tt-1] + mod2_ar2 * z[tt-2] + mod2_ar3 * z[tt-3] + errors[tt-1, 2]
    } else if (randstate[tt]==3) {
      z[tt] <- mod3_ar1 * z[tt-1] + mod3_ar2 * z[tt-2] + mod3_ar3 * z[tt-3] + errors[tt-1, 3]
    }
  }
  Z <- z[-1]
  if (PLOT) {
    par(mfrow=c(2,1))
    plot(1:TT, y1, type="l", col="blue")
    lines(y2, col="red")
    lines(y3, col="orange"); grid()
    plot(Z, type="l"); grid()
  }
  S <- matrix(0, nrow=nrow(G), ncol=ncol(G))
  for (rr in seq_len(TT)) {
    S[rr, randstate[rr + 1]] <- 1
  }
  LL <- mapply(function(x) dnorm(errors[,x], sd=sqrt(sigs[[x]]), log=T),
               seq_len(3))
  attr(Z, "log-lik") <- sum(LL * S)
  return(Z)
}


respecify <- function(obj, maxar=7, constant=TRUE, criterion=c("AIC", "BIC"),
                      display=FALSE)
{
  unembed <- function(M)
  {
    stopifnot(is.matrix(M))
    matrix(c(M[1,-1], M[,1]))
  }
  criterion <- match.arg(criterion)
  if (criterion=="AIC")
    crit <- AIC
  else
    crit <- BIC
  nr <- nrow(obj$Y)
  maxlag <- max(sapply(obj$AR, function(x) max(x[[1]])))
  Y <- as.data.frame(embed(unembed(obj$Y), maxar + 1))
  names(Y) <- c("y", paste0("yL", 1:maxar))
  WTS <- sweep(obj$weights$GD, 1, rowSums(obj$weights$GD), `/`)[-(1:(maxar-maxlag)),]
  out <- vector("list", obj$num.experts)
  for (ee in 1:obj$num.experts) {
    Y$wts <- WTS[,ee]
    critvalues <- matrix(NA_real_, nrow=maxar, ncol=1+constant)
    for (ll in 1:maxar) {
      form <- paste("y ~ 0", paste(names(Y)[1:ll + 1], collapse=" + "), sep=" + ")
      lmmod <- lm(form, data=Y, weights=wts)
      critvalues[ll, 1] <- crit(lmmod)
      if (constant) {
        form <- paste("y ~", paste(names(Y)[1:ll + 1], collapse=" + "))
        lmmod <- lm(form, data=Y, weights=wts)
        critvalues[ll, 2] <- crit(lmmod)
      }
      rownames(critvalues) <- paste0("ar", 1:maxar)
      colnames(critvalues) <- c("No Constant", if (constant) "Constant")
    }
    out[[ee]] <- critvalues
  }
  if (display) {
    for (ee in seq_along(out)) {
      cat(sprintf("Expert %d \n", ee))
      print(out[[ee]])
      cat(sprintf("minimum = %.3f \n\n", min(out[[ee]])))
    }
  }
  invisible(out)
}


concat_respec <- function(idx, XX)
{
  which(XX[[idx]]==min(XX[[idx]]), arr.ind=TRUE)
}

################################################################################
# inputs
DD <- 2
first_split <- 2
nexperts <- 3
ARspec <- list(list(1,FALSE,FALSE,FALSE),
               list(1,FALSE,FALSE,FALSE),
               list(1,FALSE,FALSE,FALSE))
sample.size <- seq(100, 1000, by=100)
NN <- 100
ffile <-"respecify_R3D2S2.RDS"


################################################################################


lst.d <- vector("list", length(sample.size) * NN)
names(lst.d) <- rep(sample.size, each=NN)
attr(lst.d, "beta") <-   c(m1ar1=0.7,  m1s=0.03,
                           m2ar1=-0.5, m2ar2=0.3, m2s=0.05,
                           m3ar1=0.05,  m3ar2=0.1, m3ar3=0.6,  m3s=0.1,
                           NN=100, from=100, to=1000, by=100)
  
ii <- 1
for (ss in sample.size) { 
  jj <- 0
  while (jj < NN) {
    ysim <- sim_3r_ar_smooth(ss, PLOT=FALSE)
    mod.d <- try(multi_call_search("y ~ .", DATA=data.frame(y=as.vector(ysim)),
                                   DEPTH=DD, NEXPERTS=nexperts, FIRST_SPLIT=first_split,
                                   AR=ARspec, TOL=1e-3, maxloops=6, maxepochs=75))
    if (inherits(mod.d, "try-error"))
      next

    respecAIC <- respecify(mod.d, constant=FALSE, criterion="AIC", display=FALSE)
    respecBIC <- respecify(mod.d, constant=FALSE, criterion="BIC", display=FALSE)
    G <- Reduce(`+`, mod.d$weights$G, init=0)
    s1 <- which.max(colMeans(head(G, nrow(G)/3)))
    s3 <- which.max(colMeans(tail(G, nrow(G)/3)))
    s2 <- (1:3)[!(1:3) %in% c(s1,s3)]
    lst.d[[ii]] <- list(AIC=list(do.call(c, lapply(c(s1,s2,s3), concat_respec, XX=respecAIC))),
                        BIC=list(do.call(c, lapply(c(s1,s2,s3), concat_respec, XX=respecBIC))))

    jj <- jj + 1
    ii <- ii + 1
  }
  print(sprintf("Sample size %d finished at %s", ss, Sys.time()))
  saveRDS(lst.d, file=ffile)
}




#####################################################################################
# Plotting
library(plotly)

lstd <- readRDS("respecify_R3D2S2.RDS")

NN <- unique(names(lstd))

cbind(attr(lstd, "beta"), names(attr(lstd, "beta")))

plot_modsel_hist_3R <- function(lst, n)
{
  lst <- lst[names(lst) %in% n]
  AIC <- as.data.frame(do.call(rbind, sapply(lst, `[[`, "AIC")), row.names=F); AIC$crit <- "AIC"
  BIC <- as.data.frame(do.call(rbind, sapply(lst, `[[`, "BIC")), row.names=F); BIC$crit <- "BIC"
  CRIT <- melt(rbind(AIC, BIC), id.vars="crit")
  CRIT$crit <- as.factor(CRIT$crit)
  require(plotly)
  plot_hist <- function()
  {
    p1 <- plot_ly(alpha=0.6, colors=c("orange", "sky blue")) %>% 
            add_histogram(x=AIC[,1], name="AIC") %>%
            add_histogram(x=BIC[,1], name="BIC") %>%
            layout(barmode = "overlay", title="AR(1)")
    p2 <- plot_ly(alpha=0.6, colors=c("orange", "sky blue")) %>% 
            add_histogram(x=AIC[,3], name="AIC", showlegend=F) %>%
            add_histogram(x=BIC[,3], name="BIC", showlegend=F) %>%
            layout(barmode = "overlay")
    p3 <- plot_ly(alpha=0.6, colors=c("orange", "sky blue")) %>% 
            add_histogram(x=AIC[,5], name="AIC", showlegend=F) %>%
            add_histogram(x=BIC[,5], name="BIC", showlegend=F) %>%
            layout(barmode = "overlay")
    subplot(p1, p2, p3, shareY=TRUE)
  }
  plot_hist()
}

p <- plot_ly(subset(CRIT, variable=="V1"), x = ~value, color=~crit, alpha=0.6, type="histogram") 



percent_correct <- function(n, lst, good)
{
  lst <- lst[names(lst) %in% n]
  AIC <- do.call(rbind ,lapply(lst, function(x) x$AIC[[1]][1:6]))
  BIC <- do.call(rbind ,lapply(lst, function(x) x$BIC[[1]][1:6]))
  out <- rbind(colSums(sweep(AIC[,c(1,3,5)], 2, good, `==`)) / nrow(AIC),
               colSums(sweep(BIC[,c(1,3,5)], 2, good, `==`)) / nrow(BIC))
  dimnames(out) <- list(c("AIC", "BIC"), c("R1", "R2", "R3"))
  out
}


simulate_and_test <- function(N, ar=list(ar=c(0.8)), sd, maxar=7, Nsim=100,
                              CRIT=c("AIC", "BIC"))
{
  library(forecast)
  CRIT <- match.arg(CRIT)
  if (CRIT=="AIC")
    crit <- AIC
  else
    crit <- BIC
  arlen <- length(ar[["ar"]])
  Y <- arima.sim(n = N + arlen, ar, sd=sd)
  Y <- as.data.frame(embed(as.matrix(Y), 1 + maxar))
  names(Y) <- c("y", paste0("yL", seq_len(maxar)))
  out <- vector("integer", Nsim)
  for (jj in seq_len(Nsim)) {
    Y <- arima.sim(n = N + arlen, ar, sd=sd)
    Y <- as.data.frame(embed(as.matrix(Y), 1 + maxar))
    names(Y) <- c("y", paste0("yL", seq_len(maxar)))
    critvalues <- matrix(NA_real_, nrow=maxar, ncol=1)
    for (ll in 1:maxar) {
      form <- paste("y ~ 0", paste(names(Y)[1:ll + 1], collapse=" + "), sep=" + ")
      lmmod <- lm(form, data=Y)
      critvalues[ll, 1] <- crit(lmmod)
    }
    rownames(critvalues) <- paste0("ar", 1:maxar)
    colnames(critvalues) <- "No Constant"
    out[jj] <- as.integer(which.min(critvalues))
  }
  sum(out==arlen)/length(out)
}

lapply(NN, simulate_and_test, ar=list(ar=c(0.05, 0.1, 0.6)), sd=0.1, maxar=7, Nsim=100,
       CRIT="BIC")



