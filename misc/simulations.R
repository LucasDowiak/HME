#################
# Three regime sharp cut

sim_3r_ar <- function(TT, BT=round(0.10*TT))
{
  nar1 <- 1
  mod1_c <- 0.05
  mod1_ar1 <- 0.6
  var1 <- 0.03
  
  mod2_c <- 0.03
  mod2_ar1 <- -0.45
  mod2_ar2 <- 0
  var2 <- 0.01
  
  mod3_c <- 0
  mod3_ar1 <- 0.3
  mod3_ar2 <- -0.1
  mod3_ar3 <- 0
  var3 <- 0.06
  #breaks <- floor(seq(0, TT, by=TT/3)[-c(1,3+1)])
  #y <- rep(NA_real_, TT + BT)
  #y[1:nar1] <- rnorm(nar1, sd=sqrt(var1))
  # for (ii in 2:(BT + TT)) {
  #   if (ii > 0 && ii <= BT + breaks[1]) {
  #     # burn-in plus first ar
  #     y[ii] <- mod1_c + mod1_ar1 * y[ii-1] + rnorm(1, sd=sqrt(var1))
  #   } else if (ii > BT + breaks[1] && ii <= BT + breaks[2]) {
  #     # second ar
  #     y[ii] <- mod2_c + mod2_ar1 * y[ii-1] + mod2_ar2 * y[ii-2] + rnorm(1, sd=sqrt(var2))
  #     #} else if (ii > BT + breaks[2] && ii <= BT + breaks[3]) {
  #     # first ar again
  #     #  y[ii] <- mod1_c + mod1_ar1 * y[ii-1] + rnorm(1, sd=sqrt(var1))
  #   } else if (ii > BT + breaks[2]) {
  #     # third ar
  #     y[ii] <- mod3_c + mod3_ar1 * y[ii-1] + mod3_ar2 * y[ii-2] + 
  #       mod3_ar3 * y[ii-3] + rnorm(1, sd=sqrt(var3))
  #   }
  # }
  # y[-(1:BT)]
  st <- ifelse(TT%%3==0, 1, 0)
}




#################
# Two Regime Smooth

sim_2r_ar_smooth <- function(TT, xrange=c(-10, 10), phi=1, PLOT=TRUE,
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



#################
# Two Regime Sharp

sim_2r_ar_smooth <- function(TT, BT=round(0.10*TT))
{
  nar1 <- 1
  mod1_c <- 0.05
  mod1_ar1 <- 0.6
  var1 <- 0.03
  
  mod2_c <- 0.03
  mod2_ar1 <- -0.45
  mod2_ar2 <- 0
  var2 <- 0.01
  
  G <- matrix(0, TT, 2)
  b1 <- round(TT/3)
  b2 <- 2 * b1
  G[1:b1, 1] <- 1
  G[b1:b2, 1] <- (seq(b2,b2-b1) / (b2-b1)) - 1
  G[(b2+1):TT, 1] <- 0
  G[,2] <- 1 - G[,1]
  
  randstate <- c(sample.int(ncol(G), replace=TRUE, BT, prob=G[1,]),
                 apply(G, 1, function(x) sample.int(ncol(G), 1, prob=x)))
  
  y <- rep(NA_real_, TT + BT)
  y[1:nar1] <- rnorm(nar1, sd=sqrt(var1))
  for (tt in 2:length(y)) {
    if (randstate[tt]==1) {
      y[tt] <- mod1_c + mod1_ar1 * y[tt-1] + rnorm(1, sd=sqrt(var1))
    } else {
      y[tt] <- mod2_c + mod2_ar1 * y[tt-1] + rnorm(1, sd=sqrt(var2))
    }
  }
  y[-(1:BT)]
}


yy_whitenoise <- rnorm(1000)


unembed <- function(M)
{
  stopifnot(is.matrix(M))
  matrix(c(M[1,-1], M[,1]))
}


respecify <- function(obj, maxar=7, constant=TRUE, criterion=c("AIC", "BIC"),
                      wtscheme=c("softmax", "select_one"))
{
  criterion <- match.arg(criterion)
  wtscheme <- match.arg(wtscheme)
  if (criterion=="AIC")
    crit <- AIC
  else
    crit <- BIC
  nr <- nrow(obj$Y)
  maxlag <- max(sapply(obj$AR, function(x) max(x[[1]])))
  Y <- as.data.frame(embed(unembed(obj$Y), maxar + 1))
  names(Y) <- c("y", paste0("yL", 1:maxar))
  WTS <- sweep(obj$weights$GD, 1, rowSums(obj$weights$GD), `/`)[-(1:(maxar-maxlag)),]
  if (wtscheme=="select_one") {
    maxreg <- apply(WTS, 1, which.max)
    WTS1 <- matrix(0, nrow=length(WTS))
    WTS1[(seq_len(nrow(WTS)) - 1) * ncol(WTS) + maxreg] <- 1
    WTS <- matrix(c(WTS1), nrow=nrow(WTS), byrow=T)
  }
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
  for (ee in seq_along(out)) {
    cat(sprintf("Expert %d \n", ee))
    print(out[[ee]])
    cat(sprintf("minimum = %.3f \n\n", min(out[[ee]])))
  }
  invisible(out)
}


y <- sim_3r_ar(300)

tst1 <- tsHME("y ~ .",
              DATA=data.frame(y=y),
              DEPTH=2,
              NEXPERTS=3,
              TOL=1e-1,
              MAXEPOCHS=50,
              AR=list(list(1,FALSE,FALSE,FALSE),
                      list(1,FALSE,FALSE,FALSE),
                      list(1,FALSE,FALSE,FALSE)),
              MAXIT=5,
              FIRST_SPLIT=3, INIT=tst1$par) #

respecify(tst1, criterion="BIC", wtscheme="select_one")

tst2 <- tsHME("y ~ .",
              DATA=data.frame(y=y),
              DEPTH=3,
              NEXPERTS=3,
              TOL=1e-1,
              MAXEPOCHS=250,
              AR=list(list(1:2,FALSE,FALSE,FALSE),
                      list(1,TRUE,FALSE,FALSE),
                      list(1:2,FALSE,FALSE,FALSE)),
              MAXIT=10,
              FIRST_SPLIT=3) # , PRIOR=tst1$weights$prior)
              #INIT=tst2$par)

respecify(tst2, criterion="BIC")


lstBootSmp <- sample_se(tst2, no_replicates = 1100)

plot(density(sapply(lstSmp, `[[`, "loglik")))

parms <- rbind(tst2$par,
               t(sapply(lstSmp, `[[`, "parm")))
epidx <- unlist(tst2$par.expert.index)
alpha <- 0.025



glmnet(tst1$model.matrix.gate, tst1$weights$posterior$`0.2`, family="binomial",
       weights=tst1$weights$posterior$`0`[,2], intercept=FALSE)



for (pp in cppckgs) {
  from <- paste0(pth, "3.2/Resources/library" pp, sep="/")
  print(from)
}


