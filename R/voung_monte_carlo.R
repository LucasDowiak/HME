setwd("~/Git/hme/")
source("R/benchmark_hme.R")


# --------------- Monte Carlo Building Blocks --------------- #

logit <- function(x, scale=1, mid=0)
{
  1 / (1 + exp(scale * (mid - x)))
}


get_means <- function(b) b["c"] + b["z"] * 0.5 


order_betas <- function(b)
{
  const <- sapply(b, `[[`, 1)
  return(order(const))
}


generate_obs <- function(n, rho=1)
{
  M <- matrix(nrow=n, ncol=2, dimnames=list(NULL, c("x", "z")))
  for (ii in seq_len(n)) {
    X <- rnorm(2)
    Z <- (X[1] + rho * X[2]) / sqrt(1 + rho**2)
    M[ii, ] <- c(X[1], pnorm(Z))
  }
  M
}


sample_experts <- function(X, scale=1, mid=0.5)
{
  # Should we switch this be dependent on the X variable?
  stopifnot(scale > 0)
  W <- X
  W[, 1] <-  logit(X[, 2], scale=scale, mid=mid)
  W[, 2] <-  1 - W[, 1]
  apply(W, 1, function(x) sample.int(2, 1, prob=x))
}


sample_multinomial <- function(X, Om)
{
  k <- length(Om)
  W <- vector("list", length=k)
  for (e in seq_len(k)) {
    W[[e]] <- exp(cbind(1, X[,2]) %*% Om[[e]])
  }
  W <- do.call(cbind, W)
  apply(W, 1, function(x) sample.int(k, 1, prob=x))
}


generate_data3 <- function(X)
{
  # Should we switch this be dependent on the X variable?
  wt1 <- logit(X[, 2], scale=-8, mid=0.25)
  wt2 <- logit(X[, 2], scale=8, mid=0.75)
  wt3 <- 1 - (wt1 + wt2)
  W <- cbind(wt1, wt2, wt3)
  return(apply(W, 1, function(x) sample.int(3, 1, prob=x)))
}


generate_data <- function(n, B, rho=1, scale=1, mid=0.5)
{
  M <- generate_obs(n, rho=rho)
  experts <- sample_experts(M, scale=scale, mid=mid)
  M <- cbind(M, NA, NA)
  colnames(M) <- c("x", "z", "y", "expert")
  
  for (ii in seq_along(experts)) {
    e <- experts[ii]
    e_sd <- sqrt(exp(B[[e]]['phi']))
    y <- B[[e]]["c"] + B[[e]]["x"] * M[ii, "x"] + B[[e]]["z"] * M[ii, "z"] + rnorm(1, sd=e_sd)
    M[ii, "y"] <- y
    M[ii, "expert"] <- e
  }
  data.frame(M)
}


generate_multinomial <- function(n, B, Om, rho=1)
{
  M <- generate_obs(n, rho=rho)
  experts <- sample_multinomial(M, Om=Om)
  M <- cbind(M, NA, NA)
  colnames(M) <- c("x", "z", "y", "expert")
  
  for (ii in seq_along(experts)) {
    e <- experts[ii]
    e_sd <- sqrt(exp(B[[e]]['phi']))
    y <- B[[e]]["c"] + B[[e]]["x"] * M[ii, "x"] + B[[e]]["z"] * M[ii, "z"] + rnorm(1, sd=e_sd)
    M[ii, "y"] <- y
    M[ii, "expert"] <- e
  }
  data.frame(M)
}


calc_KL_divergence <- function(B, plot_dist=FALSE)
{
  for (i in 1:length(B)) {
    B[[i]]["sd"] <- sqrt(exp(B[[i]]["sd"]))
  }
  means <- lapply(B, get_means)
  if (plot_dist) {
    xx <- seq(-2, 2, 0.005)
    d1 <- dnorm(xx, mean=means[[1]], sd=B[[1]]["sd"])
    d2 <- dnorm(xx, mean=means[[2]], sd=B[[2]]["sd"])
    lb <- min(qnorm(0.005, mean=means[[1]], sd=B[[1]]["sd"]),
              qnorm(0.005, mean=means[[2]], sd=B[[2]]["sd"]))
    ub <- max(qnorm(0.995, mean=means[[1]], sd=B[[1]]["sd"]),
              qnorm(0.995, mean=means[[2]], sd=B[[2]]["sd"]))
    my <- max(c(d1, d2))
    plot(xx, d1, ylim=c(-0.05, my + 0.05), xlim=c(lb, ub), type="l")
    lines(xx, d2)
    grid()
  }
  
  nn <- length(B)
  out <- c()
  nms <- c()
  for (b in 1:(nn - 1)) {
    for (d in (b + 1):nn) {
      t1 <- log(B[[d]]["sd"] / B[[b]]["sd"])
      t2 <- (B[[b]]["sd"]**2 + (means[[b]] - means[[d]])**2) / (2 * B[[d]]["sd"]**2)
      kl <- sum(t1, t2, -0.5)
      out <- c(out, kl)
      nms <- c(nms, paste(c(b, d), collapse="-"))
    }
  }
  names(out) <- nms
  out
}


simulate_2me <- function(seed, n, B, rho, scale, mid)
{
  set.seed(seed)
  dat <- generate_data(n=n, B=B, rho=rho, scale=scale, mid=mid)
  mod <- hme(c("0", "0.1", "0.2"), y ~ x + z | z, data=dat, trace=1L)
  mod$full.vcv <- calculate_sandwich_vcov(mod)
  mod[c("Y", "X", "Z", "scores", "parM")] <- NA
  return(list(seed=seed, model=mod))
}


run_monte_carlo <- function(S=1000, f, output_file, ...)
{
  out <- vector("list", S)
  seeds <- sample.int(1e8, S)
  cnt <- 0L
  for (s in seq_along(seeds)) {
    lst <- try(f(seeds[s], ...))
    if (!inherits(lst, "try-catch")) {
      out[[s]] <- lst
    } else {
      out[[s]] <- list(seed=seeds[s], model=NULL)
    }
    cnt <- cnt + 1L
    if (cnt %% 25 == 0) {
      saveRDS(out, file=output_file)
    }
  }
}


report_result <- function(r)
{
  if (r$variance_test$reject_null) {
    if (r$voung_lr_test$reject_null) {
      return(sign(r$voung_lr_test$pars$LRn))
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}

# --------------- Simulate Data for 2-experts --------------- #

# Sample size
N <- 1000

# Default Parameters
# from the paper phi = log(var) which implies that sd = sqrt(exp(phi))
Betas <- list(c(c=1.2, x=1.5, z=0.5, phi=-3.218876), # 0.2 
              c(c=0.8, x=2.0, z=1.1, phi=-3.429597)) # 0.18

calc_KL_divergence(Betas, TRUE)

DEF_rho <- 1
DEF_scale <- 8
DEF_mid <- 0.5


Mids <- c(0, 0.1, 0.2, 0.3, 0.4)
Scales <- c(125, 25, 4, 2, 1, 0.01)
Rhos <- c(0.125, 0.25, 0.5, 2, 4)
for (r in Rhos) {
  filenm <- paste0(sprintf("models/simulations/Rhos_%d.RDS", as.integer(r * 1000)))
  run_monte_carlo(S=1000, f=simulate_2me, output_file=filenm,
                  N=1000, B=Betas, rho=r, scale=DEF_scale, mid=DEF_mid)
}



r <- 4
file_nm <- paste0(sprintf("models/simulations/Rhos_%d.RDS", as.integer(r * 1000)))
mods <- readRDS(file_nm)
loglikes <- sapply(mods, function(x) logLik(x$model))
hist(loglikes)

seeds <- sapply(mods, `[[`, "seed")

thresh <- 1100
for (ii in seq_along(mods)) {
  m <- mods[[ii]]
  LL <- logLik(m$model)
  if (LL < thresh) {
    newLL <- -Inf
    while (newLL < thresh) {
      newseed <- sample.int(1e8, 1)
      if (newseed %in% seeds)
        next
      newmod <- simulate_2me(seed=newseed, N=N, B=Betas, rho=r, scale=DEF_scale, mid=DEF_mid)
      newLL <- logLik(newmod[["model"]])
      if (newLL > thresh)
        mods[[ii]] <- newmod
    }
  }
}

saveRDS(mods, file=file_nm)

summary(tst)
plot(tst[, 2:1], ylim=c(-4,4), xlim=c(-0.05, 1))

lapply(split(DT, DT$expert), summary)
plot(density(DT[DT$expert == 1, "y"]))
lines(density(DT[DT$expert == 2, "y"]))
lines(density(DT[DT$expert == 3, "y"]))
lines(density(DT[DT$expert == 4, "y"]))


# --------------------------- Standard Error Means --------------------------- #

grab_errs_and_pars <- function(obj, vcv=c("sandwich", "OPG", "hessian"))
{
  vcv <-match.arg(vcv)

  # Grab the specified variance-covariance matrix
  if (vcv == "OPG") {
    std_errs <- sqrt(diag(solve(obj$full.vcv[[vcv]])))
  } else if (vcv == "hessian") {
    std_errs <- sqrt(diag(-solve(obj$full.vcv[[vcv]])))
  } else {
    std_errs <- sqrt(diag(obj$full.vcv[[vcv]]))
  }
  
  
  # the std_errs vector is not broken up as a list of gate or expert nodes
  # so we must do it
  
  # Gate std_errs
  depth <- max(sapply(strsplit(obj[["expert.nms"]], "\\."), length))
  lgpn <- length(obj[["gate.pars.nms"]])
  
  # What kind of MoE is it
  wide_hme <- depth == 2 && length(obj[["expert.nms"]]) > 2
  
  # Wide HME
  if (wide_hme) {
    lgp <- length(obj[["gate.pars"]][[1]])
    ngates <- lgp * lgpn
    gate.errs <- std_errs[1:ngates]
    gate.errs <- list(split(gate.errs, ((seq_along(gate.errs) - 1) %/% lgpn) + 1))
  } else {
    lgp <- length(obj[["gate.pars"]])
    ngates <- lgp * lgpn
    gate.errs <- std_errs[1:ngates]
    gate.errs <- split(gate.errs, ((seq_along(gate.errs) - 1) %/% lgpn) + 1)
  }
  names(gate.errs) <- names(obj[["gate.pars"]])
  
  # Expert std_errs
  lep <- length(obj[["expert.pars"]])
  lepn <- length(obj[["expert.pars.nms"]])
  expert.errs <- std_errs[(ngates + 1):length(std_errs)]
  expert.errs <- split(expert.errs, ((seq_along(expert.errs) - 1) %/% lepn))
  names(expert.errs) <- names(obj[["expert.pars"]])
  
  
  pars <- vector("list", length(obj[["expert.pars"]]))
  gpps <- vector("numeric", length(pars))
  names(gpps) <- names(pars) <- obj[["expert.nms"]]
  
  for (e in names(gpps)) {
    #gpp <- gate_path_product("0", e, obj[["list_priors"]])
    pars[[e]] <- obj[["expert.pars"]][[e]]
    gpps[e] <- obj[["expert.pars"]][[e]][1]
  }
  expert_order <- order_betas(pars)
  pars <- unlist(pars[expert_order])
  expert.errs <- unlist(expert.errs[expert_order])
  names(pars) <- names(expert.errs) <- c("c1", "x1", "z1", "sd1",
                                         "c2", "x2", "z2", "sd2")
  return(list(pars=pars, std_errs=expert.errs))
}


grab_all <- function(lst, vcv=c("sandwich", "OPG", "hessian"))
{
  vcv <- match.arg(vcv)
  values <- lapply(lst, function(x) grab_errs_and_pars(x$model, vcv=vcv))
  pars <- do.call(rbind, lapply(values, `[[`, "pars"))
  pars <- cbind(t(data.frame(apply(pars, 2, sd))), data.frame(type="Par Std Dev"))
  std_errs <- do.call(rbind, lapply(values, `[[`, "std_errs"))
  std_errs <- cbind(t(data.frame(colMeans(std_errs))), data.frame(type="Avg Std Err"))
  return(rbind(pars, std_errs))
}


VCV <- c("sandwich", "OPG", "hessian")
Metric <- c("Mid", "Scales", "Rhos")
dir_ <- "models/simulations/monte_carlo/"
collect <- data.frame()
for (m in Metric) {
  for (v in VCV) {
    files_ <- list.files(dir_, pattern=m)
    lst_sims <- lapply(paste0(dir_, files_), readRDS)
    names(lst_sims) <- files_
    tst <- lapply(lst_sims, grab_all,  vcv=v)
    normalization <- switch(m, "Rhos"=1000, "Mid"=10, "Scales"=100)
    values <- as.numeric(sapply(strsplit(files_, split="_|\\."), `[[`, 2))
    tst <- tst[order(values)]
    tst <- do.call(rbind, tst)
    tst$par_value <- rep(values[order(values)], each=2)
    tst$norm <- normalization
    tst$Metric <- m
    tst$VCV <- v
    row.names(tst) <- NULL
    collect <- rbind(collect, tst)
  }
}
collect2 <- data.table(collect)
collect2[, VCV := as.factor(VCV)]
collect2[, type := as.factor(type)]
collect2[, x_label := par_value / norm]
collect2[, x_label := ifelse(x_label == as.integer(x_label),
                             as.character(as.integer(x_label)),
                             as.character(x_label))]

# Rhos -----------------------------------------------
collect_rho <- collect2[Metric == "Rhos"]
collect_rho[, x_label := factor(x_label, levels=collect_rho[, unique(x_label)])]
collect_rho[, norm := NULL]
collect_rho[, Metric := NULL]
collect_rho[, par_value := NULL]
collect_rho <- melt(collect_rho)
collect_rho <- unique(collect_rho, by="value")

ggplot(collect_rho, aes(x=x_label, y=value, col=type)) +
  geom_point(aes(shape=VCV), position="jitter") +
  facet_wrap(~ variable, nrow=2, ncol=4) +
  scale_shape_manual(values=c(0, 1, 2)) +
  labs(title="Average Simulated Standard Error -vs- Simulated Parameter Standard Deviation",
       subtitle=expression(paste("Association Parameter: ", rho)),
       y="", x=expression(paste(rho)), col="", shape="")


# Midpoints ------------------------------
collect_mid <- collect2[Metric == "Mid"]
collect_mid[, x_label := factor(x_label, levels=collect_mid[, unique(x_label)])]
collect_mid[, norm := NULL]
collect_mid[, Metric := NULL]
collect_mid[, par_value := NULL]
collect_mid <- melt(collect_mid)
collect_mid <- unique(collect_mid, by="value")

ggplot(collect_mid, aes(x=x_label, y=value, col=type)) +
  geom_point(aes(shape=VCV), position="jitter") +
  facet_wrap(~ variable, nrow=2, ncol=4) +
  scale_shape_manual(values=c(0, 1, 2)) +
  labs(title="Average Simulated Standard Error -vs- Parameter Standard Deviation",
       subtitle="Probability Parameter: c",
       y="", x="c", shape="", col="")


# Slope ------------------------------
collect_scale <- collect2[Metric == "Scales"]
collect_scale[, x_label := factor(x_label, levels=collect_scale[, unique(x_label)])]
collect_scale[, norm := NULL]
collect_scale[, Metric := NULL]
collect_scale[, par_value := NULL]
collect_scale <- melt(collect_scale)
collect_scale <- unique(collect_scale, by="value")

ggplot(collect_scale, aes(x=x_label, y=value, col=type)) +
  geom_point(aes(shape=VCV), position="jitter") +
  facet_wrap(~ variable, nrow=2, ncol=4) +
  scale_shape_manual(values=c(0, 1, 2)) +
  labs(title="Average Simulated Standard Error -vs- Parameter Standard Deviation",
       subtitle=expression(paste("Probability Parameter: ", alpha)),
       y="", x=expression(paste(alpha)), shape="", col="")




# ------------------------------ Coverage Ratio ------------------------------ #

grab_coverage_ratio <- function(lst, B, vcv=c("sandwich", "OPG", "hessian"))
{
  vcv <- match.arg(vcv)
  values <- lapply(lst, function(x) grab_errs_and_pars(x$model, vcv=vcv))
  pars <- do.call(rbind, lapply(values, `[[`, "pars"))
  std_errs <- do.call(rbind, lapply(values, `[[`, "std_errs"))
  B <- unlist(B[order_betas(B)])
  
  out <- matrix(FALSE, nrow=nrow(pars), ncol=ncol(pars))
  for (ii in seq_len(nrow(pars))) {
    ub <- B + qnorm(0.975) * std_errs[ii, ]
    lb <- B - qnorm(0.975) * std_errs[ii, ]
    out[ii, ] <- pars[ii, ] < lb | pars[ii, ] > ub
  }
  out <- data.frame(t(colMeans(out)))
  names(out) <- c("c1", "x1", "z1", "sd1", "c2", "x2", "z2", "sd2")
  return(out)
}

VCV <- c("sandwich", "OPG", "hessian")
Metric <- c("Mid", "Scales", "Rhos")
dir_ <- "models/simulations/monte_carlo/"
collect <- data.frame()
for (m in Metric) {
  for (v in VCV) {
    files_ <- list.files(dir_, pattern=m)
    lst_sims <- lapply(paste0(dir_, files_), readRDS)
    names(lst_sims) <- files_
    tst <- lapply(lst_sims, grab_coverage_ratio,  B=Betas, vcv=v)
    normalization <- switch (m, "Rhos"=1000, "Mid"=10, "Scales"=100)
    values <- as.numeric(sapply(strsplit(files_, split="_|\\."), `[[`, 2))
    tst <- tst[order(values)]
    tst <- do.call(rbind, tst)
    tst$par_value <- values[order(values)]
    tst$norm <- normalization
    tst$Metric <- m
    tst$VCV <- v
    row.names(tst) <- NULL
    collect <- rbind(collect, tst)
  }
}



# --------------- Graphs --------------- #

par(mfrow=c(1,2))
# Node Purity: Scale parameter affect the probabilities of 
xx <- seq(-0.1, 1.1, .005)
Alphas <- c(125, 25, 8, 4, 2, 1, 0.01)
plot(xx, logit(xx, scale=Alphas[1], mid=0.5), type="l",
     xlab="z", ylab="P(z)", main="Slope parameter")
for (jj in Alphas[-1]) {
  lwd_ <- ifelse(jj == DEF_scale, 2, 1)
  lines(xx, logit(xx, scale=jj, mid=0.5), lwd=lwd_)
}
grid()


# Relative Expert Size: Mid paramter can be shifted to create groups of 
# different sizes
Mids <- seq(0, 0.5, 0.1)
plot(xx, logit(xx, scale=DEF_scale, mid=Mids[length(Mids)]), type="l",
     xlab="z", ylab="P(z)", main="Mid-point", lwd=2)
for (m in Mids[-length(Mids)]) {
  lines(xx, logit(xx, scale=DEF_scale, mid=m))
}
grid()



# Correlation Strength: Rho parameter effects tightness of relationship between
# x and z
par(mfrow=c(2, 3))
Rhos <- c(0.125, 0.25, 0.5, 1, 2, 4)
for (r in Rhos) {
  tst <- generate_data(N, B=Betas, rho=r, scale=DEF_scale, mid=DEF_mid)
  plot(tst[, 2:1], ylim=c(-4,4), xlim=c(-0.05, 1), main=bquote(rho == .(r)))
  grid()
}


# -------------- Simulate Data for Expert Over-expert Evaluation -------------- #

Betas <- list(c(c=1.2, x=1.5, z=0.5, phi=-3.218876), # 0.2 
              c(c=0.8, x=2.0, z=1.1, phi=-3.429597), # 0.18
              c(c=0.4, x=0.8, z=2.0, phi=-3.543914), # 0.17
              c(c=1.5, x=0.5, z=0.8, phi=-3.321462)) # 0.19
Omega <- list(c(-4, 6),
              c(2, -6),
              c(6, -18),
              c(-15, 20))

N <- 4
DT <- generate_multinomial(2000, Betas[1:N], Omega[1:N], rho=DEF_rho)
with(DT, plot(z, x, col=expert))
legend("bottomright", legend=paste0("E", 1:4), col=1:4, pch=1)
table(DT$expert)


tstHME2 <- hme(c("0", "0.1", "0.2"),  y ~ x + z | z, data=DT)
tstHME2$full.vcv <- calculate_sandwich_vcov(tstHME2)

st_betas <- Betas[c(1, 3, 4)]
names(st_betas) <- c("0.1", "0.2", "0.3")
tstHME3w <- hme(c("0", "0.1", "0.2", "0.3"), y ~ x + z | z, data=DT,
                init_expert_pars = st_betas, maxiter=10)
tstHME3w$full.vcv <- calculate_sandwich_vcov(tstHME3w)

names(st_betas) <- c("0.2", "0.1.1", "0.1.2")
tstHME3d <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"), y ~ x + z | z, data=DT,
                init_expert_pars = st_betas)
tstHME3d$full.vcv <- calculate_sandwich_vcov(tstHME3d)


st_betas <- Betas[c(1, 2, 4, 3)]
names(st_betas) <- c("0.1", "0.2", "0.3", "0.4")
tstHME4w <- hme(c("0", "0.1", "0.2", "0.3", "0.4"),  y ~ x + z | z, data=DT,
                init_expert_pars = st_betas, maxiter=3)
tstHME4w$full.vcv <- calculate_sandwich_vcov(tstHME4w)

names(st_betas) <- c("0.1.1", "0.1.2", "0.2.1", "0.2.2")
tstHME4d <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2", "0.2.1", "0.2.2"),
                y ~ x + z | z, data=DT, init_expert_pars = st_betas)
tstHME4d$full.vcv <- calculate_sandwich_vcov(tstHME4d)


modnms <- ls(pattern="tstHME")
R <- matrix(NA, nrow=5, ncol=5, dimnames=list(modnms, modnms))

for (m in modnms) {
  mod1 <- get(m)
  for (n in modnms) {
    if (m == n)
      next
    mod2 <- get(n)
    vs <- voung_selection(mod1, mod2)
    r <- report_result(vs)
    R[m, n] <- r
    R[n, m] <- r * -1
  }
}

R

for (m in ls(pattern="HME")) {
  saveRDS(get(m), file=paste0("models/simulations/cross_experts/", sprintf("%s_%de.RDS", m, N)))
}


# -------------- Graph the cross-expert results -------------- #

modnms <- c("tstHME2_2e.RDS",  "tstHME2_3e.RDS", "tstHME2_4e.RDS")
mods <- lapply(modnms, function(x) readRDS(paste0("models/simulations/cross_experts/", x)))
names(mods) <- modnms

plot_densities <- function(X, tit="")
{
  dts <- split(X, f=X$expert)
  D <- lapply(dts, function(x) density(x$y))
  N <- length(D)
  yrange <- range(unlist(lapply(D, function(x) x$y)))
  xrange <- range(unlist(lapply(D, function(x) x$x)))
  yrange <- c(0, 1.02 * yrange[2])
  plot(density(X[X$expert == 1, "y"]), ylim=yrange, xlim=xrange, main=tit,
       xlab=sprintf("%d - Sub-populations", N), ylab="", col=1)
  for (ii in 2:length(D)) {
    lines(density(X[X$expert == ii, "y"]), col=ii)
  }
  grid()
}

par(mfcol=c(2, 3))
for (n in 2:4) {
  DT <- generate_multinomial(2000, Betas[1:n], Omega[1:n], rho=DEF_rho)
  with(DT, plot(z, x, col=expert, pch=1:n, cex=0.25, main="(X, Z) Joint Dist.")); grid()
  plot_densities(DT, "Density: y")
}



legend("bottomright", legend=paste0("E", 1:4), col=1:4, pch=1)



