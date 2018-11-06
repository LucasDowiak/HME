source("~/hme/R/benchmark_hme.R")
setwd("~/hme/")

library(data.table)


dtf <- fread("wage_equation.csv")

hme_wage_1D_2E <- readRDS("wage/hme_wage_1D_2E.RDS")
hme_wage_2D_3E <- readRDS("wage/hme_wage_2D_3E.RDS")
hme_wage_2D_4E <- readRDS("wage/hme_wage_2D_4E.RDS")

lapply(list(hme_wage_1D_2E, hme_wage_2D_3E, hme_wage_2D_4E, hme3), criterion, "aic")
saveRDS(hme3, file="wage/hme_1D_3E.RDS")
form <- "lnwage ~ age16 + age16sq + yreduc | age16 + age16sq + black + indian + asian + hisp + yreduc"
c("0", "0.1", "0.2")
c("0", "0.1", "0.2", "0.1.1", "0.1.2")
c("0", "0.1", "0.2", "0.1.1", "0.1.2", "0.2.1", "0.2.2")
c("0",
  "0.1", "0.2",
  "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
  "0.1.1.1", "0.1.1.2")
c("0",
  "0.1", "0.2",
  "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
  "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2")

"lnwage ~ age16 + age16sq | age16 + age16sq + black + indian + asian + hisp + yreduc"

lm_mod <- lm("lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp",
             data=dtf[sex==0])


debugonce(hme)
hme1 <- hme(c("0",
              "0.1", "0.2",
              "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2"),
            "lnwage ~ age16 + age16sq + black + indian + asian + hisp + yreduc | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
            #init_gate_pars = hme1$gate.pars, init_expert_pars = hme1$expert.pars)
ME1 <- marginal_effects(hme1)
hme2 <- hme(c("0",
              "0.1", "0.2",
              "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2"),
            "lnwage ~ age16 + age16sq + black + indian + asian + hisp + yreduc | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
            #init_gate_pars = hme2$gate.pars, init_expert_pars = hme2$expert.pars)
ME2 <- marginal_effects(hme2)
hme3 <- hme(c("0",
              "0.1", "0.2",
              "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2"),
            "lnwage ~ age16 + age16sq + black + indian + asian + hisp + yreduc | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
            # init_gate_pars = hme3$gate.pars, init_expert_pars = hme3$expert.pars)
ME3 <- marginal_effects(hme3)
hme4 <- hme(c("0",
              "0.1", "0.2",
              "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2"),
            "lnwage ~ age16 + age16sq + black + indian + asian + hisp + yreduc | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
            # init_gate_pars = hme4$gate.pars, init_expert_pars = hme4$expert.pars)
ME4 <- marginal_effects(hme4)
hme5 <- hme(c("0",
              "0.1", "0.2",
              "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2"),
            "lnwage ~ age16 + age16sq + black + indian + asian + hisp + yreduc | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
            # init_gate_pars = hme5$gate.pars, init_expert_pars = hme5$expert.pars)
ME5 <- marginal_effects(hme5)


hme2s <- readRDS("hme_wage_1D_2E.RDS")
hme3d <- readRDS("hme_wage_2D_3E.RDS")
hme3w <- readRDS("hme_1D_3E.RDS")
hme4s <- readRDS("hme_wage_2D_4E.RDS")



sapply(list(hme2, hme3d, hme3w, hme4), logLik)
sapply(list(hme1, hme2s, hme3d, hme3w, hme4s), criterion, "bic")

sapply(list(hme1, hme2, hme3, hme4, hme5), logLik)

sapply(list(hme1, hme2, hme3, hme4, hme5), logLik, "bic")
lapply(list(hme1, hme2, hme3, hme4, hme5), `[[`, "expert.pars")

lapply(list(ME1, ME2, ME3, ME4, ME5), `[[`, "gate_expert_margins")


debugonce(grow_the_tree)
tst <- grow_the_tree(hme_wage_3D)



grid_ <- seq(0, 5.64, 0.05)
expert011 <- hme_wage_2D_3E$expert.pars[["0.1.1"]][[1]] + hme_wage_2D_3E$expert.pars[["0.1.1"]][[2]] * grid_ +  hme_wage_2D_3E$expert.pars[["0.1.1"]][[3]] * grid_**2
expert012 <- hme_wage_2D_3E$expert.pars[["0.1.2"]][[1]] + hme_wage_2D_3E$expert.pars[["0.1.2"]][[2]] * grid_ +  hme_wage_2D_3E$expert.pars[["0.1.2"]][[3]] * grid_**2
expert02 <- hme_wage_2D_3E$expert.pars[["0.2"]][[1]] + hme_wage_2D_3E$expert.pars[["0.2"]][[2]] * grid_ +  hme_wage_2D_3E$expert.pars[["0.2"]][[3]] * grid_**2

plot(grid_, expert011, lwd=2, type="l", ylim=c(0, 3), col="blue")
lines(grid_, expert012, lwd=2, col="orange")
lines(grid_, expert02, lwd=2, col="green")
grid()



for (f_ in dir(pattern="min")) {
  rds <- readRDS(f_)
  cat(f_, "\n", logLik(rds), "\n", criterion(rds, "aic") / rds$N, "\n")
  #print(marginal_effects(rds)) 
}


plot_histograms <- function(obj)
{
  nms <- obj$expert.nms
  for (nn in nms) {
    gpp <- gate_path_product("0", nn, obj$list_priors)
    share <- sum(gpp) / obj$N
    tits <- sprintf("%s - Share: %.3f", nn, share)
    hist(gpp, main=tits, ylab="")
  }
}

par(mfrow=c(3,2))
plot_histograms(readRDS("wages/hme_3D_5E_mid.RDS"))

par(mfrow=c(2,2))
plot_histograms(readRDS("wages/hme_2D_4E_mid.RDS"))

par(mfrow=c(2,2))
plot_histograms(readRDS("wages/hme_2D_3E_mid.RDS"))

hme3e <- readRDS("wages/hme_2D_3E_mid.RDS")
age_ <- seq(0, 49, .25)
p1 <- hme3e$expert.pars[["0.1.1"]]
p2 <- hme3e$expert.pars[["0.1.2"]]
p3 <- hme3e$expert.pars[["0.2"]]
proj1 <- p1[1] + p1[2] * age_ + p1[3] * age_**2 + p1[4] * 14
proj2 <- p2[1] + p2[2] * age_ + p2[3] * age_**2 + p2[4] * 14
proj3 <- p3[1] + p3[2] * age_ + p3[3] * age_**2 + p3[4] * 14

plot(age_, proj2, ylim=c(0, 4), type="l", col="blue", lwd=2,
     main="Log Wage Equation (Yrs Educ = 14)", ylab="", xlab="age16")
lines(age_, proj1, col="red", lwd=2)
lines(age_, proj3, col="purple", lwd=3)
grid()
for (ee in hme3e$expert.nms) {
  p <- hme3e$expert.pars[[ee]]
  proj <- p[1] + p[2] * age_ + p[3] * age_**2 + p[4] * 14
  print(summary(proj))
}
