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
c("0", "0.1", "0.2", "0.1.1", "0.1.2", "0.2.1", "0.2.2")

debugonce(hme)
hme1 <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"),
            "lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1, #)
            init_gate_pars = hme1$gate.pars, init_expert_pars = hme1$expert.pars)
ME1 <- marginal_effects(hme1)
hme2 <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"),
            "lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
ME2 <- marginal_effects(hme2)
hme3 <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"),
            "lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
ME3 <- marginal_effects(hme3)
hme4 <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"),
            "lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
ME4 <- marginal_effects(hme4)
hme5 <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"),
            "lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp | age16 + age16sq + black + indian + asian + hisp + yreduc",
            data=dtf[sex==0], maxiter=400, tolerance=1e-3, trace=1)
ME5 <- marginal_effects(hme5)


napply(list(hme1, hme_wage_3D, hme_wage_2D), criterion)


sapply(list(hme1, hme2, hme3, hme4, hme5), logLik)
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
