setwd("~/Git/hme/")
source("R/benchmark_hme.R")
source("R/wage_data_prep.R")



all_vars <- c("lnwage", "age16", "age16sq", "sex", "black", "indian", "asian", "hisp",
              "yreduc", "Creativity", "Design", "Analytics", "Perceptive")

boolcc <- complete.cases(dtf[, .SD, .SDcols=all_vars[-10]])
booldd <- complete.cases(dtf[, .SD, .SDcols=all_vars])
set.seed(433)
booltest <- as.logical(rbinom(dtf[, .N], 1, 0.1))
dtftest <- dtf[booltest & booldd]
dtftrain <- dtf[!booltest & booldd]
# seperate into test and validation set


form_ols <- "lnwage ~ age16 + age16sq + yreduc + sex + black + indian + asian + hisp + Creativity + Design + Analytics + Perceptive"

form_full <- "lnwage ~ age16 + age16sq + yreduc + sex + black + indian + asian + hisp + Creativity + Design + Analytics + Perceptive | age16 + age16sq + yreduc + sex + black + indian + asian + hisp + Creativity + Design + Analytics + Perceptive"

form_mid <- "lnwage ~ age16 + age16sq + yreduc | age16 + age16sq + yreduc + sex + black + indian + asian + hisp + Creativity + Design + Analytics + Perceptive"

form_min <- "lnwage ~ age16 + age16sq | age16 + age16sq + yreduc + sex + black + indian + asian + hisp + Creativity + Design + Analytics + Perceptive"

logLik(lm(form_ols, data=dtftrain))

dir("models/")
mod3d <- readRDS("models/hme_3d_full.RDS")
B <- 199
b <- 1
ME_mat <- matrix(0, B, 13)
while (b <= B) {
  t1 <- Sys.time()
  print(sprintf("Starting %s out of %s. Time is %s", b, B, Sys.time()))
  dtfshuffle <- dtftrain[sample.int(.N, replace=TRUE)]
  mod <- try(
    hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"), formula=form_full,
        data=dtfshuffle, holdout=dtftest, maxiter=150, tolerance=1e-4, trace=0)
  )
  if (inherits(mod, "try-error")) {
    next
  } else {
    me <- marginal_effects(mod)
    ME_mat[b,] <- me$margins
    t2 <- Sys.time()
    print(sprintf("Ending %s out of %s. It took %.3f minutes.", b, B, t2 - t1))
    write.csv(ME_mat, "margins/min_spec_3w_me.csv", row.names=FALSE)
    b <- b + 1
  }
}
colnames(ME_mat) <- names(me$margins)
write.csv(ME_mat, "margins/min_spec_3w_me.csv", row.names=FALSE)



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

lm_mod <- lm("lnwage ~ age16 + age16sq + yreduc + black + indian + asian + hisp",
             data=dtftrain)


debugonce(hme)

# Wage
hme_2d <- hme(c("0", "0.1", "0.2"), form=form_full,
              data=iris, holdout=NULL, maxiter=1, tolerance=1e-5, trace=0,
              init_gate_pars=hme_iris_2d$gate.pars, init_expert_pars=hme_iris_2d$expert.pars)

hme_3d <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"), form=form_full,
              data=dtftrain, holdout=dtftest, maxiter=100, tolerance=1e-5, trace=0)

hme_3w <- hme(c("0", "0.1", "0.2", "0.3"), form=form_full,
              data=dtftrain, holdout=dtftest, maxiter=100, tolerance=1e-5, trace=0)


# Iris
hme_iris_2d <- hme(c("0", "0.1", "0.2"), form="Sepal.Width ~ Petal.Width | Petal.Width",
                   data=iris, maxiter=75, tolerance=1e-5, trace=0)
hme_iris_2d$full.vcv <- calculate_sandwich_vcov(hme_iris_2d)

hme_iris_3d <- hme(c("0", "0.1", "0.2", "0.1.1", "0.1.2"), form="Sepal.Width ~ Petal.Width | Petal.Width",
                   data=iris, maxiter=250, tolerance=1e-6, trace=0)
hme_iris_3d$full.vcv <- calculate_sandwich_vcov(hme_iris_3d)

hme_iris_3w <- hme(c("0", "0.1", "0.2", "0.3"), form="Sepal.Width ~ Petal.Width | Petal.Width",
                   data=iris, maxiter=250, tolerance=1e-6, trace=0)
hme_iris_3w$full.vcv <- calculate_sandwich_vcov(hme_iris_3w)

ME1 <- marginal_effects(hme1)
# saveRDS(hme1, file="models/*d.RDS")
hme2 <- hme(c("0", "0.1", "0.2", "0.3"),
            formula=form_full,
            data=dtftrain, holdout=dtftest, maxiter=1000, tolerance=1e-3, trace=1)
            #init_gate_pars = hme2$gate.pars, init_expert_pars = hme2$expert.pars)
ME2 <- marginal_effects(hme2)


hme3 <- hme(c("0",
             "0.1", "0.2", "0.3"),
             formula = form_full,
             data=dtftrain, holdout=dtftest, maxiter=1000, tolerance=1e-4, trace=0)
             # init_gate_pars = hme3$gate.pars, init_expert_pars = hme3$expert.pars)
hme3$full.vcv <- calculate_sandwich_vcov(hme3)
vs <- voung_selection(hme3, readRDS("models/hme_5d_full_v2.RDS"))

str(vs)
str(hme3$call_)
saveRDS(hme3, file="models/hme_3w_mid_v2.RDS")
ME3 <- marginal_effects(hme3)

hme4c <- readRDS("models/hme_4w_full_v2.RDS")
hme4 <- hme(c("0",
              "0.1", "0.2", "0.3", "0.4"),
            form_full,
            data=dtftrain, holdout=dtftest, maxiter=75, tolerance=1e-5, trace=0)
            # init_gate_pars = hme4c$gate.pars, init_expert_pars = hme4c$expert.pars)
hme4$full.vcv <- calculate_sandwich_vcov(hme4)
vs <- voung_selection(hme4, hme4c)
str(vs)
str(hme4$call_)
saveRDS(hme4, file="models/hme_4w_full_v2.RDS")

ME4 <- marginal_effects(hme4)


hme5 <- hme(c("0",
              "0.1", "0.2", "0.3", "0.4", "0.5"),
            form_min,
            data=dtftrain, holdout=dtftest, maxiter=100, tolerance=1e-5, trace=0)
hme5$full.vcv <- calculate_sandwich_vcov(hme5)

ME5 <- marginal_effects(hme5)




