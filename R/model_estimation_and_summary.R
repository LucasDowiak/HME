Trees <- list(#D1_E2=c("0", "0.1", "0.2"),
              D1_E3=c("0", "0.1", "0.2", "0.3"),
              D1_E3=c("0", "0.1", "0.2", "0.3", "0.4")
              #D2_E3=c("0", "0.1", "0.2", "0.1.1", "0.1.2"),
              #D2_E4=c("0", "0.1", "0.2", "0.1.1", "0.1.2", "0.2.1", "0.2.2"),
              #3_E5=c("0",
              #        "0.1", "0.2",
              #        "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              #        "0.1.1.1", "0.1.1.2"),
              #D3_E6=c("0",
              #        "0.1", "0.2",
              #        "0.1.1", "0.1.2", "0.2.1", "0.2.2", 
              #        "0.1.1.1", "0.1.1.2", "0.1.2.1", "0.1.2.2")
)


for (tree in names(Trees)) {
  
  Hmes <- vector("list", 5)
  for (ii in seq_along(Hmes)) {
    
    again <- TRUE
    while(again) {
      mod <- try(hme(c("0", "0.1", "0.2", "0.3"),# Trees[[tree]],
                     "lnwage ~ age16 + age16sq + yreduc | age16 + age16sq + yreduc + sex + black + indian + asian + hisp + Creativity + Design + Analytics + Perseptive",
                     data=dtftrain, holdout=dtftest, maxiter=20, tolerance=1e-3, trace=1))
      
      if (inherits(mod, "hme")) {
        again <- FALSE
      }
    }
    
    Hmes[[ii]] <- mod
    
  }
  
  # Which of the models maximizes the likelihood function
  #    >  If it's been refactored, make sure the methods still work accrtly
  #        > Save the file
  idxLL <- which.max(sapply(Hmes, function(x) max(x$logL)))
  mod <- Hmes[[idxLL]]
  if (which.max(mod[["logL"]]) != length(mod[["logL"]])) {
    mod <- refactor_hme(mod, "loglik")
  }
  file_name <- sprintf("models/hme_%s_Onet_full_ll.RDS", tree)
  saveRDS(mod, file=file_name)
  
  # Which of the models mimimize the likelihood function
  #    >  If it's been refactored, make sure the methods still work accrtly
  #        > Save the file
  idxMSE <- which.max(sapply(Hmes, function(x) min(x$MSE)))
  mod <- Hmes[[idxMSE]]
  if (which.min(mod[["MSE"]]) != length(mod[["MSE"]])) {
    mod <- refactor_hme(mod, "mse")
  }
  file_name <- sprintf("models/hme_%s_Onet_full_mse.RDS", tree)
  saveRDS(mod, file=file_name)
}



model_metric <- function(file_)
{
  bkts <- unlist(strsplit(file_, "_"))
  mod <- readRDS(file_)
  modtyp <- bkts[5]
  crityp <- gsub("\\.RDS", "\\1", bkts[6])
  depth <- as.numeric(gsub("^.*([0-9]+)", "\\1", bkts[2]))
  experts <- as.numeric(gsub("^.*([0-9]+)", "\\1", bkts[3]))
  metrics <- unlist(napply(c("aic", "bic", "mse"),
                           function(x) criterion(mod, x)))
  
  data.frame(model=modtyp,
             depth=depth,
             experts=experts,
             criterion=crityp,
             logl=logLik(mod) / mod$N,
             aic=metrics[1] / mod$N,
             bic=metrics[2] / mod$N,
             mse=metrics[3]
             )
}

yrsedu_coef <- function(file_)
{
  bkts <- unlist(strsplit(file_, "_"))
  mod <- readRDS(file_)
  modtyp <- bkts[5]
  crityp <- gsub("\\.RDS", "\\1", bkts[6])
  depth <- as.numeric(gsub("^.*([0-9]+)", "\\1", bkts[2]))
  experts <- as.numeric(gsub("^.*([0-9]+)", "\\1", bkts[3]))
  
  if (modtyp == "min") {
    margeffs <- marginal_effects(mod)$gate_margins
  } else {
    margeffs <- marginal_effects(mod)$full_margins
  }
  
  beta <- margeffs[which(names(margeffs) == "yreduc")]
  if (length(beta) == 0)
    beta <- NA_real_
  
  data.frame(model=modtyp,
             depth=depth,
             experts=experts,
             criterion=crityp,
             beta=beta
  )
}


tsts <- dir("models/", pattern="^hme_", full.names = T)
mods <- lapply(tsts, readRDS)

M <- do.call(rbind, lapply(tsts, yrsedu_coef))
M <- with(M, M[order(model, criterion), ])
#MM <- with(M, M[M$criterion == "ll", ][order(model, criterion), ])
split(M, M$model)

tt <- readRDS(tsts[1])
