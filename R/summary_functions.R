summary.hme <- function(obj)
{  
  # cat(sprintf("Log-Likelihood: %f", log_Like(obj)))
  for (g in names(obj[["gate.pars"]])) {
    
    pars <- obj[["gate.pars"]][[g]]
    if (!is(pars, "list"))
      pars <- list(pars)
    info <- obj[["gate.info.matrix"]]
    if (!is(info, "list"))
      info <- list(info)
    
    nn <- length(pars)
    cat(sprintf("\n--------------------------\ngate-%s\n", g))
    for (i in seq_len(nn)) {
      p <- pars[[i]]
      se <- sqrt(diag(info[[i]][["sandwich"]]))
      zstat <- p / se
      tbl <- cbind(p, se, zstat, 2 * pnorm(-abs(zstat)))
      colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
      printCoefmat(tbl, digits = 3)
      
      # zscore <- 
      #colnames(p) <- list("Estimate Std.", "Std. Error") # "t-value", "Pr(>|t|)")
      #print(p)
      cat("\n")
    }
  }
  
  cat("\n\n----------------------------------------------------\n\n")
  
  for (e in names(obj[["expert.pars"]])) {
    p <- obj[["expert.pars"]][[e]]
    
    p <- as.matrix(round(c(p, exp(p[length(p)])), 4))
    dimnames(p) <- list(c(obj[["expert.pars.nms"]], ""), "")
    cat(sprintf("\n--------------------------\nexpert-%s\n", e))
    print(p)
  }
  
  cat("\n\n----------------------------------------------------\n\n")
  
  for (e in names(obj[["expert.pars"]])) {
    p <- obj[["gate.margins"]][[e]]
    
    p <- as.matrix(round(p, 4))
    dimnames(p) <- list(c(obj[["gate.pars.nms"]]), "")
    cat(sprintf("\n--------------------------\nexpert-%s marginal effects\n", e))
    print(p)
  }
}

#debugonce(summary.hme)
#summary(tst)