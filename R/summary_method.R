summary.hme <- function(obj, type=c("all", "gates", "experts"))
{
  type <- match.arg(type)
  gatebool <- expertbool <- FALSE
  if (type == "all") {
    gatebool <- expertbool  <- TRUE
  } else if (type == "gates") {
    gatebool <- TRUE
  } else if (type == "experts") {
    expertbool <- TRUE
  }
  
  cat(sprintf("Log-Likelihood: %f", logLik(obj)))
  if (expertbool) {
    info <- obj[["expert.info.matrix"]]
    for (e in names(obj[["expert.pars"]])) {
      cat(sprintf("\n--------------------------\nexpert-%s\n\n", e))
      p <- obj[["expert.pars"]][[e]]
      se <- sqrt(diag(info[[e]][["sandwich"]]))
      zstat <- p / se
      tbl <- cbind(p, se, zstat, 2 * pnorm(-abs(zstat)))
      colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
      printCoefmat(tbl, digits = 3)
    }
  }

  cat("\n\n----------------------------------------------------\n\n")
  
  if (gatebool) {
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
        cat("\n")
      }
    }
  }
  
  cat("\n\n----------------------------------------------------\n\n")
  
  # Seperate the calculation of the marginal effects outside the main function call
  #if (margbool) {
  #  for (e in names(obj[["expert.pars"]])) {
  #    p <- as.matrix(round(obj[["gate.margins"]][[e]], 4))
  #    dimnames(p) <- list(c(obj[["gate.pars.nms"]]), "")
  #    cat(sprintf("\n--------------------------\nexpert-%s marginal effects\n", e))
  #    print(p)
  #  }
  # }
}
