summary.hme <- function(obj, type=c("experts", "gates", "all"))
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
  
  # Grab the robust standard errors
  std_errs <- sqrt(diag(obj$full.vcv$sandwich))
  
  # the std_errs vector is not broken up as a list of gate or expert nodes so we must do it
  
  # Gate std_errs
  depth <- max(sapply(strsplit(obj[["expert.nms"]], "\\."), length))
  lgpn <- length(obj[["gate.pars.nms"]])
  if (depth == 2 && length(obj[["expert.nms"]]) > 2) {
    lgp <- length(mod[["gate.pars"]][[1]])
    ngates <- lgp * lgpn
    gate.errs <- std_errs[1:ngates]
    gate.pars <- list(split(gate.errs, ((seq_along(gate.errs) - 1) %/% lgpn) + 1))
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
  
  
  cat(sprintf("Log-Likelihood: %f", logLik(obj)))
  if (expertbool) {
    #info <- obj[["expert.info.matrix"]]
    
    #lep <- length(obj[["expert.pars"]])
    #lepn <- length(obj[["expert.pars.nms"]])
    #expert.pars <- pars[(ngates + 1):length(pars)]
    #expert.pars <- split(expert.pars, ((seq_along(expert.pars) - 1) %/% lepn))
    
    for (e in names(obj[["expert.pars"]])) {
      
      gpp <- gate_path_product("0", e, obj$list_priors)
      share <- sum(gpp) / obj$N
      
      cat(sprintf("\n--------------------------\nexpert-%s \t share:  %.3f\n\n", e, share))
      
      p <- obj[["expert.pars"]][[e]]
      se <- expert.errs[[e]]
      zstat <- p / se
      tbl <- cbind(p, se, zstat, 2 * pnorm(-abs(zstat)))
      dimnames(tbl) <- list(obj[['expert.pars.nms']],
                            c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)"))
      printCoefmat(tbl, digits = 3)
    }
  }

  cat("\n\n----------------------------------------------------\n\n")
  
  if (gatebool) {
    for (g in names(obj[["gate.pars"]])) {
      
      pars <- obj[["gate.pars"]][[g]]
      if (!is(pars, "list"))
        pars <- list(pars)
      
      nn <- length(pars)
      cat(sprintf("\n--------------------------\ngate-%s\n", g))
      for (i in seq_len(nn)) {
        p <- pars[[i]]
        se <- gate.errs[[i]]
        zstat <- p / se
        tbl <- cbind(p, se, zstat, 2 * pnorm(-abs(zstat)))
        dimnames(tbl) <- list(obj[['gate.pars.nms']],
                              c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)"))
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
