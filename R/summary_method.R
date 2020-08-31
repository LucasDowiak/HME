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
  
  
  cat(sprintf("Log-Likelihood: %f", logLik(obj)))
  if (expertbool) {
    
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
    
    if (wide_hme) {
      pars <- obj[["gate.pars"]][["0"]]
      cat(sprintf("\n--------------------------\ngate-%s\n", "0"))
      for (g in seq_along(pars)) {
        p <- pars[[g]]
        se <- gate.errs[["0"]][[g]]
        zstat <- p / se
        tbl <- cbind(p, se, zstat, 2 * pnorm(-abs(zstat)))
        dimnames(tbl) <- list(obj[['gate.pars.nms']],
                              c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)"))
        printCoefmat(tbl, digits = 3)
        cat("\n")
      }
    } else {
      pars <- obj[["gate.pars"]]
      for (g in names(pars)) {
        cat(sprintf("\n--------------------------\ngate-%s\n", g))
        p <- pars[[g]]
        se <- gate.errs[[g]]
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
}
