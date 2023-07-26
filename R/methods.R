
#' Get fixed effects
#' @param object Any fitted model object from which fixed effects estimates can be extracted.
#' @param \dots	Additional arguments. Currently none.
#' @importFrom nlme fixef
#' @aliases fixef fixef.snaTMB
#' @export fixef
#' @export

fixef.snaTMB <- function(object, ...) {
  xlabels <- colnames(object$tmbArg$dataArg$X)
  b_table <- TMB::summary.sdreport(object$sdr, "report")[, 1L]
  b_table <- b_table[which(names(b_table) == "b")]
  names(b_table) <- xlabels
  return(b_table)
}


#' Get variance-covariance for fixed effects
#' @param object A fitted model object, typically.
#' @param \dots Additional arguments. Currently none.
#' @aliases vcov vcov.snaTMB
#' @export

vcov.snaTMB <- function(object, ...) {
  xlabels <- colnames(object$tmbArg$dataArg$X)
  m <- object$sdr$cov.fixed
  b_index <- which(colnames(m) == "b")
  m <- data.matrix(m[b_index, b_index])
  if (sum(dim(m)) != 0) rownames(m) <- colnames(m) <- xlabels
  return(m)
}


#' Print snaTMB output
#' @param x Object class snaTMB
#' @aliases print print.snaTMB
#' @export

print.snaTMB <- function(x,
                         digits = max(3L, getOption("digits") - 3L),
                         ...) {

  # fixed effects
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(fixef(x))) {
    cat("Fixed effects:\n")
    print.default(fixef(x),
                  digits = digits,
                  print.gap = 2L,
                  quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")

  # random effects
  if (!is.null(x$tmbArg$rlist)) {
    cat("Random effects:\n")
    print.data.frame(snaTMB::format_ranef(x$fit, x$sdr, x$tmbArg),
                     digits = digits,
                     row.names = FALSE)
    cat("\n")
  }

  # spatial random effects
  if (!is.null(x$tmbArg$spatial)) {
    cat("Spatial random effects:\n")
    print.data.frame(snaTMB::format_spatial(x$sdr, x$tmbArg),
                     digits = digits,
                     row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}


#' Summarize snaTMB output
#' @param object An object for which a summary is desired.
#' @param \dots Additional arguments affecting the summary produced.
#' @aliases summary summary.snaTMB
#' @export

summary.snaTMB <- function(object,
                           digits = max(3L, getOption("digits") - 3L),
                           ...) {

  # model information -------------------------------------------------------
  fit <- object$fit
  sdr <- object$sdr
  p_table <- TMB::summary.sdreport(sdr, "report")
  tmbArg <- object$tmbArg
  rlist <- object$tmbArg$rlist

  # n observations
  nobs <- nrow(tmbArg$dataArg$X)

  # n parameters
  fixef_np <- ncol(tmbArg$dataArg$X)
  ranef_np <- ifelse(is.null(rlist),
                     yes = 0,
                     no = sum(sapply(tmbArg$dataArg$term, function(x) x$n_param)))
  spatial_np <- ifelse(is.null(tmbArg$spatial),
                       yes = 0,
                       no = 2) # lambda and phi

  np <- fixef_np + ranef_np + spatial_np

  if (tmbArg$fam$family %in% "gaussian" & is.null(tmbArg$spatial))
    np <- np + 1

  # degrees of freedom
  DF <- ifelse(nobs > np, nobs - np, NaN)

  # fixed effects -----------------------------------------------------------

  # format_fixef(): utility function
  df_b <- snaTMB::format_fixef(sdr, tmbArg, DF)

  # random effects ----------------------------------------------------------

  # format_ranef(): utility function
  if (!is.null(rlist)) df_re <- snaTMB::format_ranef(fit, sdr, tmbArg)

  # spatial effects ---------------------------------------------------------

  if (!is.null(tmbArg$spatial)) df_spatial <- snaTMB::format_spatial(sdr, tmbArg)

  # likelihood + others -----------------------------------------------------

  loglik <- -as.numeric(object$fit$objective)
  dev <- -2 * loglik
  aic <- dev + 2 * np
  bic <- dev + np * log(nobs)
  info_table <- c(`AIC` = aic,
                  `BIC` = bic,
                  `logLik` = loglik,
                  `Deviance` = dev,
                  `df.resid` = DF)

  # print -------------------------------------------------------------------

  # information criteria
  print.default(info_table, digits = digits, ...)
  cat("\n")

  # fixed effects
  cat("Fixed effects: \n");
  if (sum(dim(tmbArg$dataArg$X)) > 0) {
    print.data.frame(df_b,
                     digits = digits,
                     row.names = FALSE,
                     ...)
    cat("\n")
  } else {
    cat("No coefficients \n")
  }

  # non-spatial random effects
  if (!is.null(rlist)) {
    cat("Random effects: \n");
    print.data.frame(df_re,
                     digits = digits,
                     row.names = FALSE,
                     na.print = "",
                     ...)
    cat("\n");
  }

  # spatial random effects
  if (!is.null(tmbArg$spatial)) {
    cat(paste0("Spatial random effects: \n"))
    print.data.frame(df_spatial,
                     digits = digits,
                     row.names = FALSE,
                     ...)
    cat("\n")
  }

  # family-specific dispersion
  if (tmbArg$fam$family %in% "gaussian") {
    cat(paste0("Family specific dispersion (", tmbArg$fam$family, "): ",
               format(p_table[which(rownames(p_table) == "sigma"), "Estimate"],
                      digits = digits)))
    cat("\n")
    if (!is.null(tmbArg$spatial)) cat("Note: family-specific dispersion parameter is fixed at zero to avoid overparameterization (non-spatial standard deviation represents a comaptible dispersion parameter)")
  }

}


