
#' Extract Fixed Effects
#' @description This function is generic; method functions can be written to handle specific classes of objects.
#' @inheritParams nlme::fixef
#' @importFrom nlme fixef
#' @aliases fixef fixef.snaTMB
#' @export fixef
#' @export

fixef.snaTMB <- function(object, ...) {
  xlabels <- colnames(object$tmb_arg$data_arg$X)
  b_table <- summary.sdreport(object$sdr, "report")[, 1L]
  b_table <- b_table[which(names(b_table) == "b")]
  names(b_table) <- xlabels
  return(b_table)
}

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#' @description Returns the variance-covariance matrix of the main parameters of a fitted model object. The “main” parameters of model correspond to those returned by \code{\link{coef}}, and typically do not contain a nuisance scale parameter (\code{\link{sigma}}).
#' @inheritParams stats::vcov
#' @aliases vcov vcov.snaTMB
#' @export

vcov.snaTMB <- function(object, ...) {
  xlabels <- colnames(object$tmb_arg$data_arg$X)
  m <- object$sdr$cov.fixed
  b_index <- which(colnames(m) == "b")
  m <- data.matrix(m[b_index, b_index])
  if (sum(dim(m)) != 0) rownames(m) <- colnames(m) <- xlabels
  return(m)
}


#' Print Values
#' @description \code{print} prints its argument and returns it invisibly (via \print{\link{invisible}}(x)). It is a generic function which means that new printing methods can be easily added for new \code{\link{class}}es.
#' @inheritParams base::print
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
  } else {
    cat("No coefficients\n")
  }
  cat("\n")

  # random effects
  if (!is.null(x$tmb_arg$rlist)) {
    cat("Random effects:\n")
    print.data.frame(snaTMB::format_ranef(x$fit, x$sdr, x$tmb_arg),
                     digits = digits,
                     row.names = FALSE)
    cat("\n")
  }

  # spatial random effects
  if (!is.null(x$tmb_arg$spatial)) {
    cat("Spatial random effects:\n")
    print.data.frame(snaTMB::format_spatial(x$sdr, x$tmb_arg),
                     digits = digits,
                     row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}

#' Object Summaries
#' @description \code{summary} is a generic function used to produce result summaries of the results of various model fitting functions. The function invokes particular \code{\link{methods}} which depend on the \code{\link{class}} of the first argument.
#' @inheritParams base::summary
#' @aliases summary summary.snaTMB
#' @export

summary.snaTMB <- function(object,
                           digits = max(3L, getOption("digits") - 3L),
                           ...) {

  # model information -------------------------------------------------------
  fit <- object$fit
  sdr <- object$sdr
  p_table <- summary.sdreport(sdr, "report")
  tmb_arg <- object$tmb_arg
  rlist <- object$tmb_arg$rlist

  # n observations
  nobs <- nrow(tmb_arg$data_arg$X)

  # n parameters
  fixef_np <- ncol(tmb_arg$data_arg$X)
  ranef_np <- ifelse(is.null(rlist),
                     yes = 0,
                     no = sum(sapply(tmb_arg$data_arg$term, function(x) x$n_param)))
  spatial_np <- ifelse(is.null(tmb_arg$spatial),
                       yes = 0,
                       no = 2) # lambda and phi

  np <- fixef_np + ranef_np + spatial_np

  if (tmb_arg$fam$family %in% "gaussian" && is.null(tmb_arg$spatial))
    np <- np + 1

  # degrees of freedom
  dfr <- ifelse(nobs > np, nobs - np, NaN)

  # fixed effects -----------------------------------------------------------

  # format_fixef(): utility function
  df_b <- snaTMB::format_fixef(sdr, tmb_arg)

  # random effects ----------------------------------------------------------

  # format_ranef(): utility function
  if (!is.null(rlist)) df_re <- snaTMB::format_ranef(fit, sdr, tmb_arg)

  # spatial effects ---------------------------------------------------------

  if (!is.null(tmb_arg$spatial)) df_spatial <- snaTMB::format_spatial(sdr, tmb_arg)

  # likelihood + others -----------------------------------------------------

  loglik <- -as.numeric(object$fit$objective)
  dev <- -2 * loglik
  aic <- dev + 2 * np
  bic <- dev + np * log(nobs)
  info_table <- c(`AIC` = aic,
                  `BIC` = bic,
                  `logLik` = loglik,
                  `Deviance` = dev,
                  `df.resid` = dfr)

  # print -------------------------------------------------------------------

  # information criteria
  print.default(info_table, digits = digits, ...)
  cat("\n")

  # fixed effects
  cat("Fixed effects: \n")
  if (sum(dim(tmb_arg$data_arg$X)) > 0) {
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
    cat("Random effects: \n")
    print.data.frame(df_re,
                     digits = digits,
                     row.names = FALSE,
                     na.print = "",
                     ...)
    cat("\n")
  }

  # spatial random effects
  if (!is.null(tmb_arg$spatial)) {
    cat(paste0("Spatial random effects: \n"))
    print.data.frame(df_spatial,
                     digits = digits,
                     row.names = FALSE,
                     ...)
    cat("\n")
  }

  # family-specific dispersion
  if (tmb_arg$fam$family %in% "gaussian") {
    cat(paste0("Family specific dispersion (", tmb_arg$fam$family, "): ",
               format(p_table[which(rownames(p_table) == "sigma"), "Estimate"],
                      digits = digits)))
    cat("\n")
    if (!is.null(tmb_arg$spatial)) cat("Note: family-specific dispersion parameter is fixed at zero to avoid overparameterization (non-spatial standard deviation represents a comaptible dispersion parameter)")
  }

}
