#' Get arguments for MakeADfun
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param family a probability distribution used as an error distribution
#' @param D distance matrix
#' @param W spatial weight matrix
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#' @export

getArg <- function(formula,
                   data,
                   family,
                   spatial = NULL,
                   D = NULL,
                   W = NULL) {

  # create base arguments for MakeADFun -------------------------------------

  # boolean to classify models
  re_term <- !is.null(lme4::findbars(formula))

  if (!re_term) {
    # class glm (no random effect)
    fr <- stats::model.frame(formula, data)
    trm <- terms(fr)

    y <- fr[, attr(trm, "response")]
    X <- stats::model.matrix(formula, data)
    Z <- as(matrix(0, nrow(fr), 1L), "TsparseMatrix")
    term <- list(list(n_fix = 1L,
                      n_group = 1L,
                      n_param = 1L))

    dataArg <- list(y = y,
                    X = X,
                    Z = Z,
                    term = term)
    parArg <- list(b = rep(1, ncol(X)),
                   log_sigma = log(1),
                   log_theta = log(1),
                   v = rep(0, ncol(Z)))
    reArg <- rlist <- NULL
    mapArg <- list(log_theta = factor(NA),
                   v = factor(rep(NA, ncol(Z))))

  } else {
    # class glmm (with non-spatial random effect)
    rlist <- lme4::lFormula(formula, data)
    fr <- rlist$fr
    trm <- terms(fr)

    y <- fr[, attr(trm, "response")]
    X <- rlist$X

    reTrms <- rlist$reTrms

    n_fix <- as.vector(sapply(reTrms$cnms, length))
    n_group <- as.vector(reTrms$nl)
    n_param <- as.vector(sapply(n_fix, function(x) choose(x, 2L) + x))

    term <- lapply(1L:length(reTrms$flist), function(i) {
      list(n_fix = n_fix[i],
           n_group = n_group[i],
           n_param = n_param[i])
    })

    # initial values for var-cov matrices
    init_log_theta <- unlist(sapply(1L:length(reTrms$flist), function(i) {
      c(rep(log(1), n_fix[i]), # diagonal elements
        rep(0, choose(n_fix[i], 2L))) # off-diagonal elements
    }))

    Zlist <- lapply(reTrms$Ztlist,
                    function(x) t(as.matrix(x)))
    Z <- as(do.call(cbind, Zlist), "TsparseMatrix")

    dataArg <- list(y = y,
                    X = X,
                    Z = Z,
                    term = term)
    parArg <- list(b = rep(0, ncol(X)),
                   log_sigma = 0,
                   log_theta = init_log_theta,
                   v = rep(0, ncol(Z)))
    reArg <- c("v")
    mapArg <- list()
  } # ifelse


  # spatial random effect ---------------------------------------------------

  if(is.null(spatial)) {
    # non-spatial model
    dataArg$D <- matrix(0, nrow = 1L, ncol = 1L)
    dataArg$W <- matrix(1, nrow = 1L, ncol = 1L)
    parArg$u <- rep(0, nrow(fr))
    mapArg$u <- factor(rep(NA, nrow(fr)))
    mapArg$log_phi <- factor(NA)
    mapArg$log_lambda <- factor(NA)
  } else {
    # spatial model
    dataArg$D <- D

    if (is.null(W)) {
      dataArg$W <- matrix(1, nrow = nrow(D), ncol = ncol(D))
    } else {
      dataArg$W <- W
    }

    parArg$u <- rep(0, nrow(fr))
    parArg$log_sigma <- log(sqrt(.Machine$double.eps))
    reArg <- c(reArg, "u")
    mapArg$log_sigma <- factor(NA)
  }

  parArg$log_phi <- 0
  parArg$log_lambda <- 0

  # family specific arguments -----------------------------------------------

  fam <- family

  .valid_link <- c(identity = 0L,
                   log = 1L)

  .valid_family <- c(gaussian = 0L,
                     poisson = 1L)

  # for no-variance family
  if (fam$family %in% c("poisson")) mapArg$log_sigma <- factor(NA)

  # family and link code
  dataArg$link <- .valid_link[fam$link]
  dataArg$family <- .valid_family[fam$family]


  # offset term -------------------------------------------------------------

  if (!is.null(attr(trm, "offset"))) {
    # column index for offset terms
    offcol <- attr(trm, "offset")
    Xi <- fr[ , offcol]

    # sum across offset columns if more than one term
    if (length(offcol) > 1) xi <- rowSums(Xi) else xi <- Xi

    dataArg$xi <- xi
  } else {
    dataArg$xi <- rep(0, nrow(fr))
  }

  # return ------------------------------------------------------------------

  list_out <- list(dataArg = dataArg,
                   parArg = parArg,
                   reArg = reArg,
                   mapArg = mapArg,
                   fam = fam,
                   rlist = rlist,
                   spatial = spatial)

  return(list_out)
}


#' Fit TMB model
#' @param tmbArg An object from \code{getArg()} containing arguments for \code{TMB::MakeADFun()}
#' @param verbose Logical. If \code{TRUE}, print maximum gradient \code{mgc} components while fitting.
#' @param control Optional control arguments for \code{nlminb()}. Supply as \code{list()}.
#' @export

fitTMB <- function(tmbArg,
                   verbose = FALSE,
                   control = list()) {

  # apply MakeADFun
  obj <- with(tmbArg,
              TMB::MakeADFun(data = dataArg,
                             parameters = parArg,
                             map = mapArg,
                             profile = NULL,
                             random = reArg,
                             DLL = "snaTMB",
                             silent = !verbose))

  # fit TMB model
  fit <- with(obj,
              stats::nlminb(start = par,
                            obj = fn,
                            gr = gr,
                            control = control))

  attr(fit, which = "inits") <- obj$par
  attr(fit, which = "method") <- obj$method
  attr(fit, which = "report") <- obj$report()
  attr(fit, which = "map") <- with(tmbArg, names(mapArg))

  # sd report
  sdr <- TMB::sdreport(obj, getJointPrecision = TRUE)

  # return
  list_out <- list(fit = fit,
                   sdr = sdr,
                   tmbArg = tmbArg)

  class(list_out) <- "snaTMB"

  return(list_out)
}


#' Report estimation tables with SEs
#' @param x Object class snaTMB
#' @export

report <- function(x) {
  if (class(x) != "snaTMB") stop("'x' must be calss 'snaTMB'")

  fit <- x$fit
  sdr <- x$sdr
  tmbArg <- x$tmbArg
  rlist <- tmbArg$rlist
  p_table <- TMB::summary.sdreport(sdr, "report")
  rows <- rownames(p_table)

  # n observations
  nobs <- nrow(tmbArg$dataArg$X)

  # exclude dummy parameters
  pmap <- stringr::str_remove_all(names(tmbArg$mapArg), "^log_")
  p_table <- p_table[-which(rows %in% pmap),]

  np <- nrow(p_table)

  if (tmbArg$fam$family %in% "gaussian" &
      is.null(tmbArg$spatial)) np <- np + 1

  DF <- nobs - np

  # return
  list_out <- list(fixef = format_fixef(sdr, tmbArg, DF))
  if (!is.null(rlist)) list_out$ranef <- format_ranef(fit, sdr, tmbArg, full = TRUE)
  if (!is.null(tmbArg$spatial)) list_out$spatial <- format_spatial(sdr, tmbArg, full = TRUE)

  return(list_out)
}


#' Format fixed effects
#' @param sdr Output from \code{TMB::sdreport()}
#' @param tmbArg Output from \code{snaTMB::getArg()}
#' @param DF Degree of freedom for Student-t distribution
#' @export

format_fixef <- function(sdr, tmbArg, DF = NULL) {

  p_table <- TMB::summary.sdreport(sdr, "report")
  df_b <- data.frame(Terms = colnames(tmbArg$dataArg$X),
                     Estimate = p_table[1:ncol(tmbArg$dataArg$X), 1L],
                     S.E. = p_table[1:ncol(tmbArg$dataArg$X), 2L],
                     row.names = NULL)

  if (tmbArg$fam$family %in% "gaussian") {
    # report t-values if gaussian
    v <- df_b$`t-value` <- ifelse(df_b$S.E. > 0,
                                  df_b$Estimate / df_b$S.E.,
                                  NaN)

    df_b$`Pr(>|t|)` <- 2 * pt(abs(v), df = DF, lower.tail = FALSE)
  } else {
    # report z-values if non-gaussian
    v <- df_b$`z-value` <- ifelse(df_b$S.E. > 0,
                                  df_b$Estimate / df_b$S.E.,
                                  NaN)

    df_b$`Pr(>|z|)` <- 2 * pnorm(abs(v), lower.tail = FALSE)
  }

  return(df_b)
}


#' Format random effects
#' @param fit Output from \code{\link{snaTMB::snaTMB}}
#' @param sdr Output from \code{\link{TMB::sdreport}}
#' @param tmbArg Output from \code{\link{snaTMB::getArg}}
#' @param full Logical. If TRUE, standard errors are printed.
#' @export

format_ranef <- function(fit, sdr, tmbArg, full = FALSE) {
  # v_theta: random effect sd, vector
  # R: correlation matrix
  p_table <- TMB::summary.sdreport(sdr, "report")
  v_theta <- p_table[which(rownames(p_table) == "theta"), 1L]
  v_theta_sd <- p_table[which(rownames(p_table) == "theta"), 2L]
  R <- attr(fit, "report")$cor
  R <- lapply(R, function(X) {
    X[] <- sprintf("%.2f", X)
    X[upper.tri(X, diag = TRUE)] <- ""
    return(X)
  })
  v_Rdim <- unlist(lapply(R, function(x) max(dim(x))))

  # list_fix: fixed effect terms which random effects apply
  # nr: # of random effects
  # v_np: # of parameters for each random effect
  # v_nf: # of fixed terms for each random effect
  # v_recol: random effect columns
  # v_f: random name assignment to theta vector
  # list_theta: named list for random sds
  list_fix <- tmbArg$rlist$reTrms$cnms
  v_recol <- names(list_fix)
  nr <- length(v_recol)
  v_np <- sapply(tmbArg$dataArg$term, function(x) x$n_param)
  v_nf <- sapply(tmbArg$dataArg$term, function(x) x$n_fix)
  v_f <- unlist(sapply(1:length(v_recol), function(i) rep(v_recol[i], v_np[i])))
  list_theta <- split(v_theta, v_f)
  list_theta_sd <- split(v_theta_sd, v_f)

  list_table <- lapply(1:nr, function(i) {
    d0 <- data.frame(Groups = rep(v_recol[i], v_nf[i]),
                     Terms = list_fix[[i]],
                     Std.Dev. = head(list_theta[[i]], n = v_nf[i]),
                     row.names = NULL)

    if (full) d0$S.E. <- head(list_theta_sd[[i]], n = v_nf[i])

    if (any(v_Rdim > 0)) {
      if (v_Rdim[i] == 0) {
        empty <- rep("Cor", max(v_Rdim))
        d0[, empty] <- ""
      } else {
        d0 <- cbind(d0, R[[i]])
      }
      names(d0) <- c("Groups",
                     "Terms",
                     "Std.Dev",
                     if (full) "S.E.",
                     "Cor",
                     rep("", max(v_Rdim) - 1))
    }

    return(d0)
  })

  df_re <- do.call(rbind, list_table)

  df_re$Groups[duplicated(df_re$Groups)] <- ""

  return(df_re)
}


#' Format spatial effects
#' @param sdr Output from \code{TMB::sdreport()}
#' @param tmbArg Output from \code{snaTMB::getArg()}
#' @param full Logical. Full report with SEs if \code{TRUE}.
#' @export

format_spatial <- function(sdr, tmbArg, full = FALSE) {

  p_table <- TMB::summary.sdreport(sdr, "report")
  rowid <- which(rownames(p_table) %in% c("lambda", "phi"))

  if (full) {
    df_spatial <- data.frame(Parameter = c("Spatial rate parameter",
                                           "Non-spatial standard deviation"),
                             Estimate = p_table[rowid, 1L],
                             S.E. = p_table[rowid, 2L],
                             row.names = NULL)
  } else {
    df_spatial <- data.frame(Parameter = c("Spatial rate parameter",
                                           "Non-spatial standard deviation"),
                             Estimate = p_table[rowid, 1L],
                             row.names = NULL)
  }

  return(df_spatial)
}
