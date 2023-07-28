#' Get arguments for MakeADfun
#'
#' @inheritParams snglmm
#' @export

get_arg <- function(formula,
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
    trm <- stats::terms(fr)

    y <- fr[, attr(trm, "response")]
    X <- stats::model.matrix(formula, data)
    Z <- methods::as(matrix(0, nrow(fr), 1L), "TsparseMatrix")
    term <- list(list(n_fix = 1L,
                      n_group = 1L,
                      n_param = 1L))

    data_arg <- list(y = y,
                     X = X,
                     Z = Z,
                     term = term)
    par_arg <- list(b = rep(1, ncol(X)),
                    log_sigma = log(1),
                    log_theta = log(1),
                    v = rep(0, ncol(Z)))
    re_arg <- rlist <- NULL
    map_arg <- list(log_theta = factor(NA),
                    v = factor(rep(NA, ncol(Z))))

  } else {
    # class glmm (with non-spatial random effect)
    rlist <- lme4::lFormula(formula, data)
    fr <- rlist$fr
    trm <- stats::terms(fr)

    y <- fr[, attr(trm, "response")]
    X <- rlist$X

    reTrms <- rlist$reTrms

    n_fix <- as.vector(sapply(reTrms$cnms, length))
    n_group <- as.vector(reTrms$nl)
    n_param <- as.vector(sapply(n_fix, function(x) choose(x, 2L) + x))

    term <- lapply(seq_len(length(reTrms$flist)),
                   function(i) {
                     list(n_fix = n_fix[i],
                          n_group = n_group[i],
                          n_param = n_param[i])
                   })

    # initial values for var-cov matrices
    init_log_theta <- unlist(sapply(seq_len(length(reTrms$flist)),
                                    function(i) {
                                      # 1st-L: diagonal elements
                                      # 2nd-L: off-diagonal elements
                                      c(rep(log(1), n_fix[i]),
                                        rep(0, choose(n_fix[i], 2L)))
                                    }))

    Zlist <- lapply(reTrms$Ztlist,
                    function(x) t(as.matrix(x)))
    Z <- methods::as(do.call(cbind, Zlist), "TsparseMatrix")

    data_arg <- list(y = y,
                     X = X,
                     Z = Z,
                     term = term)
    par_arg <- list(b = rep(0, ncol(X)),
                    log_sigma = 0,
                    log_theta = init_log_theta,
                    v = rep(0, ncol(Z)))
    re_arg <- c("v")
    map_arg <- list()
  } # ifelse


  # spatial random effect ---------------------------------------------------

  if (is.null(spatial)) {
    # non-spatial model
    data_arg$D <- matrix(0, nrow = 1L, ncol = 1L)
    data_arg$W <- matrix(1, nrow = 1L, ncol = 1L)
    par_arg$u <- rep(0, nrow(fr))
    map_arg$u <- factor(rep(NA, nrow(fr)))
    map_arg$log_phi <- factor(NA)
    map_arg$log_lambda <- factor(NA)
  } else {
    # spatial model
    data_arg$D <- D

    if (is.null(W)) {
      data_arg$W <- matrix(1, nrow = nrow(D), ncol = ncol(D))
    } else {
      data_arg$W <- W
    }

    par_arg$u <- rep(0, nrow(fr))
    par_arg$log_sigma <- log(sqrt(.Machine$double.eps))
    re_arg <- c(re_arg, "u")
    map_arg$log_sigma <- factor(NA)
  }

  par_arg$log_phi <- 0
  par_arg$log_lambda <- 0

  # family specific arguments -----------------------------------------------

  fam <- family

  .valid_link <- c(identity = 0L,
                   log = 1L)

  .valid_family <- c(gaussian = 0L,
                     poisson = 1L)

  # for no-variance family
  if (fam$family %in% c("poisson")) map_arg$log_sigma <- factor(NA)

  # family and link code
  data_arg$link <- .valid_link[fam$link]
  data_arg$family <- .valid_family[fam$family]

  # offset term -------------------------------------------------------------

  if (!is.null(attr(trm, "offset"))) {
    # column index for offset terms
    offcol <- attr(trm, "offset")
    Xi <- fr[, offcol]

    # sum across offset columns if more than one term
    if (length(offcol) > 1) xi <- rowSums(Xi) else xi <- Xi

    data_arg$xi <- xi
  } else {
    data_arg$xi <- rep(0, nrow(fr))
  }

  # return ------------------------------------------------------------------

  list_out <- list(data_arg = data_arg,
                   par_arg = par_arg,
                   re_arg = re_arg,
                   map_arg = map_arg,
                   fam = fam,
                   rlist = rlist,
                   spatial = spatial)

  return(list_out)
}


#' Fit TMB model
#'
#' @param tmb_arg An object from \code{\link{get_arg}} containing arguments for \code{TMB::MakeADFun}
#' @param verbose Logical. If \code{TRUE}, print maximum gradient \code{mgc} components while fitting.
#' @param control Optional control arguments for \code{\link{nlminb}}. Supply as \code{list()}.
#' @export

fitTMB <- function(tmb_arg,
                   verbose = FALSE,
                   control = list()) {

  # apply MakeADFun
  obj <- with(tmb_arg,
              MakeADFun(data = data_arg,
                        parameters = par_arg,
                        map = map_arg,
                        profile = NULL,
                        random = re_arg,
                        DLL = "snglmm",
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
  attr(fit, which = "map") <- with(tmb_arg, names(map_arg))

  # sd report
  sdr <- sdreport(obj, getJointPrecision = TRUE)

  # return
  list_out <- list(fit = fit,
                   sdr = sdr,
                   tmb_arg = tmb_arg)

  class(list_out) <- "snglmm"

  return(list_out)
}


#' Report estimation tables with SEs
#'
#' @param x Object class snglmm
#' @export

report <- function(x) {
  if (!inherits(x, "snglmm")) stop("'x' must be calss 'snglmm'")

  fit <- x$fit
  sdr <- x$sdr
  tmb_arg <- x$tmb_arg
  rlist <- tmb_arg$rlist
  p_table <- summary.sdreport(sdr, "report")
  rows <- rownames(p_table)

  # n observations
  nobs <- nrow(tmb_arg$data_arg$X)

  # exclude dummy parameters
  pmap <- stringr::str_remove_all(names(tmb_arg$map_arg), "^log_")
  p_table <- p_table[-which(rows %in% pmap), ]

  np <- nrow(p_table)

  if (tmb_arg$fam$family %in% "gaussian" && is.null(tmb_arg$spatial))
    np <- np + 1

  dfr <- nobs - np

  # return
  list_out <- list(fixef = format_fixef(sdr, tmb_arg))
  if (!is.null(rlist)) list_out$ranef <- format_ranef(fit, sdr, tmb_arg, full = TRUE)
  if (!is.null(tmb_arg$spatial)) list_out$spatial <- format_spatial(sdr, tmb_arg, full = TRUE)

  return(list_out)
}


#' Format fixed effects
#'
#' @inheritParams format_ranef
#' @export

format_fixef <- function(sdr, tmb_arg) {

  p_table <- summary.sdreport(sdr, "report")
  df_b <- data.frame(Terms = colnames(tmb_arg$data_arg$X),
                     Estimate = p_table[seq_len(ncol(tmb_arg$data_arg$X)), 1L],
                     S.E. = p_table[seq_len(ncol(tmb_arg$data_arg$X)), 2L],
                     row.names = NULL)

  v <- df_b$`z-value` <- ifelse(df_b$S.E. > 0,
                                df_b$Estimate / df_b$S.E.,
                                NaN)

  df_b$`Pr(>|z|)` <- 2 * stats::pnorm(abs(v), lower.tail = FALSE)

  return(df_b)
}


#' Format random effects
#'
#' @param fit Output from \code{snglmm::snglmm}
#' @param sdr Output from \code{TMB::sdreport}
#' @param tmb_arg Output from \code{snglmm::get_arg}
#' @param full Logical. If TRUE, standard errors are printed.
#' @export

format_ranef <- function(fit, sdr, tmb_arg, full = FALSE) {
  # v_theta: random effect sd, vector
  # R: correlation matrix
  p_table <- summary.sdreport(sdr, "report")
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
  list_fix <- tmb_arg$rlist$reTrms$cnms
  v_recol <- names(list_fix)
  nr <- length(v_recol)
  v_np <- sapply(tmb_arg$data_arg$term, function(x) x$n_param)
  v_nf <- sapply(tmb_arg$data_arg$term, function(x) x$n_fix)
  v_f <- unlist(sapply(seq_len(length(v_recol)),
                       function(i) rep(v_recol[i], v_np[i])))
  list_theta <- split(v_theta, v_f)
  list_theta_sd <- split(v_theta_sd, v_f)

  list_table <- lapply(seq_len(nr), function(i) {
    d0 <- data.frame(Groups = rep(v_recol[i], v_nf[i]),
                     Terms = list_fix[[i]],
                     Std.Dev. = utils::head(list_theta[[i]], n = v_nf[i]),
                     row.names = NULL)

    if (full) d0$S.E. <- utils::head(list_theta_sd[[i]], n = v_nf[i])

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
#'
#' @inheritParams format_ranef
#' @export

format_spatial <- function(sdr, tmb_arg, full = FALSE) {

  p_table <- summary.sdreport(sdr, "report")
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
