#' Spatial Network GLMM with TMB
#' @aliases NULL snglmm-package
#' @useDynLib snglmm
#' @import RcppEigen TMB
"_PACKAGE"


#' Fit a spatial model with TMB
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).
#' @param family Family assumed for an error distribution. Default \code{gaussian()}.
#' @param spatial Model type used for a spatial model (either \code{"exp"} or \code{"gaussian"}). Disabled by default.
#' @param D An optional distance matrix for a spatial model.
#' @param W An optional spatial weight matrix for a spatial model.
#' @param verbose Logical. If \code{TRUE}, print maximum gradient \code{mgc} components while fitting.
#' @param inits A list of initial parameter values.
#' @param control Optional arguments passed to \code{\link{nlminb}}.
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#' @export

snglmm <- function(formula,
                   data,
                   family = stats::gaussian(),
                   spatial,
                   D,
                   W,
                   verbose = FALSE,
                   inits,
                   control = list()) {


  # argument check ----------------------------------------------------------

  # call
  cl <- match.call()
  tmb_arg <- match.call(expand.dots = FALSE)

  # family check
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (is.null(family)) stop("'family' not recognized")

  # distance and sptial weight matrix
  if (!missing(W) && missing(D)) stop("Argument 'D' (distance matrix) is required to build a spatial model")

  if (!missing(D)) {
    if (!is.matrix(D)) {
      stop("'D' must be a matrix")
    } else {
      if (ncol(D) != nrow(D)) stop("Invalid dimensions in 'D'")
      if (!missing(W)) {
        if (!identical(dim(D), dim(W))) stop("Dimensions in 'D' and 'W' do not match")
      }
    }
  }

  # get arguments for fitTMB ------------------------------------------------

  index <- match(x = c("formula", "data", "spatial", "D", "W"),
                 table = names(tmb_arg),
                 nomatch = 0L)
  tmb_arg <- tmb_arg[c(1L, index)]

  tmb_arg$family <- family

  if (!missing(spatial)) {
    if (missing(D)) stop("Argument 'D' (distance matrix) is required to build a spatial model")

    if (!(spatial %in% c("exp", "gaussian"))) stop(paste("Method",
                                                         sQuote(spatial),
                                                         "is not implemented."))
  } else {
    if (!missing(D)) {
      message("Note: distance matrix D is supplied - a spatial exponential model is assumed by default. Specify 'spatial' argument to choose different types of models")
      tmb_arg$spatial <- "exp"
    }
  }

  tmb_arg[[1L]] <- quote(get_arg)
  tmb_arg <- eval(tmb_arg, envir = parent.frame())


  # start values ------------------------------------------------------------

  if (!missing(inits)) {
    par_arg <- tmb_arg$par_arg
    parnames <- names(par_arg)
    parsmap <- names(tmb_arg$map_arg)
    initspars <- parnames[!parnames %in% parsmap]

    bl_par <- !(names(inits) %in% parnames)
    bl_map <- names(inits) %in% parsmap
    n_elm_par <- sapply(par_arg, length)
    n_elm_inits <- sapply(inits, length)

    if (!is.list(inits)) {
      # stop if inits is not a list
      stop("'inits' must be supplied as a named list")
    } else {
      # check if inits values are substitutable
      if (any(bl_par)) {
        stop(paste(sQuote(names(inits)[bl_par]),
                   "not allowed to set initial values.",
                   "The following parameters are substitutable: ",
                   paste(sQuote(initspars), collapse = ",")))
      }

      if (any(bl_map)) {
        stop(paste(sQuote(names(inits)[bl_map]),
                   "not allowed to set initial values.",
                   "The following parameters are substitutable: ",
                   paste(sQuote(initspars), collapse = ",")))
      } else {

        # boolean: match in parameter vector lengths
        bl_parl <- n_elm_par[parnames %in% names(inits)] != n_elm_inits

        if (any(bl_parl)) {
          mismatch <- names(inits)[bl_parl]
          stop(paste("Vector lengths in the following parameters do not match:",
                     paste(sQuote(mismatch), collapse = ","),
                     "\n",
                     "The expected lengths are:\n",
                     paste0("length(", mismatch, ")",
                            " = ",
                            n_elm_par[names(n_elm_par) %in% mismatch],
                            collapse = ", ")))
        } else {
          index <- intersect(parnames, names(inits))
          tmb_arg$par_arg <- utils::modifyList(x = par_arg,
                                               val = inits[index])
        }

      }
    }
  }

  # distribution compatibility ----------------------------------------------

  if (grepl("nbinom|pois", family$family)) {
    eta <- tmb_arg$data_arg$y
    if (any(abs(eta - round(eta)) > 0.001)) {
      stop(sprintf("non-integer counts in a %s model",
                   family$family))
    }
  }

  # fit TMB model -----------------------------------------------------------

  # fitTMB passes arguments to TMB::MakeADFun() and stats::nlminb()
  opt <- fitTMB(tmb_arg, verbose, control)
  opt$call <- cl

  opt
}


#' Perform Regression/Universal Kriging
#'
#' @param object Object class \code{snglmm}
#' @param newdata A required data frame containing predictor values at new sites.
#' @param cD A required cross-distance matrix for a spatial model.
#' @param cW An optional cross-weight matrix for a spatial model.
#' @export

kriging <- function(object,
                    newdata,
                    cD,
                    cW) {

  # input check -------------------------------------------------------------

  if (!inherits(object, "snglmm"))
    stop("object class must be 'snglmm'")

  if (is.null(object$tmb_arg$spatial))
    stop("the model contains no spatial term;
         kriging methods are not applicable.")

  if (missing(cD))
    stop("A cross-distance matrix is required for regression/universal kriging.")

  if (!is.matrix(cD))
    stop("Supply 'cD' as a matrix.")

  # y: n observed values
  # X: n x p predictor matrix at observed sites
  # X0: m x p predictor matrix at new sites
  # D: n x n distance matrix
  # W: n x n weight matrix
  y <- object$tmb_arg$data_arg$y
  X <- object$tmb_arg$data_arg$X
  X0 <- newdata
  D <- object$tmb_arg$data_arg$D
  W <- object$tmb_arg$data_arg$W

  sdr <- summary(object$sdr, "report")
  rid <- rownames(sdr)
  beta <- sdr[rid == "b", 1L]
  phi <- sdr[rid == "phi", 1L]
  lambda <- sdr[rid == "lambda", 1L]

  # dimension check cD cW
  ## row
  vrows <- c(`cD` = nrow(cD),
             `cW` = if (!missing(cW)) nrow(cW))

  if (any(vrows != nrow(X)))
    stop(paste("Dimension mismatch.",
               sQuote(names(vrows)[vrows != nrow(X)]),
               "must have",
               nrow(X),
               "rows (the number of observations in the original model)."))

  ## column
  vcols <- c(`cD` = ncol(cD),
             `cW` = if (!missing(cW)) ncol(cW))

  if (any(vcols != nrow(X0)))
    stop(paste("Dimension mismatch.",
               sQuote(names(vcols)[vcols != nrow(X0)]),
               "must have",
               nrow(X0),
               "columns (the number of prediction sites)."))


  # SIGMA: n x n vcov matrix at obseved sites
  # SIGMA0: n x m vcov matrix of at new sites
  # TAU: n x n precision matrix at observed sites
  SIGMA <- W * (phi^2 * exp(-lambda * D))
  TAU <- Matrix::solve(SIGMA)

  if (missing(cW)) {
    SIGMA0 <- phi^2 * exp(-lambda * cD)
  } else {
    if (!is.matrix(cW))
      stop("Supply 'cW' as a matrix.")

    SIGMA0 <- cW * phi^2 * exp(-lambda * cD)
  }

  # u: n vector of residuals at observed sites
  u <- attr(object$fit, "report")$u

  # nu: m x 1 weighted residuals
  # z: predicted values at m new sites
  nu <- drop(t(SIGMA0) %*% TAU %*% u)
  z <- drop(X0 %*% beta) + nu

  return(z)
}
