#' Fit a spatial error model (SEM) with TMB
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#' @useDynLib semtmb
#' @export
#'

semtmb <- function(formula,
                   data) {

  y <- model.frame(formula, data)[,1]
  X <- model.matrix(formula, data)

  obj <- TMB::MakeADFun(data = list(y = y,
                                    D = D,
                                    X = X),
                        parameters = list(b = rep(0, ncol(X)),
                                          log_sigma = log(0.1),
                                          log_theta = log(10)),
                        DLL = "semtmb",
                        silent = TRUE)

  opt <- nlminb(start = obj$par,
                obj = obj$fn,
                gr = obj$gr)

  rep <- TMB::sdreport(obj)

  return(summary(rep, "report", p.value = TRUE))

}
