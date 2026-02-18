#' Print the output from a `surviv` model fit.
#'
#' This method will print a better-looking output for summary of a `surviv` model object.
#'
#' @param object a `summary.surviv` object.
#'
#' @examples
#' print(summary(CoxIVeq(Surv(time, censor) ~ x, data = dat)))
#'
#' @export

print.summary.surviv <- function(object, digits=max(1L, getOption("digits") - 3L), ...) {
  # model
  cat("Summary of Cox IV MOM model Result\n")
  cat("=============================\n")
  cat("Coefficients: \n")
  printCoefmat(object$coefficients)
  cat("_____________________________\n")
  cat("\n")
  cat("\n")
  cat("Marginal Hazard Ratios: \n")
  print(object$HR)
  cat("______________________________________________\n")
  invisible(x)
}
