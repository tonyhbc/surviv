#' Print the output from a `coxivomom` model fit.
#'
#' Provide succinct and organized output for a `coxivomom` model, such as `coxiv_omom()`.
#'
#' @param x an object of the class "`coxivomom`", usually created by `coxiv_omom()`.
#' @param digits the number of significant digits to display for coefficients and hazard ratios.
#' @param ... additional arguments (currently unused).
#'
#' @examples
#' print(coxiv_omom(Surv(time, censor) ~ x, data = dat))
#'
#' @export

print.coxivomom <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  # Verify the class of the input object coefficients
  if (!inherits(x, "coxivomom")) {
    stop("Input must be of class 'coxivomom'.")
  }

  # Header
  cat("Orthogonality method-of-moment IV estimation of Cox model\n")
  cat(strrep("=", 57), "\n")

  # Model Call
  cat("Call:\n")
  print(x$call)
  cat(strrep("-", 57), "\n")

  # Sample Size and Events
  cat(sprintf("Sample size: %d\n", x$n))
  cat(sprintf("Number of events: %d\n", x$events))
  cat(strrep("-", 57), "\n")

  # Treatment effect
  cat(sprintf("Marginal treatment effect (%s):\n", x$input$trt))
  cat(sprintf("- Estimate: %.3f [%.3f, %.3f] \n", x$est_coef$Estimate[1], x$conf_int$`2.5%`[1], x$conf_int$`97.5%`[1]))
  cat(sprintf("- Hazard ratio: %.3f [%.3f, %.3f] \n", exp(x$est_coef$Estimate[1]),
              exp(x$conf_int$`2.5%`[1]), exp(x$conf_int$`97.5%`[1])))
  cat(strrep("-", 57), "\n")

  # Coefficients and Hazard Ratios
  est <- x$est_coef$Estimate
  se <- x$est_coef$Std.Error
  HR <- exp(est)
  coef_summary <-
    round(as.matrix(data.frame(
      est = est,
      se = se,
      HR = HR
  )), digits = 3)
  rownames(coef_summary) <- rownames(x$est_coef)

  cat("Coefficients:\n")
  print(coef_summary)
  cat(strrep("-", 57), "\n")
  cat("Nested exclusion test for weak IV:\n")
  cat(sprintf("- F-statistics: %.2f%s\n", x$iv_diag$f_stat, ifelse(x$iv_diag$f_stat > 10, " [> 10.0]", " [< 10.0]")))
  cat(strrep("-", 57), "\n")

  # Footer
  cat("A `coxivomom` object\n")
  invisible(x)
}
