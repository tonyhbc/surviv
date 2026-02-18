#' Summarizing Cox Regression using Instrumental Variable
#'
#' Provide summary for a `surviv` object.
#'
#' @param object an object of the class "`surviv`", usually a call from "`coxiv_mom`".
#'
#' @examples
#' summary(coxiv_mom(Surv(time, censor) ~ x, data = dat))
#'
#' @export

summary.surviv = function(object, ...){
  coef_table = object$Result.Matrix
  colnames(coef_table) = c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")
  HR_table = cbind(exp(coef_table[,1]), exp(coef_table[,1] - qnorm(0.975)*coef_table[,2]),
                   exp(coef_table[,1] + qnorm(0.975)*coef_table[,2]))
  colnames(HR_table) = c("exp(coef)", "lower 95", "upper 95")
  summary_res = list(coefficients = coef_table, HR = HR_table)
  class(summary_res) <- "summary.surviv"
  return(summary_res)
}
