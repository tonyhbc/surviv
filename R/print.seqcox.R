#' Print method for seqcox objects
#'
#' @param x An object of class `seqcox` returned by [seqcox()].
#' @param digits Number of decimal places for printed coefficient summaries.
#' @param ... Ignored.
#'
#' @return The input object `x`, invisibly. Called for its side effect of printing a summary to the console.
#'
#' @export
print.seqcox <- function(x, digits = 4, ...) {
  if (!inherits(x, "seqcox")) {
    return(NextMethod())
  }

  title <- "Sequential Cox analysis via sequential trials emulation"
  cat(title, "\n", sep = "")
  cat(strrep("=", nchar(title)), "\n", sep = "")

  cat(sprintf("Composite analysis across %d emulated trials \n", x$n_trials))
  cat(sprintf("Stratified by trial : %s\n", if (isTRUE(x$stratify)) "yes" else "no"))
  cat(sprintf("Number of individuals : %d\n", x$n_ids))
  cat(sprintf("Number of events: %d\n", x$n_event))
  cat(sprintf("SE type: %s\n", x$se_type))

  n_fit_by_trial <- 0L
  if (!is.null(x$fit_by_trial)) {
    n_fit_by_trial <- sum(!vapply(x$fit_by_trial, is.null, logical(1)))
  }
  cat(sprintf("Trial-specific fits: %d model(s) available in [ $fit_by_trial ]\n",
              n_fit_by_trial))

  cat(strrep("~", nchar(title)), "\n", sep = "")
  cat("Coefficients (composite fit):\n")

  tab <- x$coef[, c("logHR", "HR", "se", "p.val"), drop = FALSE]
  tab_chr <- as.matrix(tab)
  tab_chr[,] <- formatC(as.matrix(tab), format = "f", digits = digits)
  print(noquote(tab_chr), quote = FALSE, right = TRUE)

  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("A `seqcox` object\n")
  invisible(x)
}
