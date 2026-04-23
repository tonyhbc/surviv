#' Print method for coxivseq objects
#'
#' @param x An object of class `coxivseq` returned by [coxiv_seq()].
#' @param digits Number of decimal places for printed coefficient summaries.
#' @param ... Ignored.
#'
#' @return The input object `x`, invisibly. Called for its side effect of printing a summary to the console.
#'
#' @export
print.coxivseq <- function(x, digits = 4, ...) {
  if (!inherits(x, "coxivseq")) {
    return(NextMethod())
  }

  title <- "Sequential 2SRI Cox analysis via sequential trials emulation"
  cat(title, "\n", sep = "")
  cat(strrep("=", nchar(title)), "\n", sep = "")

  cat(sprintf(" * Time-varying treatment: %s\n", x$tvtrt))
  cat(sprintf(" * Baseline IV: %s\n", x$iv))
  cat(strrep("~", nchar(title)), "\n", sep = "")

  cat(sprintf("Composite analysis across %d emulated trials \n", x$n_trials))
  cat(sprintf("Stratified by trial: %s\n", if (isTRUE(x$stratify)) "yes" else "no"))
  cat(sprintf("Number of individuals: %d\n", x$n_ids))
  cat(sprintf("Number of rows (composite dataset): %d\n", x$n))
  cat(sprintf("Number of events: %d\n", x$n_event))

  n_trtfit <- 0L
  if (!is.null(x$trtfit_by_trial)) {
    n_trtfit <- sum(!vapply(x$trtfit_by_trial, is.null, logical(1)))
  }
  n_fit <- 0L
  if (!is.null(x$fit_by_trial)) {
    n_fit <- sum(!vapply(x$fit_by_trial, is.null, logical(1)))
  }

  cat(sprintf("First-stage models: %d fit(s) available in [ $trtfit_by_trial ]\n", n_trtfit))
  cat(sprintf("Second-stage models: %d fit(s) available in [ $fit_by_trial ]\n", n_fit))
  cat(strrep("~", nchar(title)), "\n", sep = "")

  cat("Coefficients (composite fit):\n")
  tab <- x$coef[, c("logHR", "HR", "se", "p.val"), drop = FALSE]
  tab_chr <- as.matrix(tab)
  tab_chr[,] <- formatC(as.matrix(tab), format = "f", digits = digits)
  print(noquote(tab_chr), quote = FALSE, right = TRUE)

  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("A `coxivseq` object\n")
  invisible(x)
}
