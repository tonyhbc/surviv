#' Print the output from a `coxivtsrif` model fit.
#'
#' Provide succinct and organized output for a `coxivtsrif` model, such as `coxiv_tsrif()`.
#'
#' @param x an object of the class "`coxivtsrif`", usually created by `coxiv_tsrif()`.
#' @param digits the number of significant digits to display for coefficients and hazard ratios.
#' @param ... additional arguments (currently unused).
#'
#' @examples
#' print(coxiv_tsrif(Surv(time, censor) ~ x, data = dat))
#'
#' @export

print.coxivtsrif <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  # Verify the class of the input object coefficients
  if (!inherits(x, "coxivtsrif")) {
    stop("Input must be of class 'coxivtsrif'.")
  }

  # Header
  cat("TSRI-Frailty instrumental variable analysis of Cox model \n")
  cat(strrep("=", 56), "\n")

  # Separate time-constant and time-varying effect case

  if (x$input$tvareff) {
    # Model Call
    cat("Call:\n")
    cat("- First-stage (treatment) model:\n")
    cat(paste(deparse(x$input$trtformula), collapse = " "), "\n")
    cat("- Second-stage (outcome) model:\n")
    cat(paste(deparse(x$input$formula), collapse = " "), "\n")
    cat(sprintf("* Time-varying treatment effect updates at t = %d\n", x$input$tchange))
    cat(strrep("-", 56), "\n")

    # Sample Size and Events
    cat(sprintf("Sample size: %d\n", x$surv_curve$n))
    cat(sprintf("Number of events: %d\n", sum(x$surv_curve$n.event)))
    cat(strrep("-", 56), "\n")

    # Treatment effect
    cat(sprintf("Conditional treatment effect (%s):\n", x$input$trt))
    cat(sprintf("# Period 1 # (t<%.2f):\n", x$input$tchange))
    cat(sprintf("  - Estimate: %.3f [%.3f, %.3f]\n", x$est_coef[1], x$conf_int[1,1], x$conf_int[1,2]))
    cat(sprintf("  - Hazard ratio: %.3f [%.3f, %.3f]\n", exp(x$est_coef[1]), exp(x$conf_int[1,1]), exp(x$conf_int[1,2])))
    cat(strrep("~", 22), "\n")
    cat(sprintf("# Period 2 # (t>%.2f):\n", x$input$tchange))
    cat(sprintf("  - Estimate: %.3f [%.3f, %.3f]\n", x$est_coef[2], x$conf_int[2,1], x$conf_int[2,2]))
    cat(sprintf("  - Hazard ratio: %.3f [%.3f, %.3f]\n", exp(x$est_coef[2]), exp(x$conf_int[2,1]), exp(x$conf_int[2,2])))
    cat(strrep("-", 56), "\n")

    # Coefficients estimates
    cat("Coefficients:\n")
    cat("# Period 1 #:\n")

    beta1 <- x$out_model1$coefficients
    se1   <- sqrt(diag(x$out_model1$var))

    # If you have a special variance for the treatment effect, overwrite SE for that row:
    if (!is.null(x$variance) && length(x$variance) >= 1 && x$input$trt %in% names(se1)) {
      se1[x$input$trt] <- sqrt(x$variance[1])
    }

    tab1  <- .coef_table(beta1, se1)
    tab1f <- .fmt_noquote(tab1, digits = digits)

    .print_text(tab1f, right = TRUE)

    cat(strrep("~", 13), "\n")
    cat("# Period 2 #:\n")

    beta2 <- setNames(x$est_coef[2], x$input$trt)
    se2   <- setNames(sqrt(x$variance[2]), x$input$trt)

    tab2  <- .coef_table(beta2, se2)
    tab2f <- .fmt_noquote(tab2, digits = digits)

    .print_text(tab2f, right = TRUE)
    cat(strrep("-", 56), "\n")
  } else {
    # Model Call
    cat("Call:\n")
    cat("- First-stage (treatment) model:\n")
    cat(paste(deparse(x$input$trtformula), collapse = " "), "\n")
    cat("- Second-stage (outcome) model:\n")
    cat(paste(deparse(x$input$formula), collapse = " "), "\n")
    cat("* Time-constant treatment effect\n")
    cat(strrep("-", 56), "\n")

    # Sample Size and Events
    cat(sprintf("Sample size: %d\n", x$surv_curve$n))
    cat(sprintf("Number of events: %d\n", sum(x$surv_curve$n.event)))
    cat(strrep("-", 56), "\n")

    # Treatment effect
    cat(sprintf("Conditional treatment effect (%s):\n", x$input$trt))
    cat(sprintf("- Estimate: %.3f [%.3f, %.3f] \n", x$est_coef, x$conf_int[1], x$conf_int[2]))
    cat(sprintf("- Hazard ratio: %.3f [%.3f, %.3f] \n", exp(x$est_coef), exp(x$conf_int[1]), exp(x$conf_int[2])))
    cat(strrep("-", 56), "\n")

    # Coefficients estimates
    cat("Coefficients:\n")
    beta <- x$out_model$coefficients
    se   <- sqrt(diag(x$out_model$var))

    tab  <- .coef_table(beta, se)
    tabf <- .fmt_noquote(tab, digits = digits)

    .print_text(tabf, right = TRUE)
    cat(strrep("-", 56), "\n")
  }

  cat("Nested exclusion test for weak IV:\n")
  cat(sprintf("- F-statistics: %.2f%s\n", x$iv_diag$f_stat, ifelse(x$iv_diag$f_stat > 10, " [> 10.0]", " [< 10.0]")))
  cat(strrep("-", 56), "\n")

  # Footer
  cat("A `coxivtsrif` object\n")
  invisible(x)

  # cat(strrep("-", 56), "\n")
  # cat("Nested IV-exclusion anova for weak IV:\n\n")
  # cat(sprintf("F-statistics: %.2f%s\n", x$iv_diag$f_stat, ifelse(x$iv_diag$f_stat > 10, " [> 10.0]", " [< 10.0]")))
  # cat(strrep("-", 56), "\n")

}

# helper: force console-style printing in Rmd + console
.print_text <- function(x, ...) {
  txt <- utils::capture.output(print(x, ...))
  cat(paste0(txt, collapse = "\n"), "\n")
}

# helper: format a numeric matrix to fixed decimals, but print WITHOUT quotes
.fmt_noquote <- function(mat, digits) {
  mat <- as.matrix(mat)  # ensure 2D
  mat_chr <- matrix(
    formatC(mat, format = "f", digits = digits),
    nrow = nrow(mat), ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  noquote(mat_chr)
}

# helper: build coef / se / HR numeric matrix
.coef_table <- function(beta, se) {
  rn <- names(beta)

  beta_num <- as.numeric(beta)
  se_num   <- as.numeric(se)

  # restore names after numeric coercion
  if (!is.null(rn)) {
    names(beta_num) <- rn
    names(se_num)   <- rn
  }

  out <- cbind(est = beta_num, se = se_num, HR = exp(beta_num))

  if (!is.null(rn)) rownames(out) <- rn
  out
}
