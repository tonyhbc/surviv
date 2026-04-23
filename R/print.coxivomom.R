#' Print the output from a `coxivomom` model fit.
#'
#' Provide succinct and organized output for a `coxivomom` model, such as
#' `coxiv_omom()`.
#'
#' @param x an object of class `"coxivomom"`, usually created by
#'   `coxiv_omom()`.
#' @param digits the number of decimal places to display for coefficient
#'   summaries and hazard ratios.
#' @param ... additional arguments (currently unused).
#'
#' @return The input object `x`, invisibly. Called for its side effect of printing a summary to the console.
#'
#' @examples
#' data("VitD", package = "surviv")
#' fit_omom <- surviv::coxiv_omom(
#'   formula = survival::Surv(time, death) ~ vitd + age,
#'   trt = "vitd",
#'   iv = "filaggrin",
#'   data = VitD
#' )
#' print(fit_omom)
#'
#' @export
print.coxivomom <- function(x, digits = 3L, ...) {
  if (!inherits(x, "coxivomom")) {
    stop("Input must be of class 'coxivomom'.")
  }

  .print_text <- function(obj, ...) {
    txt <- utils::capture.output(print(obj, ...))
    cat(paste0(txt, collapse = "\n"), "\n")
  }

  .fmt_noquote <- function(mat, digits) {
    mat <- as.matrix(mat)
    out <- matrix(
      formatC(mat, format = "f", digits = digits),
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
    noquote(out)
  }

  .standardize_coef_df <- function(obj) {
    coef_df <- obj$coef
    if (is.null(coef_df)) {
      coef_df <- obj$est_coef
    }
    if (is.null(coef_df)) {
      stop("Cannot find coefficient summary in `x$coef` or `x$est_coef`.")
    }

    coef_df <- as.data.frame(coef_df, check.names = FALSE)

    nm <- names(coef_df)
    map <- c(
      "Estimate" = "logHR",
      "Std.Error" = "se",
      "Z.value" = "Z",
      "P(>|z|)" = "pval",
      "logHR" = "logHR",
      "se" = "se",
      "Z" = "Z",
      "pval" = "pval"
    )
    nm2 <- unname(map[nm])
    names(coef_df) <- ifelse(is.na(nm2), nm, nm2)

    needed <- c("logHR", "se")
    if (!all(needed %in% names(coef_df))) {
      stop("Coefficient summary must contain `logHR` and `se` columns.")
    }

    coef_df$logHR <- as.numeric(coef_df$logHR)
    coef_df$se <- as.numeric(coef_df$se)
    if ("Z" %in% names(coef_df)) coef_df$Z <- as.numeric(coef_df$Z)
    if ("pval" %in% names(coef_df)) coef_df$pval <- as.numeric(coef_df$pval)

    coef_df
  }

  .standardize_confint <- function(obj, rn = NULL) {
    ci <- obj$conf_int
    if (is.null(ci)) {
      stop("Cannot find confidence intervals in `x$conf_int`.")
    }

    ci <- as.data.frame(ci, check.names = FALSE)
    if (ncol(ci) < 2L) {
      stop("`x$conf_int` must have two columns.")
    }
    ci <- ci[, 1:2, drop = FALSE]
    names(ci) <- c("2.5%", "97.5%")
    ci[[1]] <- as.numeric(ci[[1]])
    ci[[2]] <- as.numeric(ci[[2]])

    if (!is.null(rn) && is.null(rownames(ci)) && nrow(ci) == length(rn)) {
      rownames(ci) <- rn
    }
    ci
  }

  coef_df <- .standardize_coef_df(x)
  conf_df <- .standardize_confint(x, rn = rownames(coef_df))

  trt_name <- if (!is.null(x$input) && !is.null(x$input$trt)) x$input$trt else NULL
  if (is.null(trt_name) || !nzchar(trt_name)) {
    trt_name <- rownames(coef_df)[1]
  }
  trt_row <- if (!is.null(rownames(coef_df)) && trt_name %in% rownames(coef_df)) {
    trt_name
  } else {
    rownames(coef_df)[1]
  }

  trt_est <- coef_df[trt_row, "logHR"]
  trt_lcl <- if (!is.null(rownames(conf_df)) && trt_row %in% rownames(conf_df)) conf_df[trt_row, 1] else conf_df[1, 1]
  trt_ucl <- if (!is.null(rownames(conf_df)) && trt_row %in% rownames(conf_df)) conf_df[trt_row, 2] else conf_df[1, 2]

  coef_tab <- cbind(
    logHR = coef_df$logHR,
    se = coef_df$se,
    HR = exp(coef_df$logHR)
  )
  rownames(coef_tab) <- rownames(coef_df)
  coef_tab <- .fmt_noquote(coef_tab, digits = digits)

  n_obs <- if (!is.null(x$n)) {
    as.integer(x$n[1])
  } else if (!is.null(x$surv_curve$n)) {
    as.integer(x$surv_curve$n[1])
  } else {
    NA_integer_
  }

  n_events <- if (!is.null(x$events)) {
    as.integer(x$events[1])
  } else if (!is.null(x$surv_curve$n.event)) {
    as.integer(sum(x$surv_curve$n.event, na.rm = TRUE))
  } else {
    NA_integer_
  }

  cat("Orthogonality method-of-moment IV estimation of Cox model\n")
  cat(strrep("=", 57), "\n")

  cat("Call:\n")
  print(x$call)
  cat(strrep("-", 57), "\n")

  if (!is.na(n_obs)) cat(sprintf("Sample size: %d\n", n_obs))
  if (!is.na(n_events)) cat(sprintf("Number of events: %d\n", n_events))
  cat(strrep("-", 57), "\n")

  cat(sprintf("Marginal treatment effect (%s):\n", trt_name))
  cat(sprintf("- Estimate: %.*f [%.*f, %.*f]\n", digits, trt_est, digits, trt_lcl, digits, trt_ucl))
  cat(sprintf(
    "- Hazard ratio: %.*f [%.*f, %.*f]\n",
    digits, exp(trt_est), digits, exp(trt_lcl), digits, exp(trt_ucl)
  ))
  cat(strrep("-", 57), "\n")

  cat("Coefficients:\n")
  .print_text(coef_tab, right = TRUE)
  cat(strrep("-", 57), "\n")

  if (!is.null(x$iv_diag$f_stat)) {
    cat("Nested exclusion test for weak IV:\n")
    cat(sprintf(
      "- F-statistics: %.*f%s\n",
      2, as.numeric(x$iv_diag$f_stat),
      ifelse(as.numeric(x$iv_diag$f_stat) > 10, " [> 10.0]", " [< 10.0]")
    ))
    cat(strrep("-", 57), "\n")
  }

  cat("A `coxivomom` object\n")
  invisible(x)
}
