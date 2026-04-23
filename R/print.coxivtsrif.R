#' Print the output from a `coxivtsrif` model fit.
#'
#' Provide succinct and organized output for a `coxivtsrif` model, such as
#' `coxiv_tsrif()`.
#'
#' @param x an object of class `"coxivtsrif"`, usually created by
#'   `coxiv_tsrif()`.
#' @param digits the number of decimal places to display for coefficient
#'   summaries and hazard ratios.
#' @param ... additional arguments (currently unused).
#'
#' @return The input object `x`, invisibly. Called for its side effect of printing a summary to the console.
#'
#' @examples
#' data("VitD", package = "surviv")
#' fit_tsrif <- surviv::coxiv_tsrif(
#'   surv = "time",
#'   cens = "death",
#'   covs = c("age"),
#'   trt = "vitd",
#'   iv = "filaggrin",
#'   tchange = NULL,
#'   data = VitD,
#'   fdist = "gaussian",
#'   bootvar = FALSE,
#'   B = 20,
#'   tvareff = FALSE
#' )
#' print(fit_tsrif)
#'
#' @export
print.coxivtsrif <- function(x, digits = 3L, ...) {
  if (!inherits(x, "coxivtsrif")) {
    stop("Input must be of class 'coxivtsrif'.")
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

  .extract_trt_coef <- function(obj) {
    trt_coef <- obj$trt_coef
    if (!is.null(trt_coef)) {
      trt_coef <- as.data.frame(trt_coef, check.names = FALSE)
      nm <- names(trt_coef)
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
      names(trt_coef) <- ifelse(is.na(nm2), nm, nm2)
      if (!all(c("logHR", "se") %in% names(trt_coef))) {
        stop("`x$trt_coef` must contain `logHR` and `se` columns.")
      }
      trt_coef$logHR <- as.numeric(trt_coef$logHR)
      trt_coef$se <- as.numeric(trt_coef$se)
      if ("Z" %in% names(trt_coef)) trt_coef$Z <- as.numeric(trt_coef$Z)
      if ("pval" %in% names(trt_coef)) trt_coef$pval <- as.numeric(trt_coef$pval)
      return(trt_coef)
    }

    # Backward compatibility with older objects
    est <- obj$est_coef
    if (is.null(est)) {
      stop("Cannot find treatment-effect summary in `x$trt_coef` or `x$est_coef`.")
    }

    ci <- obj$conf_int
    if (is.null(ci)) {
      stop("Cannot find confidence intervals in `x$conf_int`.")
    }

    if (is.matrix(ci) || is.data.frame(ci)) {
      ci_mat <- as.matrix(ci)
    } else {
      ci_mat <- matrix(ci, nrow = length(est), byrow = TRUE)
    }

    est_num <- as.numeric(est)
    if (is.null(names(est)) && !is.null(rownames(ci_mat))) names(est_num) <- rownames(ci_mat)

    se_num <- rep(NA_real_, length(est_num))
    if (!is.null(obj$variance)) {
      se_num <- sqrt(as.numeric(obj$variance))
    }

    z_num <- est_num / se_num
    p_num <- 2 * stats::pnorm(abs(z_num), lower.tail = FALSE)

    trt_coef <- data.frame(
      logHR = est_num,
      se = se_num,
      Z = z_num,
      pval = p_num,
      row.names = if (!is.null(names(est_num))) names(est_num) else NULL,
      check.names = FALSE
    )
    trt_coef
  }

  .extract_confint <- function(obj, rn = NULL) {
    ci <- obj$conf_int
    if (is.null(ci)) {
      stop("Cannot find confidence intervals in `x$conf_int`.")
    }

    if (is.vector(ci) && !is.list(ci)) {
      if (length(ci) %% 2L != 0L) stop("`x$conf_int` has invalid length.")
      ci <- matrix(ci, ncol = 2, byrow = TRUE)
    }

    ci <- as.data.frame(ci, check.names = FALSE)
    if (ncol(ci) < 2L) stop("`x$conf_int` must have two columns.")
    ci <- ci[, 1:2, drop = FALSE]
    names(ci) <- c("2.5%", "97.5%")
    ci[[1]] <- as.numeric(ci[[1]])
    ci[[2]] <- as.numeric(ci[[2]])

    if (!is.null(rn) && is.null(rownames(ci)) && nrow(ci) == length(rn)) {
      rownames(ci) <- rn
    }
    ci
  }

  .compact_cox_table <- function(model, trt_name = NULL, trt_override = NULL) {
    beta <- stats::coef(model)
    if (is.null(beta) || length(beta) == 0L) {
      return(NULL)
    }
    beta <- as.numeric(beta)
    names(beta) <- names(stats::coef(model))

    se <- tryCatch(
      sqrt(diag(stats::vcov(model))),
      error = function(e) sqrt(diag(model$var))
    )
    se <- as.numeric(se)
    names(se) <- names(stats::coef(model))

    tab <- data.frame(
      logHR = beta,
      se = se,
      HR = exp(beta),
      row.names = names(beta),
      check.names = FALSE
    )

    if (!is.null(trt_name) && !is.null(trt_override) && trt_name %in% rownames(tab)) {
      tab[trt_name, "logHR"] <- as.numeric(trt_override[["logHR"]])
      tab[trt_name, "se"] <- as.numeric(trt_override[["se"]])
      tab[trt_name, "HR"] <- exp(as.numeric(trt_override[["logHR"]]))
    }

    tab
  }

  .pick_n_events <- function(obj) {
    n_obs <- NA_integer_
    n_evt <- NA_integer_

    if (!is.null(obj$surv_curve$n)) {
      n_obs <- as.integer(obj$surv_curve$n[1])
    } else if (!is.null(obj$out_model$n)) {
      n_obs <- as.integer(obj$out_model$n[1])
    } else if (!is.null(obj$out_model1$n)) {
      n_obs <- as.integer(obj$out_model1$n[1])
    }

    if (!is.null(obj$surv_curve$n.event)) {
      n_evt <- as.integer(sum(obj$surv_curve$n.event, na.rm = TRUE))
    } else if (!is.null(obj$out_model$nevent)) {
      n_evt <- as.integer(obj$out_model$nevent[1])
    } else if (!is.null(obj$out_model1$nevent) && !is.null(obj$out_model2$nevent)) {
      # For time-varying models, `surv_curve$n.event` is preferred; otherwise use summed events.
      n_evt <- as.integer(obj$out_model1$nevent[1] + obj$out_model2$nevent[1])
    }

    list(n = n_obs, events = n_evt)
  }

  trt_name <- if (!is.null(x$input) && !is.null(x$input$trt)) x$input$trt else NULL
  trt_coef <- .extract_trt_coef(x)
  if (is.null(trt_name) || !nzchar(trt_name)) {
    trt_name <- rownames(trt_coef)[1]
    trt_name <- sub("\\.t[12]$", "", trt_name)
  }
  conf_df <- .extract_confint(x, rn = rownames(trt_coef))
  n_info <- .pick_n_events(x)

  cat("TSRI-Frailty instrumental variable analysis of Cox model\n")
  cat(strrep("=", 56), "\n")

  if (isTRUE(x$input$tvareff)) {
    cat("Call:\n")
    cat("- First-stage (treatment) model:\n")
    cat(paste(deparse(x$input$trtformula), collapse = " "), "\n")
    cat("- Second-stage (outcome) model:\n")
    cat(paste(deparse(x$input$formula), collapse = " "), "\n")
    cat(sprintf("* Time-varying treatment effect updates at t = %s\n", format(x$input$tchange)))
    cat(strrep("-", 56), "\n")

    if (!is.na(n_info$n)) cat(sprintf("Sample size: %d\n", n_info$n))
    if (!is.na(n_info$events)) cat(sprintf("Number of events: %d\n", n_info$events))
    cat(strrep("-", 56), "\n")

    rn <- rownames(trt_coef)
    if (is.null(rn) || nrow(trt_coef) < 2L) {
      rn <- c(paste0(trt_name, ".t1"), paste0(trt_name, ".t2"))
      rownames(trt_coef) <- rn
      if (nrow(conf_df) == 2L && is.null(rownames(conf_df))) rownames(conf_df) <- rn
    }

    cat(sprintf("Conditional treatment effect (%s):\n", trt_name))
    cat(sprintf("# Period 1 # (t<%.*f):\n", digits, as.numeric(x$input$tchange)))
    cat(sprintf(
      "- Estimate: %.*f [%.*f, %.*f]\n",
      digits, trt_coef[1, "logHR"], digits, conf_df[1, 1], digits, conf_df[1, 2]
    ))
    cat(sprintf(
      "- Hazard ratio: %.*f [%.*f, %.*f]\n",
      digits, exp(trt_coef[1, "logHR"]), digits, exp(conf_df[1, 1]), digits, exp(conf_df[1, 2])
    ))
    cat(strrep("~", 22), "\n")
    cat(sprintf("# Period 2 # (t>%.*f):\n", digits, as.numeric(x$input$tchange)))
    cat(sprintf(
      "- Estimate: %.*f [%.*f, %.*f]\n",
      digits, trt_coef[2, "logHR"], digits, conf_df[2, 1], digits, conf_df[2, 2]
    ))
    cat(sprintf(
      "- Hazard ratio: %.*f [%.*f, %.*f]\n",
      digits, exp(trt_coef[2, "logHR"]), digits, exp(conf_df[2, 1]), digits, exp(conf_df[2, 2])
    ))
    cat(strrep("-", 56), "\n")

    cat("Coefficients:\n")
    cat("# Period 1 #:\n")
    tab1 <- .compact_cox_table(x$out_model1, trt_name = trt_name, trt_override = trt_coef[1, ])
    if (!is.null(tab1)) .print_text(.fmt_noquote(tab1, digits = digits), right = TRUE)

    cat(strrep("~", 13), "\n")
    cat("# Period 2 #:\n")
    tab2 <- .compact_cox_table(x$out_model2, trt_name = trt_name, trt_override = trt_coef[2, ])
    if (!is.null(tab2)) .print_text(.fmt_noquote(tab2, digits = digits), right = TRUE)
    cat(strrep("-", 56), "\n")
  } else {
    cat("Call:\n")
    cat("- First-stage (treatment) model:\n")
    cat(paste(deparse(x$input$trtformula), collapse = " "), "\n")
    cat("- Second-stage (outcome) model:\n")
    cat(paste(deparse(x$input$formula), collapse = " "), "\n")
    cat("* Time-constant treatment effect\n")
    cat(strrep("-", 56), "\n")

    if (!is.na(n_info$n)) cat(sprintf("Sample size: %d\n", n_info$n))
    if (!is.na(n_info$events)) cat(sprintf("Number of events: %d\n", n_info$events))
    cat(strrep("-", 56), "\n")

    cat(sprintf("Conditional treatment effect (%s):\n", trt_name))
    cat(sprintf(
      "- Estimate: %.*f [%.*f, %.*f]\n",
      digits, trt_coef[1, "logHR"], digits, conf_df[1, 1], digits, conf_df[1, 2]
    ))
    cat(sprintf(
      "- Hazard ratio: %.*f [%.*f, %.*f]\n",
      digits, exp(trt_coef[1, "logHR"]), digits, exp(conf_df[1, 1]), digits, exp(conf_df[1, 2])
    ))
    cat(strrep("-", 56), "\n")

    cat("Coefficients:\n")
    tab <- .compact_cox_table(x$out_model, trt_name = trt_name, trt_override = trt_coef[1, ])
    if (!is.null(tab)) .print_text(.fmt_noquote(tab, digits = digits), right = TRUE)
    cat(strrep("-", 56), "\n")
  }

  if (!is.null(x$iv_diag$f_stat)) {
    cat("Nested exclusion test for weak IV:\n")
    cat(sprintf(
      "- F-statistics: %.2f%s\n",
      as.numeric(x$iv_diag$f_stat),
      ifelse(as.numeric(x$iv_diag$f_stat) > 10, " [> 10.0]", " [< 10.0]")
    ))
    cat(strrep("-", 56), "\n")
  }

  cat("A `coxivtsrif` object\n")
  invisible(x)
}
