#' Sequential Cox analysis via sequential trials emulation
#'
#' @description
#' Fit a composite Cox model to a stacked sequential trial emulation object
#' produced by [seqem()]. The analysis follows the sequential Cox approach of
#' Gran et al. by using stabilized inverse probability of artificial censoring
#' weights (IPACW) to adjust for the informative artificial censoring induced by
#' restricting each emulated trial to individuals who continue to follow their
#' treatment assignment at the trial start.
#'
#' @param formula a survival formula for the composite Cox model, typically of
#'   the form `Surv(start.new, stop.new, event) ~ trt.base + L.base + ...`.
#'   The *first* `.base` term is treated as the trial-baseline treatment variable.
#' @param data a `seqem` object returned by [seqem()]. A raw data frame is not
#'   accepted.
#' @param id a string name of the subject identifier column in the sequential trials data.
#' @param trial_id a string name of the trial index column. Defaults to `"trial"`.
#' @param stratify logical; if `TRUE`, the composite Cox model is stratified by
#'   `trial_id` so each emulated trial has its own baseline hazard.
#' @param se_type standard error estimator: one of `"robust"`, `"jackknife"`,
#'   `"bootstrap"`, or `"none"`.
#' @param B number of bootstrap replicates when `se_type = "bootstrap"`.
#' @param seed optional random seed for the bootstrap.
#' @param ipacw_trunc optional truncation quantile for cumulative stabilized
#'   IPACW. Defaults to `0.95`. Set to `NULL` to disable truncation.
#' @param ipacw_den_covs optional character vector of underlying covariate names
#'   to use in the denominator IPACW model through their lead versions (for
#'   example `c("endoleak")` uses `endoleak.lead1`). If `NULL`, the underlying
#'   variables corresponding to non-treatment `.base` terms in `formula` are
#'   used.
#' @param ipacw_num_covs optional character vector of baseline `.base`
#'   covariates to use in the numerator stabilizer model. If `NULL`, an
#'   intercept-only numerator model is used.
#' @param by_trial logical; if `TRUE` (default), fit trial-specific Cox models
#'   and return them in `$fit_by_trial`.
#' @param ... additional arguments passed to [survival::coxph()].
#'
#' @return A `seqcox` object with components including:
#' * `call` : the matched call.
#' * `formula` : the composite Cox formula actually fitted.
#' * `fit` : a composite [survival::coxph()] fit.
#' * `coef` : a `data.frame` of estimated coefficient summaries with columns
#'          `logHR`, `HR`, `se`, `Z`, and `p.val`.
#' * `conf_int` : a `data.frame` of 95% confidence limits with columns
#'         `2.5%` and `97.5%`.
#' * `n_ids`, `n_trials`, `n_event`: the analysis size summaries
#'         for sample size, number of emulated trials aggregated, and number
#'         of events.
#' * `fit_by_trial` : a list of trial-specific Cox fits when
#'         `by_trial = TRUE`.
#' * `data_ste` : a `data.frame` of stacked analysis dataset actually used in
#'         the composite fit, including the computed IPACW column.
#' * `jackknife`, `bootstrap`: a vector of jackknife or bootstrap resampled
#'         estimates if `se_type` set to `jackknife` or `boostrap`.
#'
#' @references
#' Gran JM, Roysland K, Wolbers M, Didelez V, Sterne JA, Ledergerber B,
#' Furrer H, von Wyl V, Aalen OO. A sequential Cox approach for estimating the
#' causal effect of treatment in the presence of time-dependent confounding
#' applied to data from the Swiss HIV Cohort Study. \emph{Statistics in
#' Medicine}. 2010;29(26):2757-2768.
#'
#' Keogh RH, Gran JM, Seaman SR, Davies G, Vansteelandt S. Causal inference in
#' survival analysis using longitudinal observational data: Sequential trials
#' and marginal structural models. \emph{Statistics in Medicine}.
#' 2023;42(13):2191-2225.
#'
#' @examples
#'   data("vascular", package = "surviv")
#'   vasc_seqem = seqem(data = vascular, start = "time.start",
#'                      stop = "time.stop", event = "event", id = "id",
#'                      tvtrt = "reint", covs = c("diameter", "endoleak"),
#'                      coarsen = "ceiling", cbin_width = 2)
#'   fit_seqcox <- seqcox(
#'     formula = survival::Surv(start.new, stop.new, event) ~
#'       reint.base + diameter.base + endoleak.base,
#'     data = vasc_seqem,
#'     id = "id",
#'     se_type = "robust"
#'   )
#'   print(fit_seqcox)
#'
#' @import survival
#' @export
seqcox <- function(formula,
                   data,
                   id,
                   trial_id = "trial",
                   stratify = FALSE,
                   se_type = c("robust", "jackknife", "bootstrap", "none"),
                   B = 50,
                   seed = NULL,
                   ipacw_trunc = 0.95,
                   ipacw_den_covs = NULL,
                   ipacw_num_covs = NULL,
                   by_trial = TRUE,
                   ...) {

  cl <- match.call()
  se_type <- match.arg(se_type)

  if (!inherits(data, "seqem")) {
    stop("`data` must be a `seqem` object produced by `seqem()`.")
  }

  dat_orig <- data$data_seqorig
  if (!is.data.frame(dat_orig)) {
    stop("`data` must contain `data_seqorig` from `seqem()`.")
  }

  needed_core <- c(id, trial_id, "rownum")
  if (!all(needed_core %in% names(dat_orig))) {
    stop("Sequential trials data are missing one or more required columns: ",
         paste(setdiff(needed_core, names(dat_orig)), collapse = ", "))
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }

  tt <- stats::terms(formula)
  rhs_terms <- attr(tt, "term.labels")
  base_terms <- rhs_terms[grepl("\\.base$", rhs_terms)]

  if (length(base_terms) == 0L) {
    stop("The RHS of `formula` must contain at least one `.base` term.")
  }

  trt_base_var <- base_terms[1L]
  strip_base <- function(x) sub("\\.base$", "", x)
  trt_var <- strip_base(trt_base_var)

  if (!all(c(trt_var, trt_base_var) %in% names(dat_orig))) {
    stop("Could not find treatment columns `", trt_var, "` and/or `",
         trt_base_var, "` in `data_seqorig`.")
  }

  # Baseline covariates appearing in the outcome model, excluding treatment.
  baseline_cov_base <- setdiff(base_terms, trt_base_var)

  # Denominator model should use time-varying lead versions.
  if (is.null(ipacw_den_covs)) {
    tv_covs <- unique(strip_base(baseline_cov_base))
  } else {
    tv_covs <- ipacw_den_covs
  }

  dat <- .compute_ipacw_seq(
    dat_seq_orig = dat_orig,
    id_var = id,
    trial_var = trial_id,
    trt_var = trt_var,
    trt_base_var = trt_base_var,
    rownum_var = "rownum",
    tv_cov_lead = tv_covs,
    baseline_covs_base = ipacw_num_covs,
    trunc_prob = ipacw_trunc
  )

  if (!"ipacw.stab.cum.trunc" %in% names(dat)) {
    stop("IPACW computation failed to produce `ipacw.stab.cum.trunc`.")
  }

  dat$.seqcox_wt <- as.numeric(dat$ipacw.stab.cum.trunc)

  add_strata <- function(formula, trial_id) {
    f_chr <- paste(deparse(formula), collapse = " ")
    if (grepl("strata(", f_chr, fixed = TRUE)) return(formula)
    parts <- strsplit(f_chr, "~", fixed = TRUE)[[1L]]
    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])
    rhs_new <- if (rhs == "" || rhs == "1") {
      sprintf("strata(%s)", trial_id)
    } else {
      paste(rhs, sprintf("strata(%s)", trial_id), sep = " + ")
    }
    stats::as.formula(paste(lhs, "~", rhs_new), env = environment(formula))
  }

  add_cluster <- function(formula, id_var) {
    f_chr <- paste(deparse(formula), collapse = " ")
    if (grepl("cluster(", f_chr, fixed = TRUE)) return(formula)
    parts <- strsplit(f_chr, "~", fixed = TRUE)[[1L]]
    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])
    rhs_new <- if (rhs == "" || rhs == "1") {
      sprintf("cluster(%s)", id_var)
    } else {
      paste(rhs, sprintf("cluster(%s)", id_var), sep = " + ")
    }
    stats::as.formula(paste(lhs, "~", rhs_new), env = environment(formula))
  }

  form_fit <- formula
  if (isTRUE(stratify)) {
    form_fit <- add_strata(form_fit, trial_id)
  }
  form_trial <- formula

  fit_main <- survival::coxph(
    formula = form_fit,
    data = dat,
    weights = .seqcox_wt,
    ...
  )

  coef_hat <- stats::coef(fit_main)
  if (is.null(coef_hat) || length(coef_hat) == 0L) {
    stop("No coefficients were estimated by the composite Cox model.")
  }

  var_hat <- fit_main$var
  fit_used <- fit_main
  jackknife_est <- NULL
  bootstrap_est <- NULL

  if (se_type == "robust") {
    form_rob <- add_cluster(form_fit, id)
    fit_rob <- survival::coxph(
      formula = form_rob,
      data = dat,
      weights = .seqcox_wt,
      robust = TRUE,
      ...
    )
    coef_hat <- stats::coef(fit_rob)
    var_hat <- fit_rob$var
    fit_used <- fit_rob
  } else if (se_type == "jackknife") {
    id_vec <- dat[[id]]
    ids_unique <- unique(id_vec)
    p <- length(coef_hat)
    jackknife_est <- matrix(NA_real_, nrow = length(ids_unique), ncol = p)
    colnames(jackknife_est) <- names(coef_hat)
    rows_by_id <- split(seq_len(nrow(dat)), id_vec)

    for (j in seq_along(ids_unique)) {
      keep_rows <- unlist(rows_by_id[setdiff(names(rows_by_id), ids_unique[j])],
                          use.names = FALSE)
      dat_j <- dat[keep_rows, , drop = FALSE]
      fit_j <- try(
        survival::coxph(form_fit, data = dat_j, weights = .seqcox_wt, ...),
        silent = TRUE
      )
      if (!inherits(fit_j, "try-error")) {
        jackknife_est[j, ] <- stats::coef(fit_j)
      }
    }

    ok <- stats::complete.cases(jackknife_est)
    if (sum(ok) > 1L) {
      jk <- jackknife_est[ok, , drop = FALSE]
      m <- nrow(jk)
      jk_centered <- sweep(jk, 2L, coef_hat, FUN = "-")
      var_hat <- (m - 1) / m * crossprod(jk_centered)
    } else {
      warning("Too few successful jackknife fits; using model-based variance.")
    }
  } else if (se_type == "bootstrap") {
    if (!is.null(seed)) set.seed(seed)
    id_vec <- dat[[id]]
    ids_unique <- unique(id_vec)
    p <- length(coef_hat)
    bootstrap_est <- matrix(NA_real_, nrow = B, ncol = p)
    colnames(bootstrap_est) <- names(coef_hat)
    rows_by_id <- split(seq_len(nrow(dat)), id_vec)

    for (b in seq_len(B)) {
      ids_b <- sample(ids_unique, size = length(ids_unique), replace = TRUE)
      idx_b <- unlist(rows_by_id[match(ids_b, ids_unique)], use.names = FALSE)
      dat_b <- dat[idx_b, , drop = FALSE]
      fit_b <- try(
        survival::coxph(form_fit, data = dat_b, weights = .seqcox_wt, ...),
        silent = TRUE
      )
      if (!inherits(fit_b, "try-error")) {
        bootstrap_est[b, ] <- stats::coef(fit_b)
      }
    }

    ok <- stats::complete.cases(bootstrap_est)
    if (sum(ok) > 1L) {
      var_hat <- stats::cov(bootstrap_est[ok, , drop = FALSE])
    } else {
      warning("Too few successful bootstrap fits; using model-based variance.")
    }
  }

  se_hat <- sqrt(diag(var_hat))
  z_hat <- coef_hat / se_hat
  p_hat <- 2 * stats::pnorm(abs(z_hat), lower.tail = FALSE)

  conf_low <- coef_hat - stats::qnorm(0.975) * se_hat
  conf_high <- coef_hat + stats::qnorm(0.975) * se_hat

  coef_df <- data.frame(
    logHR = as.numeric(coef_hat),
    HR = exp(as.numeric(coef_hat)),
    se = as.numeric(se_hat),
    Z = as.numeric(z_hat),
    p.val = as.numeric(p_hat),
    row.names = names(coef_hat),
    check.names = FALSE
  )

  conf_int <- data.frame(
    "2.5%" = as.numeric(conf_low),
    "97.5%" = as.numeric(conf_high),
    row.names = names(coef_hat),
    check.names = FALSE
  )

  fit_by_trial <- vector("list", 0L)
  if (isTRUE(by_trial)) {
    trial_vals <- sort(unique(dat[[trial_id]]))
    fit_by_trial <- vector("list", length(trial_vals))
    names(fit_by_trial) <- paste0("fit_et", seq_along(trial_vals))

    for (k in seq_along(trial_vals)) {
      dat_k <- dat[dat[[trial_id]] == trial_vals[k], , drop = FALSE]
      fit_k <- try(
        survival::coxph(form_trial, data = dat_k, weights = .seqcox_wt, ...),
        silent = TRUE
      )
      fit_by_trial[[k]] <- if (inherits(fit_k, "try-error")) NULL else fit_k
    }
  }

  out <- list(
    call = cl,
    formula = form_fit,
    fit = fit_used,
    coef = coef_df,
    conf_int = conf_int,
    n = fit_used$n,
    n_event = fit_used$nevent,
    n_ids = length(unique(dat[[id]])),
    n_trials = length(unique(dat[[trial_id]])),
    trial_id = trial_id,
    jackknife = jackknife_est,
    bootstrap = bootstrap_est,
    fit_by_trial = fit_by_trial,
    data_ste = dat
  )

  class(out) <- "seqcox"
  out
}


# Internal helper: stabilized IPACW for sequential trials data
.compute_ipacw_seq <- function(dat_seq_orig,
                               id_var,
                               trial_var,
                               trt_var,
                               trt_base_var,
                               rownum_var = "rownum",
                               tv_cov_lead = NULL,
                               baseline_covs_base = NULL,
                               trunc_prob = 0.95) {

  stopifnot(id_var %in% names(dat_seq_orig),
            trial_var %in% names(dat_seq_orig),
            trt_var %in% names(dat_seq_orig),
            trt_base_var %in% names(dat_seq_orig),
            rownum_var %in% names(dat_seq_orig))

  ord <- order(dat_seq_orig[[id_var]],
               dat_seq_orig[[trial_var]],
               dat_seq_orig[[rownum_var]])
  dat_seq_orig <- dat_seq_orig[ord, , drop = FALSE]

  g <- interaction(dat_seq_orig[[id_var]], dat_seq_orig[[trial_var]], drop = TRUE)

  dat_seq_orig$A.lead1 <- ave(dat_seq_orig[[trt_var]], g,
                              FUN = function(x) c(x[-1], NA))

  dat_seq_orig$Anext.equal.to.baseline <- ifelse(
    is.na(dat_seq_orig$A.lead1),
    NA_integer_,
    ifelse(dat_seq_orig$A.lead1 == dat_seq_orig[[trt_base_var]], 1L, 0L)
  )

  lead_cov_names <- character(0L)
  if (!is.null(tv_cov_lead) && length(tv_cov_lead) > 0L) {
    for (v in tv_cov_lead) {
      if (!v %in% names(dat_seq_orig)) next
      newname <- paste0(v, ".lead1")
      dat_seq_orig[[newname]] <- ave(dat_seq_orig[[v]], g,
                                     FUN = function(x) c(x[-1], NA))
      lead_cov_names <- c(lead_cov_names, newname)
    }
  }

  keep_idx <- dat_seq_orig[[trt_var]] == dat_seq_orig[[trt_base_var]]
  dat_seq <- dat_seq_orig[keep_idx, , drop = FALSE]

  dat_seq$wt.denom <- 1
  dat_seq$wt.num <- 1

  mask_est <- (dat_seq[[trt_base_var]] == 0) & !is.na(dat_seq$Anext.equal.to.baseline)

  if (any(mask_est)) {
    est_data <- dat_seq[mask_est, , drop = FALSE]

    # Denominator uses lead versions of time-varying covariates when available.
    if (length(lead_cov_names) > 0L) {
      form_denom <- stats::as.formula(
        paste("Anext.equal.to.baseline ~", paste(lead_cov_names, collapse = " + "))
      )
    } else {
      form_denom <- Anext.equal.to.baseline ~ 1
    }

    wt.mod.denom <- stats::glm(form_denom,
                               family = stats::binomial(),
                               data = est_data)

    dat_seq$wt.denom[mask_est] <- stats::predict(
      wt.mod.denom,
      newdata = est_data,
      type = "response"
    )

    # By default use an intercept-only numerator stabilizer.
    base_covs_in_data <- NULL
    if (!is.null(baseline_covs_base)) {
      base_covs_in_data <- baseline_covs_base[baseline_covs_base %in% names(est_data)]
    }

    if (is.null(base_covs_in_data) || length(base_covs_in_data) == 0L) {
      form_num <- Anext.equal.to.baseline ~ 1
    } else {
      form_num <- stats::as.formula(
        paste("Anext.equal.to.baseline ~", paste(base_covs_in_data, collapse = " + "))
      )
    }

    wt.mod.num <- stats::glm(form_num,
                             family = stats::binomial(),
                             data = est_data)

    dat_seq$wt.num[mask_est] <- stats::predict(
      wt.mod.num,
      newdata = est_data,
      type = "response"
    )
  }

  dat_seq$wt.denom <- ifelse(dat_seq[[trt_base_var]] == 1, 1, dat_seq$wt.denom)
  dat_seq$wt.num <- ifelse(dat_seq[[trt_base_var]] == 1, 1, dat_seq$wt.num)

  dat_seq$ipacw.stab <- dat_seq$wt.num / dat_seq$wt.denom

  g2 <- interaction(dat_seq[[id_var]], dat_seq[[trial_var]], drop = TRUE)
  dat_seq$ipacw.stab.lag <- ave(dat_seq$ipacw.stab, g2,
                                FUN = function(x) c(1, x[-length(x)]))
  dat_seq$ipacw.stab.cum <- ave(dat_seq$ipacw.stab.lag, g2, FUN = cumprod)

  if (!is.null(trunc_prob)) {
    pct <- stats::quantile(dat_seq$ipacw.stab.cum,
                           probs = trunc_prob,
                           na.rm = TRUE)
    dat_seq$ipacw.stab.cum.trunc <- pmin(dat_seq$ipacw.stab.cum, pct)
  } else {
    dat_seq$ipacw.stab.cum.trunc <- dat_seq$ipacw.stab.cum
  }

  dat_seq
}
