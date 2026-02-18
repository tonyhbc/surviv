#' Sequential Cox analysis via sequential trial emulation (with IPACW)
#'
#' @description
#' Fit a Cox proportional hazards model to sequential trials data created by
#' [seqem()], using stabilized inverse probability of artificial censoring
#' weights (IPACW) in the style of Gran et al. and Keogh et al.
#'
#' @param formula Survival formula, e.g.
#'   `Surv(start.new, stop.new, event) ~ A.base + L1.base + L2.base`,
#'   where:
#'   \itemize{
#'     \item the first term ending in `.base` is treated as the baseline
#'           treatment, and
#'     \item remaining `.base` terms are baseline covariates used in the
#'           IPACW numerator model by default.
#'   }
#' @param data A `seqem` object returned by [seqem()]. It must contain both
#'   `data_seqorig` (pre-censoring stacked trials) and `data_seqem` (post-
#'   censoring stacked trials), as constructed in the current implementation.
#' @param id Name of the subject identifier.
#' @param trial_id Name of the trial identifier (default `"trial"`).
#' @param stratify Logical; if `TRUE` (default), adds `strata(trial_id)` to the
#'   Cox model so each emulated trial has its own baseline hazard.
#' @param se_type Type of standard error: `"jackknife"`, `"bootstrap"`,
#'   `"robust"`, or `"none"`.
#' @param B Number of bootstrap replicates if `se_type = "bootstrap"`.
#' @param seed Optional random seed for bootstrap.
#' @param ipcw Logical; if `TRUE` (default), estimate stabilized IPACW using
#'   Keogh-style logistic models. If `FALSE`, no IPACW is estimated unless
#'   `weight_var` is supplied.
#' @param ipcw_trunc Quantile at which to truncate cumulative IPACW
#'   (default `0.95`). Set to `NULL` to disable truncation.
#' @param ipcw_den_covs Optional character vector of **time-varying** covariate
#'   names whose lead values (e.g. `L.lead1`) should be used in the denominator
#'   model. If `NULL`, the underlying versions of `.base` covariates from
#'   `formula` are used (e.g. from `L.base` we use `L.lead1`).
#' @param ipcw_num_covs Optional character vector of **baseline** `.base`
#'   covariates to use in the numerator model. If `NULL`, all `.base` covariates
#'   except the baseline treatment are used.
#' @param weight_var Optional name of a weight column already in the stacked
#'   data to use in the Cox fit. If supplied, no IPACW is estimated.
#' @param ... Additional arguments passed to [survival::coxph()].
#'
#' @return An object of class `"seqcox"` with components:
#' \itemize{
#'   \item `call` – the matched call.
#'   \item `formula` – formula used in the main Cox fit (including strata).
#'   \item `coef`, `var`, `se`, `z`, `p`, `ci` – point estimates and Wald
#'         inference for the Cox coefficients.
#'   \item `n`, `nevent` – number of observations and events in the main fit.
#'   \item `fit` – the [survival::coxph()] object from the main fit (possibly
#'         with robust variance if `se_type = "robust"`).
#'   \item `weights_name` – name of the weight column used in the Cox model.
#'   \item `se_type`, `B`, `jackknife`, `bootstrap` – information about the
#'         resampling SEs if requested.
#'   \item `ipcw` – logical; whether IPACW was estimated inside `seqcox()`.
#' }
#'
#' @import survival
#' @export
seqcox <- function(formula,
                   data,
                   id,
                   trial_id = "trial",
                   stratify = TRUE,
                   se_type = c("jackknife", "bootstrap", "robust", "none"),
                   B = 200,
                   seed = NULL,
                   ipcw = TRUE,
                   ipcw_trunc = 0.95,
                   ipcw_den_covs = NULL,
                   ipcw_num_covs = NULL,
                   weight_var = NULL,
                   by_trial = FALSE,
                   ...) {

  cl      <- match.call()
  se_type <- match.arg(se_type)

  ## ---------------------- data handling ---------------------- ##

  if (!inherits(data, "seqem")) {
    stop("`data` must be a `seqem` object produced by `seqem()`.")
  }

  dat_orig <- data$data_seqorig
  dat_cens <- data$data_seqem

  if (!is.data.frame(dat_orig) || !is.data.frame(dat_cens)) {
    stop("`seqem` object must contain `data_seqorig` and `data_seqem` data.frames.")
  }

  if (!all(c(id, trial_id) %in% names(dat_cens))) {
    stop("`id` and `trial_id` must be columns of the sequential data.")
  }

  ## ---------------------- parse formula ---------------------- ##

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a survival formula.")
  }

  tt        <- stats::terms(formula)
  rhs_terms <- attr(tt, "term.labels")

  is_base    <- grepl("\\.base$", rhs_terms)
  base_terms <- rhs_terms[is_base]

  if (length(base_terms) == 0L) {
    stop("The RHS of `formula` must contain at least one `.base` variable. ",
         "The first such term is treated as baseline treatment.")
  }

  trt_base_var <- base_terms[1L]
  strip_base   <- function(x) sub("\\.base$", "", x)
  trt_var      <- strip_base(trt_base_var)

  if (!trt_var %in% names(dat_orig)) {
    stop("Time-varying treatment column `", trt_var, "` not found in `data_seqorig`.")
  }
  if (!trt_base_var %in% names(dat_orig)) {
    stop("Baseline treatment column `", trt_base_var, "` not found in `data_seqorig`.")
  }

  # Baseline covariates for numerator model, unless overridden
  if (is.null(ipcw_num_covs)) {
    baseline_cov_base <- setdiff(base_terms, trt_base_var)
  } else {
    baseline_cov_base <- ipcw_num_covs
  }

  # Time-varying covariates for denominator model (leads)
  if (is.null(ipcw_den_covs)) {
    tv_covs <- unique(strip_base(baseline_cov_base))
  } else {
    tv_covs <- ipcw_den_covs
  }

  ## ---------------------- prepare analysis data & weights ---------------- ##

  wt_col   <- ".seqcox_wt"   # column name for weights inside `dat`
  use_ipcw <- FALSE

  if (!is.null(weight_var)) {

    if (!weight_var %in% names(dat_cens)) {
      stop("Weight column '", weight_var, "' not found in `data_seqem`.")
    }
    dat           <- dat_cens
    dat[[wt_col]] <- as.numeric(dat[[weight_var]])

  } else if (isTRUE(ipcw)) {

    ipcw_obj <- .compute_ipcw_seq(
      dat_seq_orig       = dat_orig,
      id_var             = id,
      trial_var          = trial_id,
      trt_var            = trt_var,
      trt_base_var       = trt_base_var,
      rownum_var         = "rownum",
      tv_cov_lead        = tv_covs,
      baseline_covs_base = baseline_cov_base,
      trunc_prob         = ipcw_trunc
    )

    dat <- ipcw_obj$dat_seq
    if (!"ipw.stab.cum.trunc" %in% names(dat)) {
      stop("IPACW computation did not produce 'ipw.stab.cum.trunc'.")
    }
    dat[[wt_col]] <- as.numeric(dat$ipw.stab.cum.trunc)
    use_ipcw      <- TRUE

  } else {

    dat           <- dat_cens
    dat[[wt_col]] <- 1
  }

  if (!all(c(id, trial_id) %in% names(dat))) {
    stop("`id` and `trial_id` must be present in the analysis data.")
  }

  ## ---------------------- build Cox formula ---------------------- ##

  add_strata <- function(formula, trial_id) {
    f_chr <- paste(deparse(formula), collapse = " ")
    if (grepl("strata(", f_chr, fixed = TRUE)) return(formula)
    parts <- strsplit(f_chr, "~", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop("`formula` must have a left and right hand side.")
    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])
    rhs_new <- if (rhs == "" || rhs == "1") {
      sprintf("strata(%s)", trial_id)
    } else {
      paste(sprintf("strata(%s)", trial_id), rhs, sep = " + ")
    }
    stats::as.formula(paste(lhs, "~", rhs_new), env = environment(formula))
  }

  add_cluster <- function(formula, id_var) {
    f_chr <- paste(deparse(formula), collapse = " ")
    if (grepl("cluster(", f_chr, fixed = TRUE)) return(formula)
    parts <- strsplit(f_chr, "~", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop("`formula` must have a left and right hand side.")
    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])
    rhs_new <- if (rhs == "" || rhs == "1") {
      sprintf("cluster(%s)", id_var)
    } else {
      paste(sprintf("cluster(%s)", id_var), rhs, sep = " + ")
    }
    stats::as.formula(paste(lhs, "~", rhs_new), env = environment(formula))
  }

  # Composite-analysis formula (possibly with strata(trial))
  form_fit <- formula
  if (isTRUE(stratify)) {
    form_fit <- add_strata(form_fit, trial_id)
  }

  # Trial-specific formula: we use the original `formula` (no strata(trial) added)
  form_trial <- formula

  ## ---------------------- id & trial meta ---------------------- ##

  id_vec     <- dat[, id, drop = TRUE]
  ids_unique <- unique(id_vec)
  n_id       <- length(ids_unique)

  trial_vec   <- dat[, trial_id, drop = TRUE]
  trial_vals  <- sort(unique(trial_vec))
  n_trials    <- length(trial_vals)

  ## ---------------------- main Cox fit ---------------------- ##

  fit_main <- survival::coxph(
    form_fit,
    data    = dat,
    weights = .seqcox_wt,
    ...
  )

  coef_hat <- stats::coef(fit_main)
  if (is.null(coef_hat) || length(coef_hat) == 0L) {
    stop("No coefficients were estimated; check your formula and data.")
  }
  p       <- length(coef_hat)
  var_hat <- fit_main$var

  jackknife_est <- NULL
  boot_est      <- NULL
  fit_used      <- fit_main

  ## ---------------------- SE: robust / jackknife / bootstrap ------------ ##

  if (se_type == "robust") {

    form_rob <- add_cluster(form_fit, id)
    fit_rob  <- survival::coxph(
      form_rob,
      data    = dat,
      weights = .seqcox_wt,
      robust  = TRUE,
      ...
    )

    coef_hat <- stats::coef(fit_rob)
    var_hat  <- fit_rob$var
    fit_used <- fit_rob

  } else if (se_type == "jackknife") {

    jackknife_est <- matrix(NA_real_, nrow = n_id, ncol = p)
    colnames(jackknife_est) <- names(coef_hat)

    rows_by_id <- split(seq_len(nrow(dat)), id_vec)

    for (j in seq_len(n_id)) {
      id_j      <- ids_unique[j]
      keep_ids  <- ids_unique[ids_unique != id_j]
      keep_rows <- unlist(rows_by_id[match(keep_ids, ids_unique)], use.names = FALSE)

      dat_j <- dat[keep_rows, , drop = FALSE]

      fit_j <- try(
        survival::coxph(form_fit,
                        data    = dat_j,
                        weights = .seqcox_wt,
                        ...),
        silent = TRUE
      )
      if (!inherits(fit_j, "try-error")) {
        jackknife_est[j, ] <- stats::coef(fit_j)
      }
    }

    valid_rows <- stats::complete.cases(jackknife_est)
    if (sum(valid_rows) > 1L) {
      jk <- jackknife_est[valid_rows, , drop = FALSE]
      m  <- nrow(jk)
      jk_centered <- sweep(jk, 2L, coef_hat, FUN = "-")
      var_hat <- (m - 1) / m * crossprod(jk_centered)
    } else {
      warning("Too few valid jackknife fits; using model-based variance.")
    }

  } else if (se_type == "bootstrap") {

    if (!is.null(seed)) set.seed(seed)

    boot_est <- matrix(NA_real_, nrow = B, ncol = p)
    colnames(boot_est) <- names(coef_hat)

    rows_by_id <- split(seq_len(nrow(dat)), id_vec)

    for (b in seq_len(B)) {
      ids_b <- sample(ids_unique, size = n_id, replace = TRUE)
      idx_b <- unlist(rows_by_id[match(ids_b, ids_unique)], use.names = FALSE)
      dat_b <- dat[idx_b, , drop = FALSE]

      fit_b <- try(
        survival::coxph(form_fit,
                        data    = dat_b,
                        weights = .seqcox_wt,
                        ...),
        silent = TRUE
      )
      if (!inherits(fit_b, "try-error")) {
        boot_est[b, ] <- stats::coef(fit_b)
      }
    }

    valid_rows <- stats::complete.cases(boot_est)
    if (sum(valid_rows) > 1L) {
      var_hat <- stats::cov(boot_est[valid_rows, , drop = FALSE])
    } else {
      warning("Too few valid bootstrap fits; using model-based variance.")
    }
  }

  ## ---------------------- optional: by-trial fits ------------------------ ##

  fit_by_trial <- NULL
  if (isTRUE(by_trial)) {
    fit_by_trial <- vector("list", n_trials)
    names(fit_by_trial) <- paste0("fit_et", seq_len(n_trials))

    for (k in seq_len(n_trials)) {
      tr_val <- trial_vals[k]
      dat_k  <- dat[trial_vec == tr_val, , drop = FALSE]

      if (nrow(dat_k) < 1L) {
        fit_by_trial[[k]] <- NULL
        next
      }

      # Trial-specific Cox model with same IPACW weights, no strata(trial)
      fit_k <- try(
        survival::coxph(
          form_trial,
          data    = dat_k,
          weights = .seqcox_wt,
          ...
        ),
        silent = TRUE
      )

      if (inherits(fit_k, "try-error")) {
        fit_by_trial[[k]] <- NULL
      } else {
        fit_by_trial[[k]] <- fit_k
      }
    }
  }

  ## ---------------------- summary quantities ---------------------------- ##

  se_hat <- sqrt(diag(var_hat))
  zval   <- coef_hat / se_hat
  pval   <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

  ci_low  <- coef_hat - stats::qnorm(0.975) * se_hat
  ci_high <- coef_hat + stats::qnorm(0.975) * se_hat
  ci_mat  <- cbind(lower = ci_low, upper = ci_high)

  res <- list(
    call         = cl,
    formula      = form_fit,
    coef         = coef_hat,
    var          = var_hat,
    se           = se_hat,
    z            = zval,
    p            = pval,
    ci           = ci_mat,
    n            = fit_used$n,
    nevent       = fit_used$nevent,
    fit          = fit_used,
    se_type      = se_type,
    B            = if (se_type == "bootstrap") B else NULL,
    jackknife    = jackknife_est,
    bootstrap    = boot_est,
    ipcw         = use_ipcw,
    ipcw_trunc   = ipcw_trunc,
    # the actual stacked analysis dataset with weights
    data_ste     = dat,
    wt_col       = wt_col,
    # STE meta
    n_ids        = n_id,
    n_trials     = n_trials,
    trial_values = trial_vals,
    trial_id     = trial_id,
    stratify     = stratify,
    # by-trial fits
    by_trial     = by_trial,
    fit_by_trial = fit_by_trial
  )

  class(res) <- "seqcox"
  res
}



#' @export
print.seqcox <- function(x,
                         digits = max(3L, getOption("digits") - 3L),
                         ...) {
  if (!inherits(x, "seqcox")) {
    return(NextMethod())
  }

  cat("Sequential Cox analysis via sequential trials emulation\n")

  # STE meta
  if (!is.null(x$n_trials) && !is.null(x$trial_id)) {
    cat(sprintf("Composite analysis across %d emulated trials (trial variable: '%s')\n",
                x$n_trials, x$trial_id))
  }

  if (!is.null(x$stratify)) {
    cat("Stratified by trial in composite fit: ",
        if (isTRUE(x$stratify)) "yes\n" else "no\n", sep = "")
  }

  if (!is.null(x$n_ids)) {
    cat(sprintf("Number of individuals: %d\n", x$n_ids))
  }
  if (!is.null(x$n)) {
    cat(sprintf("Number of rows (composite dataset): %d\n", x$n))
  }
  if (!is.null(x$nevent)) {
    cat(sprintf("Number of events: %d\n", x$nevent))
  }

  if (!is.null(x$ipcw)) {
    if (isTRUE(x$ipcw)) {
      cat(sprintf("IPACW: yes (truncation quantile = %s)\n",
                  if (!is.null(x$ipcw_trunc)) as.character(x$ipcw_trunc) else "none"))
    } else {
      cat("IPACW: no\n")
    }
  }

  cat(sprintf("SE type: %s\n", x$se_type))

  if (isTRUE(x$by_trial) && !is.null(x$fit_by_trial)) {
    k_avail <- sum(!vapply(x$fit_by_trial, is.null, logical(1)))
    cat(sprintf("Trial-specific fits: %d model(s) available in $fit_by_trial\n",
                k_avail))
    cat("  (Access as fit_by_trial$fit_et1, ..., fit_et", length(x$fit_by_trial),
        "; mapping to actual trial values in $trial_values)\n", sep = "")
  }

  cat("\nCoefficients (composite fit):\n")

  tab <- cbind(
    coef     = x$coef,
    exp_coef = exp(x$coef),
    se       = x$se,
    lower_95 = exp(x$ci[, "lower"]),
    upper_95 = exp(x$ci[, "upper"]),
    p_value  = x$p
  )
  print(round(tab, digits = digits))

  invisible(x)
}


# Internal: compute stabilized IPACW for sequential trials
#
# - dat_seq_orig: stacked trials BEFORE artificial censoring, with columns:
#     id_var, trial_var, trt_var, trt_base_var, rownum_var
# - returns:
#     $dat_seq      : stacked trials AFTER artificial censoring, with
#                     ipw.stab, ipw.stab.cum, ipw.stab.cum.trunc
#     $dat_seq_orig : input with A.lead1 etc. added
.compute_ipcw_seq <- function(dat_seq_orig,
                              id_var,
                              trial_var,
                              trt_var,
                              trt_base_var,
                              rownum_var = "rownum",
                              tv_cov_lead = NULL,        # e.g. c("L") -> L.lead1 in denom
                              baseline_covs_base = NULL, # e.g. c("L.base") in numerator
                              trunc_prob = 0.95) {

  stopifnot(id_var %in% names(dat_seq_orig),
            trial_var %in% names(dat_seq_orig),
            trt_var %in% names(dat_seq_orig),
            trt_base_var %in% names(dat_seq_orig),
            rownum_var %in% names(dat_seq_orig))

  # Order for lead/lag
  ord <- order(dat_seq_orig[[id_var]],
               dat_seq_orig[[trial_var]],
               dat_seq_orig[[rownum_var]])
  dat_seq_orig <- dat_seq_orig[ord, , drop = FALSE]

  g <- interaction(dat_seq_orig[[id_var]], dat_seq_orig[[trial_var]], drop = TRUE)

  ## A.lead1 and Anext.equal.to.baseline ----------------------------------
  dat_seq_orig$A.lead1 <- ave(dat_seq_orig[[trt_var]], g,
                              FUN = function(x) c(x[-1], NA))

  dat_seq_orig$Anext.equal.to.baseline <- ifelse(
    is.na(dat_seq_orig$A.lead1),
    NA_integer_,
    ifelse(dat_seq_orig$A.lead1 == dat_seq_orig[[trt_base_var]], 1L, 0L)
  )

  ## Time-varying covariate leads (if requested) ---------------------------

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

  ## Impose artificial censoring -------------------------------------------

  keep_idx <- dat_seq_orig[[trt_var]] == dat_seq_orig[[trt_base_var]]
  dat_seq  <- dat_seq_orig[keep_idx, , drop = FALSE]

  dat_seq$wt.denom <- 1
  dat_seq$wt.num   <- 1

  # Rows contributing to weight models: baseline untreated and non-missing Anext
  mask_est <- (dat_seq[[trt_base_var]] == 0) &
    !is.na(dat_seq$Anext.equal.to.baseline)

  if (any(mask_est)) {
    est_data <- dat_seq[mask_est, , drop = FALSE]

    ## Denominator: P(Anext==1 | covs at next visit, A.base=0)
    if (length(lead_cov_names) > 0L) {
      denom_rhs  <- paste(lead_cov_names, collapse = " + ")
      form_denom <- stats::as.formula(
        paste("Anext.equal.to.baseline ~", denom_rhs)
      )
    } else {
      form_denom <- Anext.equal.to.baseline ~ 1
    }

    wt.mod.denom <- stats::glm(form_denom,
                               family = stats::binomial(),
                               data   = est_data)

    dat_seq$wt.denom[mask_est] <- stats::predict(
      wt.mod.denom,
      newdata = est_data,
      type    = "response"
    )

    ## Numerator: P(Anext==1 | as.factor(rownum)*baseline_covs, A.base=0)
    base_covs_in_data <- baseline_covs_base[baseline_covs_base %in% names(dat_seq)]

    if (length(base_covs_in_data) > 0L) {
      num_rhs <- paste0(
        "as.factor(", rownum_var, ")*(",
        paste(base_covs_in_data, collapse = " + "),
        ")"
      )
    } else {
      num_rhs <- paste0("as.factor(", rownum_var, ")")
    }

    form_num <- stats::as.formula(
      paste("Anext.equal.to.baseline ~", num_rhs)
    )

    wt.mod.num <- stats::glm(form_num,
                             family = stats::binomial(),
                             data   = est_data)

    dat_seq$wt.num[mask_est] <- stats::predict(
      wt.mod.num,
      newdata = est_data,
      type    = "response"
    )
  }

  # For initiators at baseline (treated at trial start), weights = 1
  dat_seq$wt.denom <- ifelse(dat_seq[[trt_base_var]] == 1, 1, dat_seq$wt.denom)
  dat_seq$wt.num   <- ifelse(dat_seq[[trt_base_var]] == 1, 1, dat_seq$wt.num)

  # Period stabilized weights
  dat_seq$ipw.stab <- dat_seq$wt.num / dat_seq$wt.denom

  ## Cumulative stabilized weights (IPACW) ---------------------------------

  g2 <- interaction(dat_seq[[id_var]], dat_seq[[trial_var]], drop = TRUE)

  dat_seq$ipw.stab.lag <- ave(dat_seq$ipw.stab, g2,
                              FUN = function(x) c(1, x[-length(x)]))
  dat_seq$ipw.stab.cum <- ave(dat_seq$ipw.stab.lag, g2,
                              FUN = cumprod)

  if (!is.null(trunc_prob)) {
    pct <- stats::quantile(dat_seq$ipw.stab.cum,
                           probs = trunc_prob,
                           na.rm = TRUE)
    dat_seq$ipw.stab.cum.trunc <- pmin(dat_seq$ipw.stab.cum, pct)
  } else {
    dat_seq$ipw.stab.cum.trunc <- dat_seq$ipw.stab.cum
  }

  list(
    dat_seq      = dat_seq,
    dat_seq_orig = dat_seq_orig,
    lead_covs    = lead_cov_names
  )
}

