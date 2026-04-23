#' Sequential 2SRI Cox analysis via sequential trials emulation
#'
#' Fits a sequential two-stage residual inclusion (2SRI) Cox
#' model for a time-varying treatment with a baseline instrumental variable,
#' using a stacked sequential trial emulation created by [seqem()].
#'
#' The analysis proceeds trial by trial. Within each emulated trial, a
#' first-stage logistic regression is fit for the trial-baseline treatment
#' `{tvtrt}.base` on the baseline IV and baseline covariates from
#' `trtformula`. A trial-specific control-function residual `R.base` is then
#' constructed and carried across all rows within that subject-trial trajectory.
#'
#' The stacked data are then artificially censored when observed treatment
#' deviates from the trial-baseline treatment status, and stabilised inverse
#' probability of artificial censoring weights (IPACW) are estimated. These
#' weights address the informative censoring induced by the sequential trial
#' emulation itself. This function does not estimate additional weights for
#' naturally occurring right censoring.
#'
#' Finally, a weighted composite Cox model is fit on the artificially censored,
#' weighted stacked data. The returned coefficient summary is organised in the
#' same style as [seqcox()], [coxiv_omom()], and [coxiv_tsrif()].
#'
#' @param formula a `formula` for the composite second-stage Cox model,
#'   typically of the form
#'   `Surv(start.new, stop.new, event) ~ tvtrt.base + x1.base + x2.base`.
#'   The right-hand side must include `paste0(tvtrt, ".base")`. The control
#'   function term `R.base` is added automatically if absent.
#' @param trtformula a `formula` for the trial-specific first-stage treatment
#'   model. Its left-hand side must equal `paste0(tvtrt, ".base")`, and its
#'   right-hand side must include `paste0(iv, ".base")`.
#' @param data a `seqem` object returned by [seqem()].
#' @param id a string name for the subject identifier variable.
#' @param tvtrt a string name for the original time-varying treatment
#'   variable.
#' @param iv a string name for the original baseline instrumental
#'   variable.
#' @param stratify logical; if `TRUE`, the composite Cox model is stratified by
#'   emulated trial so that each trial has its own baseline hazard.
#' @param se_type standard error estimator: one of `"robust"`, `"jackknife"`,
#'   `"bootstrap"`, or `"none"`.
#' @param B the number of bootstrap replicates when `se_type = "bootstrap"`.
#' @param seed optional random seed for bootstrap resampling.
#' @param ipacw_den_covs optional string vector of underlying covariate names
#'   to use in the denominator IPACW model through their lead versions (for
#'   example `c("endoleak")` uses `endoleak.lead1`). If `NULL`, the underlying
#'   variables corresponding to non-treatment `.base` terms in `formula` and
#'   `trtformula` are used whenever the original variables are available in the
#'   stacked data.
#' @param ipacw_num_covs optional string vector of baseline `.base`
#'   covariates to use in the numerator stabiliser model. If `NULL`, an
#'   intercept-only numerator model is used.
#' @param by_trial logical; if `TRUE`, fit and return separate second-stage Cox
#'   models within each emulated trial.
#' @param ... additional arguments passed to [survival::coxph()].
#'
#' @return A `coxivseq` object with components including:
#' \itemize{
#'   \item `call` : the matched call.
#'   \item `formula` : the composite second-stage Cox formula actually fitted.
#'   \item `trtformula` : the sequential first-stage treatment models formula.
#'   \item `fit` : the fitted composite [survival::coxph()] model.
#'   \item `coef` : a `data.frame` with columns `logHR`, `HR`, `se`, `Z`, `p.val`.
#'   \item `conf_int` : a `data.frame` of 95% confidence limits with columns
#'         `2.5%` and `97.5%` for estimated log hazard ratios.
#'   \item `vcov_mat` : variance-covariance matrix corresponding to the
#'         selected `se_type`.
#'   \item `trtfit_by_trial` : list of first-stage treatment models.
#'   \item `fit_by_trial` : list of second-stage trial-specific Cox fits when
#'         `by_trial = TRUE`.
#'   \item `data_ste` : final stacked, artificially censored analysis dataset,
#'         including `R.base` and the computed IPACW column.
#'   \item `jackknife`, `bootstrap`: a `matrix` of jackknife or bootstrap resampled
#'         estimates if `se_type` set to `jackknife` or `boostrap`.
#'   \item `tvtrt`, `iv`, `n`, `n_ids`, `n_trials`, `n_event`, `stratify`,
#'         `se_type`, `B`, `by_trial`: miscellaneous string tags for model print
#'         summary.
#' }
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
#' data("vascular", package = "surviv")
#'
#' vasc_seqem <- seqem(
#'   data = vascular,
#'   start = "time.start",
#'   stop = "time.stop",
#'   event = "event",
#'   id = "id",
#'   tvtrt = "reint",
#'   covs = c("iv", "diameter", "endoleak"),
#'   coarsen = "ceiling",
#'   cbin_width = 2
#' )
#'
#' fit_seqcox <- coxiv_seq(
#'   formula = survival::Surv(start.new, stop.new, event) ~
#'     reint.base + diameter.base + endoleak.base,
#'   trtformula = reint.base ~ iv.base + diameter.base + endoleak.base,
#'   data = vasc_seqem,
#'   id = "id",
#'   tvtrt = "reint",
#'   iv = "iv",
#'   se_type = "robust",
#'   by_trial = TRUE
#' )
#'
#' print(fit_seqcox)
#'
#' @import survival
#' @export
coxiv_seq <- function(formula,
                      trtformula,
                      data,
                      id,
                      tvtrt,
                      iv,
                      stratify = FALSE,
                      se_type = c("robust", "jackknife", "bootstrap", "none"),
                      B = 50,
                      seed = NULL,
                      ipacw_den_covs = NULL,
                      ipacw_num_covs = NULL,
                      by_trial = TRUE,
                      ...) {

  cl <- match.call()
  se_type <- match.arg(se_type)
  trial_id <- "trial"

  if (!inherits(data, "seqem")) {
    stop("`data` must be a `seqem` object produced by `seqem()`.")
  }

  dat_orig <- data$data_seqorig
  if (!is.data.frame(dat_orig)) {
    stop("`data` must contain `data_seqorig` from `seqem()`.")
  }

  needed_core <- c(id, trial_id, tvtrt, "rownum")
  missing_core <- setdiff(needed_core, names(dat_orig))
  if (length(missing_core) > 0L) {
    stop("Sequential trials data are missing required columns: ",
         paste(missing_core, collapse = ", "), ".")
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }
  if (!inherits(trtformula, "formula")) {
    stop("`trtformula` must be a formula.")
  }

  trt_base_var <- paste0(tvtrt, ".base")
  iv_base_var  <- paste0(iv, ".base")

  tt_out <- stats::terms(formula)
  rhs_terms_out <- attr(tt_out, "term.labels")
  rhs_vars_out  <- all.vars(formula[[3L]])

  if (!trt_base_var %in% rhs_terms_out && !trt_base_var %in% rhs_vars_out) {
    stop("The RHS of `formula` must include `", trt_base_var, "`.")
  }
  if (tvtrt %in% rhs_vars_out) {
    stop("`formula` uses `", tvtrt, "` without the required `.base` suffix. ",
         "Use `", trt_base_var, "` in the second-stage model.")
  }

  trt_lhs_var <- all.vars(trtformula[[2L]])
  if (length(trt_lhs_var) != 1L || !identical(trt_lhs_var, trt_base_var)) {
    stop("The left-hand side of `trtformula` must be exactly `", trt_base_var, "`.")
  }

  rhs_vars_trt <- all.vars(trtformula[[3L]])
  if (!iv_base_var %in% rhs_vars_trt) {
    stop("The RHS of `trtformula` must include `", iv_base_var, "`.")
  }
  if (iv %in% rhs_vars_trt) {
    stop("`trtformula` uses `", iv, "` without the required `.base` suffix. ",
         "Use `", iv_base_var, "` in the first-stage model.")
  }
  if (tvtrt %in% rhs_vars_trt) {
    warning("`trtformula` contains `", tvtrt, "` without `.base`. ",
            "First-stage predictors should generally use trial-baseline variables.")
  }

  remind_nonbase <- function(vars, context) {
    bad <- vars[vars %in% names(dat_orig) &
                  !grepl("\\.base$", vars) &
                  !vars %in% c(tvtrt, iv, id, trial_id, "R", "R.base")]
    bad <- unique(bad)
    if (length(bad) > 0L) {
      warning(
        context, " contains variables without a `.base` suffix: ",
        paste(bad, collapse = ", "),
        ". In sequential 2SRI, users should usually supply the trial-baseline ",
        "versions of covariates."
      )
    }
  }

  remind_nonbase(rhs_vars_out, "`formula`")
  remind_nonbase(rhs_vars_trt, "`trtformula`")

  base_terms_out <- rhs_terms_out[grepl("\\.base$", rhs_terms_out)]
  strip_base <- function(x) sub("\\.base$", "", x)

  ## ---- first-stage treatment models and control function -----------------
  if ("R.base" %in% names(dat_orig)) {
    warning("Overwriting existing `R.base` column in `data_seqorig`.")
  }
  dat_orig$R.base <- NA_real_

  trial_vec_orig <- dat_orig[[trial_id]]
  trial_vals <- sort(unique(trial_vec_orig))
  n_trials <- length(trial_vals)

  trtfit_by_trial <- vector("list", n_trials)
  names(trtfit_by_trial) <- paste0("trtfit_et", seq_len(n_trials))

  for (k in seq_len(n_trials)) {
    tr_val <- trial_vals[k]
    idx_tr <- trial_vec_orig == tr_val
    dat_k  <- dat_orig[idx_tr, , drop = FALSE]

    if (nrow(dat_k) == 0L) {
      trtfit_by_trial[[k]] <- NULL
      next
    }

    base_rows <- !duplicated(dat_k[[id]])
    dat_k_base <- dat_k[base_rows, , drop = FALSE]

    trt_fit_k <- try(
      stats::glm(trtformula, data = dat_k_base, family = stats::binomial()),
      silent = TRUE
    )

    if (inherits(trt_fit_k, "try-error")) {
      trtfit_by_trial[[k]] <- NULL
      next
    }

    trtfit_by_trial[[k]] <- trt_fit_k

    A_base_k <- dat_k_base[[trt_base_var]]
    p_hat_k  <- stats::fitted(trt_fit_k)
    R_k      <- A_base_k - p_hat_k

    ids_base_k <- dat_k_base[[id]]
    ids_all_k  <- dat_k[[id]]
    map_k      <- match(ids_all_k, ids_base_k)

    if (anyNA(map_k)) {
      stop("Failed to map first-stage residuals back to the stacked data for trial ",
           tr_val, ".")
    }

    dat_orig$R.base[idx_tr] <- R_k[map_k]
  }

  if (anyNA(dat_orig$R.base)) {
    warning("Some rows have missing `R.base` after first-stage estimation.")
  }

  dat_orig$R <- dat_orig$R.base

  ## ---- IPACW for artificial censoring only -------------------------------
  if (is.null(ipacw_den_covs)) {
    den_base_terms <- unique(c(base_terms_out, rhs_vars_trt[grepl("\\.base$", rhs_vars_trt)]))
    den_base_terms <- setdiff(den_base_terms, c(trt_base_var, iv_base_var, "R.base"))
    den_underlying <- strip_base(den_base_terms)
    den_underlying <- den_underlying[den_underlying %in% names(dat_orig)]
    tv_covs <- unique(c(den_underlying, iv, "R"))
    tv_covs <- tv_covs[tv_covs %in% names(dat_orig)]
  } else {
    tv_covs <- unique(c(ipacw_den_covs, iv, "R"))
    tv_covs <- tv_covs[tv_covs %in% names(dat_orig)]
  }

  baseline_covs_num <- NULL
  if (!is.null(ipacw_num_covs)) {
    baseline_covs_num <- ipacw_num_covs[ipacw_num_covs %in% names(dat_orig)]
  }

  dat <- .compute_ipacw_seq_iv(
    dat_seq_orig = dat_orig,
    id_var = id,
    trial_var = trial_id,
    trt_var = tvtrt,
    trt_base_var = trt_base_var,
    rownum_var = "rownum",
    tv_cov_lead = tv_covs,
    baseline_covs_base = baseline_covs_num
  )

  wt_col <- ".coxivseq_wt"
  dat[[wt_col]] <- as.numeric(dat$ipacw.stab.cum)

  ## ---- build Cox formulas ------------------------------------------------
  add_Rbase <- function(formula) {
    f_chr <- paste(deparse(formula), collapse = " ")
    if (grepl("R.base", f_chr, fixed = TRUE)) return(formula)
    parts <- strsplit(f_chr, "~", fixed = TRUE)[[1L]]
    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])
    rhs_new <- if (rhs == "" || rhs == "1") "R.base" else paste(rhs, "R.base", sep = " + ")
    stats::as.formula(paste(lhs, "~", rhs_new), env = environment(formula))
  }

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

  form_fit <- add_Rbase(formula)
  if (isTRUE(stratify)) {
    form_fit <- add_strata(form_fit, trial_id)
  }
  form_trial <- add_Rbase(formula)

  ## ---- composite fit -----------------------------------------------------
  fit_main <- survival::coxph(
    formula = form_fit,
    data = dat,
    weights = .coxivseq_wt,
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
      weights = .coxivseq_wt,
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
        survival::coxph(form_fit, data = dat_j, weights = .coxivseq_wt, ...),
        silent = TRUE
      )
      if (!inherits(fit_j, "try-error")) {
        jackknife_est[j, ] <- stats::coef(fit_j)
      }
    }

    valid_rows <- stats::complete.cases(jackknife_est)
    if (sum(valid_rows) > 1L) {
      jk <- jackknife_est[valid_rows, , drop = FALSE]
      m <- nrow(jk)
      jk_centered <- sweep(jk, 2L, coef_hat, FUN = "-")
      var_hat <- (m - 1) / m * crossprod(jk_centered)
    } else {
      warning("Too few valid jackknife refits; reverting to model-based variance.")
    }
  } else if (se_type == "bootstrap") {
    if (!is.null(seed)) set.seed(seed)

    id_vec <- dat[[id]]
    ids_unique <- unique(id_vec)
    n_id <- length(ids_unique)
    p <- length(coef_hat)
    bootstrap_est <- matrix(NA_real_, nrow = B, ncol = p)
    colnames(bootstrap_est) <- names(coef_hat)
    rows_by_id <- split(seq_len(nrow(dat)), id_vec)

    for (b in seq_len(B)) {
      ids_b <- sample(ids_unique, size = n_id, replace = TRUE)
      idx_b <- unlist(rows_by_id[match(ids_b, ids_unique)], use.names = FALSE)
      dat_b <- dat[idx_b, , drop = FALSE]
      fit_b <- try(
        survival::coxph(form_fit, data = dat_b, weights = .coxivseq_wt, ...),
        silent = TRUE
      )
      if (!inherits(fit_b, "try-error")) {
        bootstrap_est[b, ] <- stats::coef(fit_b)
      }
    }

    valid_rows <- stats::complete.cases(bootstrap_est)
    if (sum(valid_rows) > 1L) {
      var_hat <- stats::cov(bootstrap_est[valid_rows, , drop = FALSE])
    } else {
      warning("Too few valid bootstrap refits; reverting to model-based variance.")
    }
  }

  ## ---- trial-specific second-stage fits ---------------------------------
  fit_by_trial <- NULL
  if (isTRUE(by_trial)) {
    trial_vec <- dat[[trial_id]]
    trial_vals2 <- sort(unique(trial_vec))
    fit_by_trial <- vector("list", length(trial_vals2))
    names(fit_by_trial) <- paste0("fit_et", seq_along(trial_vals2))

    for (k in seq_along(trial_vals2)) {
      tr_val <- trial_vals2[k]
      dat_k <- dat[trial_vec == tr_val, , drop = FALSE]
      if (nrow(dat_k) < 1L) {
        fit_by_trial[[k]] <- NULL
        next
      }
      fit_k <- try(
        survival::coxph(form_trial, data = dat_k, weights = .coxivseq_wt, ...),
        silent = TRUE
      )
      if (inherits(fit_k, "try-error")) {
        fit_by_trial[[k]] <- NULL
      } else {
        fit_by_trial[[k]] <- fit_k
      }
    }
  }

  ## ---- summaries --------------------------------------------------------
  se_hat <- sqrt(diag(var_hat))
  z_val  <- coef_hat / se_hat
  p_val  <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)

  conf_int <- data.frame(
    `2.5%` = coef_hat - stats::qnorm(0.975) * se_hat,
    `97.5%` = coef_hat + stats::qnorm(0.975) * se_hat,
    check.names = FALSE
  )
  rownames(conf_int) <- names(coef_hat)

  coef_df <- data.frame(
    logHR = coef_hat,
    HR = exp(coef_hat),
    se = se_hat,
    Z = z_val,
    p.val = p_val,
    check.names = FALSE
  )
  rownames(coef_df) <- names(coef_hat)

  n_ids <- length(unique(dat[[id]]))
  n_trials2 <- length(unique(dat[[trial_id]]))

  res <- list(
    call = cl,
    formula = form_fit,
    trtformula = trtformula,
    tvtrt = tvtrt,
    iv = iv,
    fit = fit_used,
    coef = coef_df,
    conf_int = conf_int,
    vcov_mat = var_hat,
    n = nrow(dat),
    n_ids = n_ids,
    n_trials = n_trials2,
    n_event = if (!is.null(fit_used$nevent)) fit_used$nevent else sum(dat[[all.vars(formula[[2L]])[3L]]], na.rm = TRUE),
    stratify = stratify,
    se_type = se_type,
    B = if (se_type == "bootstrap") B else NULL,
    jackknife = jackknife_est,
    bootstrap = bootstrap_est,
    data_ste = dat,
    trtfit_by_trial = trtfit_by_trial,
    by_trial = by_trial,
    fit_by_trial = fit_by_trial
  )

  class(res) <- "coxivseq"
  res
}


.compute_ipacw_seq_iv <- function(dat_seq_orig,
                                  id_var,
                                  trial_var,
                                  trt_var,
                                  trt_base_var,
                                  rownum_var = "rownum",
                                  tv_cov_lead = NULL,
                                  baseline_covs_base = NULL) {

  stopifnot(id_var %in% names(dat_seq_orig),
            trial_var %in% names(dat_seq_orig),
            trt_var %in% names(dat_seq_orig),
            trt_base_var %in% names(dat_seq_orig),
            rownum_var %in% names(dat_seq_orig))

  ord <- order(dat_seq_orig[[id_var]], dat_seq_orig[[trial_var]], dat_seq_orig[[rownum_var]])
  dat_seq_orig <- dat_seq_orig[ord, , drop = FALSE]

  g <- interaction(dat_seq_orig[[id_var]], dat_seq_orig[[trial_var]], drop = TRUE)

  dat_seq_orig$A.lead1 <- ave(dat_seq_orig[[trt_var]], g,
                              FUN = function(x) c(x[-1L], NA))

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
                                     FUN = function(x) c(x[-1L], NA))
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

    form_denom <- if (length(lead_cov_names) > 0L) {
      stats::as.formula(paste("Anext.equal.to.baseline ~", paste(lead_cov_names, collapse = " + ")))
    } else {
      Anext.equal.to.baseline ~ 1
    }

    wt.mod.denom <- stats::glm(form_denom, family = stats::binomial(), data = est_data)
    dat_seq$wt.denom[mask_est] <- stats::predict(wt.mod.denom, newdata = est_data, type = "response")

    if (is.null(baseline_covs_base) || length(baseline_covs_base) == 0L) {
      form_num <- Anext.equal.to.baseline ~ 1
    } else {
      base_covs_in_data <- baseline_covs_base[baseline_covs_base %in% names(dat_seq)]
      if (length(base_covs_in_data) == 0L) {
        form_num <- Anext.equal.to.baseline ~ 1
      } else {
        rhs_num <- paste0(
          "as.factor(", rownum_var, ")*(",
          paste(base_covs_in_data, collapse = " + "),
          ")"
        )
        form_num <- stats::as.formula(paste("Anext.equal.to.baseline ~", rhs_num))
      }
    }

    wt.mod.num <- stats::glm(form_num, family = stats::binomial(), data = est_data)
    dat_seq$wt.num[mask_est] <- stats::predict(wt.mod.num, newdata = est_data, type = "response")
  }

  dat_seq$wt.denom <- ifelse(dat_seq[[trt_base_var]] == 1, 1, dat_seq$wt.denom)
  dat_seq$wt.num   <- ifelse(dat_seq[[trt_base_var]] == 1, 1, dat_seq$wt.num)

  dat_seq$ipacw.stab <- dat_seq$wt.num / dat_seq$wt.denom

  g2 <- interaction(dat_seq[[id_var]], dat_seq[[trial_var]], drop = TRUE)
  dat_seq$ipacw.stab.lag <- ave(dat_seq$ipacw.stab, g2,
                                FUN = function(x) c(1, x[-length(x)]))
  dat_seq$ipacw.stab.cum <- ave(dat_seq$ipacw.stab.lag, g2, FUN = cumprod)

  dat_seq
}
