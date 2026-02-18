#' Sequential 2SRI Cox model via trial emulation with a baseline IV
#'
#' @description
#' `coxiv_seq()` implements a sequential two–stage residual inclusion (2SRI)
#' Cox model for time‐to‐event outcomes with a \emph{time–varying} treatment
#' and a \emph{baseline} instrumental variable (IV), using the sequential trial
#' emulation framework of Gran et al. and Keogh et al.
#'
#' The function takes as input a `seqem` object (constructed by [seqem()])
#' and proceeds in four key steps:
#'
#' \enumerate{
#'   \item \strong{Sequential trial emulation.}
#'   Using the preprocessed data from [seqem()], each individual contributes
#'   person–time to a series of emulated trials defined by distinct treatment
#'   initiation times. In each trial, only covariates at the trial baseline
#'   (variables with suffix \code{.base}) are allowed to enter the models.
#'
#'   \item \strong{First–stage treatment models and control function.}
#'   For each emulated trial, a logistic regression specified by \code{trtformula}
#'   is fit for the trial–baseline treatment status (e.g. \code{A.base}) on the
#'   baseline IV and baseline covariates. The fitted probabilities define a
#'   trial–specific control function residual
#'   \eqn{R\_k = A\_k^{\text{base}} - \hat{\Pr}(A\_k^{\text{base}} = 1 \mid Z, X)},
#'   which is then mapped back to all rows within that trial as a baseline–only
#'   covariate \code{R.base}. The collection of first–stage models is returned
#'   in \code{$trtfit_by_trial} for IV diagnostics.
#'
#'   \item \strong{Artificial censoring and inverse probability weights.}
#'   Within each emulated trial, follow–up is artificially censored at the time
#'   of treatment switching away from the baseline treatment status. Stabilised
#'   inverse probability artificial censoring weights (IPACW) are then estimated
#'   via logistic models for the probability of remaining under the same
#'   treatment, and accumulated over follow–up to form trial– and subject–specific
#'   cumulative weights. The truncated cumulative weights are used as Cox
#'   regression weights and are stored in the returned stacked dataset.
#'
#'   \item \strong{Second–stage Cox model and composite estimation.}
#'   A weighted Cox proportional hazards model is fitted on the stacked,
#'   artificially censored data, with the trial–baseline treatment, baseline
#'   covariates and the control–function covariate \code{R.base} on the right–hand
#'   side (plus optional stratification by trial). The resulting hazard–ratio
#'   estimates are obtained under a composite–likelihood interpretation across
#'   emulated trials. Standard errors can be obtained via robust (clustered on
#'   subject), jackknife (leave–one–subject–out) or bootstrap resampling over
#'   subjects. Optionally, if \code{by_trial = TRUE}, separate 2SRI–adjusted Cox
#'   models are also fitted within each emulated trial and returned in
#'   \code{$fit_by_trial}.
#' }
#'
#' @details
#' The method assumes:
#' \itemize{
#'   \item a single binary time–varying treatment \code{tvtrt} that may switch
#'         from 0 to 1 over follow–up;
#'   \item a baseline‐only instrumental variable \code{iv} that affects treatment
#'         assignment but, conditional on baseline covariates, has no direct
#'         effect on the event time and is independent of unmeasured
#'         prognostic factors;
#'   \item all measured confounders included in the models enter as baseline
#'         covariates (i.e. as \code{*.base} variables), and treatment is the
#'         only time–varying predictor;
#'   \item artificial censoring and IPACW correctly adjust for the informative
#'         censoring induced by restricting follow–up within each emulated trial
#'         to periods in which the observed treatment remains equal to the
#'         trial‐baseline treatment status.
#' }
#'
#' The function is designed to mirror the structure of [seqcox()], but with an
#' additional first–stage IV step and inclusion of the trial–specific control
#' function \code{R.base} in the second–stage Cox model.
#'
#' @param formula A survival formula of the form
#'   \code{Surv(start, stop, event) ~ A.base + X1.base + ...}, specifying the
#'   second–stage Cox model. The right–hand side must include the baseline
#'   treatment term corresponding to \code{tvtrt.base}. Any additional
#'   \code{*.base} covariates are treated as baseline confounders. The control
#'   function \code{R.base} is added automatically if not present.
#' @param trtformula A formula for the first–stage treatment model, with left–hand
#'   side equal to \code{tvtrt.base} (e.g. \code{A.base}) and right–hand side
#'   including the IV and baseline covariates (e.g. \code{A.base ~ Z + X1.base + X2.base}).
#'   This model is fit separately within each emulated trial to construct the
#'   control–function residuals.
#' @param data A `seqem` object created by [seqem()], containing the original
#'   stacked pre–censoring data (\code{data_seqorig}) and the artificially
#'   censored, sequential–trial–emulated data (\code{data_seqem}).
#' @param id Character string giving the name of the subject identifier variable.
#' @param trial_id Character string giving the name of the emulated trial
#'   identifier (typically \code{"trial"} as produced by [seqem()]).
#' @param tvtrt Character string giving the name of the time–varying treatment
#'   variable in the original data (e.g. \code{"A"}). The corresponding
#'   trial–baseline term is assumed to be \code{paste0(tvtrt, ".base")}.
#' @param iv Character string giving the name of the baseline instrumental variable.
#' @param stratify Logical; if \code{TRUE}, include \code{strata(trial_id)} in the
#'   second–stage Cox model, allowing a common treatment effect across trials
#'   but trial–specific baseline hazards.
#' @param se_type Character string specifying the standard error estimator:
#'   one of \code{"robust"}, \code{"jackknife"}, \code{"bootstrap"}, or \code{"none"}.
#' @param B Integer; number of bootstrap replicates when \code{se_type = "bootstrap"}.
#' @param seed Optional integer seed for reproducible bootstrap sampling.
#' @param ipcw Logical; must be \code{TRUE}. Artificial censoring weights are
#'   required for valid 2SRI estimation under the sequential trial emulation.
#' @param ipcw_trunc Numeric in (0, 1]; percentile at which to truncate the
#'   cumulative IPACW weights (e.g. \code{0.95} for 95th percentile truncation).
#' @param ipcw_den_covs Optional character vector naming time–varying covariates
#'   whose leads are used in the denominator weight model. By default, derived
#'   from the baseline confounders included in \code{formula}.
#' @param ipcw_num_covs Optional character vector naming baseline covariates
#'   (as \code{*.base}) to use in the numerator weight model. By default, the
#'   baseline confounders used in \code{formula} (excluding the treatment and IV).
#' @param by_trial Logical; if \code{TRUE}, in addition to the composite Cox fit,
#'   fit separate 2SRI–adjusted Cox models within each emulated trial and store
#'   them in \code{$fit_by_trial}.
#' @param ... Additional arguments passed to [survival::coxph()] in the
#'   second–stage fits.
#'
#' @return
#' An object of class \code{"coxiv_seq"} with components including:
#' \itemize{
#'   \item \code{coef}, \code{se}, \code{var}, \code{ci}: composite second–stage
#'         log–hazard ratio estimates, standard errors, covariance matrix, and
#'         95\% confidence intervals.
#'   \item \code{fit}: the fitted composite Cox model object.
#'   \item \code{data_ste}: the final stacked, artificially censored dataset
#'         used in the second–stage analysis, including the control function
#'         \code{R.base} and the weight column indicated by \code{wt_col}.
#'   \item \code{trtfit_by_trial}: list of first–stage treatment model fits
#'         (one per emulated trial) for IV diagnostics.
#'   \item \code{fit_by_trial}: if \code{by_trial = TRUE}, list of trial–specific
#'         second–stage Cox fits.
#'   \item \code{n_ids}, \code{n_trials}, \code{trial_values}: meta–information
#'         on the number of subjects and emulated trials.
#' }
#'
#' @seealso [seqem()], [seqcox()], [survival::coxph()]
#'
#' @export
coxiv_seq <- function(formula,
                      trtformula,
                      data,
                      id,
                      trial_id = "trial",
                      tvtrt,
                      iv,
                      stratify = TRUE,
                      se_type = c("jackknife", "bootstrap", "robust", "none"),
                      B = 200,
                      seed = NULL,
                      ipcw = TRUE,
                      ipcw_trunc = 0.95,
                      ipcw_den_covs = NULL,
                      ipcw_num_covs = NULL,
                      by_trial = FALSE,
                      ...) {

  cl      <- match.call()
  se_type <- match.arg(se_type)

  ## ---------------------- sanity checks ---------------------- ##

  if (!inherits(data, "seqem")) {
    stop("`data` must be a `seqem` object produced by `seqem()`.")
  }
  if (!isTRUE(ipcw)) {
    stop("For sequential 2SRI, artificial censoring weights (IPACW) are required; ",
         "please use `ipcw = TRUE`.")
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a survival formula.")
  }
  if (!inherits(trtformula, "formula")) {
    stop("`trtformula` must be a formula for the first-stage treatment model.")
  }

  # identify LHS of treatment model and check it matches tvtrt.base
  trt_lhs_chr     <- all.vars(trtformula)[1L]
  trt_base_expected <- paste0(tvtrt, ".base")
  if (!identical(trt_lhs_chr, trt_base_expected)) {
    stop("The left-hand side of `trtformula` must be '", trt_base_expected,
         "' corresponding to the trial-baseline treatment status.")
  }

  dat_orig <- data$data_seqorig
  if (!is.data.frame(dat_orig)) {
    stop("`seqem` object must contain `data_seqorig` as a data.frame.")
  }
  if (!all(c(id, trial_id, tvtrt) %in% names(dat_orig))) {
    stop("`id`, `trial_id`, and `tvtrt` must be columns of the sequential data.")
  }

  ## ---------------------- parse outcome formula ---------------------- ##

  tt        <- stats::terms(formula)
  rhs_terms <- attr(tt, "term.labels")

  is_base    <- grepl("\\.base$", rhs_terms)
  base_terms <- rhs_terms[is_base]

  trt_base_var <- paste0(tvtrt, ".base")
  if (!trt_base_var %in% rhs_terms) {
    stop("The RHS of `formula` must contain the treatment baseline term '",
         trt_base_var, "'.")
  }

  iv_base_var <- paste0(iv, ".base")

  # baseline confounders used for IPACW by default:
  # all *.base terms except treatment, IV, and R.base
  base_confounders_base <- setdiff(base_terms,
                                   c(trt_base_var, iv_base_var, "R.base"))

  # helper: strip ".base"
  strip_base <- function(x) sub("\\.base$", "", x)

  ## ---------------------- first-stage 2SRI residuals (R.base) -------- ##

  # create container for control function
  if ("R.base" %in% names(dat_orig)) {
    warning("Overwriting existing 'R.base' column in `data_seqorig`.")
  }
  dat_orig$R.base <- NA_real_

  trial_vec_orig <- dat_orig[, trial_id, drop = TRUE]
  trial_vals     <- sort(unique(trial_vec_orig))
  n_trials       <- length(trial_vals)

  id_vec_orig     <- dat_orig[, id, drop = TRUE]
  ids_unique_orig <- unique(id_vec_orig)
  n_id            <- length(ids_unique_orig)

  trtfit_by_trial <- vector("list", n_trials)
  names(trtfit_by_trial) <- paste0("trtfit_et", seq_len(n_trials))

  for (k in seq_len(n_trials)) {
    tr_val <- trial_vals[k]
    idx_tr <- (trial_vec_orig == tr_val)
    dat_k  <- dat_orig[idx_tr, , drop = FALSE]

    if (nrow(dat_k) == 0L) {
      trtfit_by_trial[[k]] <- NULL
      next
    }

    # baseline row per id for this trial (first row per id)
    id_vec_k   <- dat_k[, id, drop = TRUE]
    base_rows  <- !duplicated(id_vec_k)
    dat_k_base <- dat_k[base_rows, , drop = FALSE]

    # fit first-stage model A.base ~ Z + (baseline covariates)
    trt_fit_k <- stats::glm(
      trtformula,
      data   = dat_k_base,
      family = stats::binomial()
    )

    trtfit_by_trial[[k]] <- trt_fit_k

    # compute residuals: A.base - E[A.base | Z, X]
    A_base_k <- dat_k_base[[trt_lhs_chr]]
    p_hat_k  <- as.numeric(trt_fit_k$fitted.values)
    if (length(p_hat_k) != length(A_base_k)) {
      stop("Length mismatch between fitted values and A.base in first-stage model (trial ",
           tr_val, ").")
    }
    R_k <- A_base_k - p_hat_k

    # map residuals back to all rows for this trial
    ids_base_k <- dat_k_base[, id, drop = TRUE]
    ids_tr     <- dat_k[, id, drop = TRUE]
    m          <- match(ids_tr, ids_base_k)
    if (any(is.na(m))) {
      stop("Unable to match IDs when assigning R.base for trial ", tr_val, ".")
    }
    dat_orig$R.base[idx_tr] <- R_k[m]
  }

  if (any(is.na(dat_orig$R.base))) {
    warning("Some rows have NA in 'R.base' after first-stage estimation.")
  }

  ## ---------------------- IPACW weights & artificial censoring -------- ##

  # After first-stage, create a per-row copy R for use in lead() inside
  # the IPACW helper (R is constant within id × trial, so R.lead1 == R.base
  # except in the last row of each trajectory).
  if ("R.base" %in% names(dat_orig)) {
    dat_orig$R <- dat_orig$R.base
  } else {
    stop("Internal error: 'R.base' not found in `data_seqorig` when constructing IPACW.")
  }

  # Determine covariates to include in the artificial censoring models.
  # Conceptually, AC is another treatment model and should use the same
  # predictors as the first-stage treatment model plus the control
  # function R (and their appropriate lead/baseline versions).
  #
  #  - Denominator: covariates at the next visit (e.g. Z.lead1, L.lead1, R.lead1).
  #  - Numerator  : baseline covariates at trial start (e.g. Z.base, L.base, R.base)
  #                 interacted with visit index via rownum.
  #
  # If the user supplies `ipcw_num_covs` / `ipcw_den_covs`, we honour those
  # but still force inclusion of the IV and R in the appropriate sets.

  # Baseline covariates for the numerator model (as *.base when available)
  if (!is.null(ipcw_num_covs)) {
    baseline_covs_base_ipcw <- unique(c(ipcw_num_covs, "R.base"))
  } else {
    # derive from RHS of trtformula (excluding the LHS A.base)
    trt_rhs_vars <- setdiff(all.vars(trtformula), trt_lhs_chr)
    baseline_covs_base_ipcw <- character(0L)

    for (v in trt_rhs_vars) {
      # Prefer a *.base version if present; otherwise use v itself.
      v_base <- NULL
      if (grepl("\\.base$", v) && v %in% names(dat_orig)) {
        v_base <- v
      } else {
        cand <- paste0(v, ".base")
        if (cand %in% names(dat_orig)) {
          v_base <- cand
        } else if (v %in% names(dat_orig)) {
          v_base <- v
        }
      }
      if (!is.null(v_base)) {
        baseline_covs_base_ipcw <- c(baseline_covs_base_ipcw, v_base)
      }
    }

    # Always include R.base explicitly
    if ("R.base" %in% names(dat_orig)) {
      baseline_covs_base_ipcw <- c(baseline_covs_base_ipcw, "R.base")
    }

    baseline_covs_base_ipcw <- unique(baseline_covs_base_ipcw)
  }

  # Time-varying covariates for the denominator model: names in the original
  # data whose leads will be constructed (e.g. Z -> Z.lead1, L -> L.lead1, R -> R.lead1).
  if (!is.null(ipcw_den_covs)) {
    tv_covs <- unique(c(ipcw_den_covs, iv, "R"))
  } else {
    base_names_for_lead <- unique(strip_base(baseline_covs_base_ipcw))
    tv_covs <- character(0L)

    for (b in base_names_for_lead) {
      if (b %in% names(dat_orig)) {
        tv_covs <- c(tv_covs, b)
      }
    }

    # Ensure IV and R are always available to the denominator model
    if (iv %in% names(dat_orig)) {
      tv_covs <- c(tv_covs, iv)
    }
    if ("R" %in% names(dat_orig)) {
      tv_covs <- c(tv_covs, "R")
    }

    tv_covs <- unique(tv_covs)
  }

  wt_col   <- ".seqcox_wt"
  use_ipcw <- FALSE

  ipcw_obj <- .compute_ipcw_seq(
    dat_seq_orig       = dat_orig,
    id_var             = id,
    trial_var          = trial_id,
    trt_var            = tvtrt,
    trt_base_var       = trt_base_var,
    rownum_var         = "rownum",
    tv_cov_lead        = tv_covs,
    baseline_covs_base = baseline_covs_base_ipcw,
    trunc_prob         = ipcw_trunc
  )

  dat <- ipcw_obj$dat_seq
  if (!"ipw.stab.cum.trunc" %in% names(dat)) {
    stop("IPACW computation did not produce 'ipw.stab.cum.trunc'.")
  }
  dat[[wt_col]] <- as.numeric(dat$ipw.stab.cum.trunc)
  use_ipcw      <- TRUE


  ## ---------------------- build Cox formulas ---------------------- ##

  add_Rbase <- function(formula) {
    f_chr <- paste(deparse(formula), collapse = " ")
    if (grepl("R.base", f_chr, fixed = TRUE)) return(formula)
    parts <- strsplit(f_chr, "~", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop("`formula` must have a left and right hand side.")
    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])
    rhs_new <- if (rhs == "" || rhs == "1") {
      "R.base"
    } else {
      paste(rhs, "R.base", sep = " + ")
    }
    stats::as.formula(paste(lhs, "~", rhs_new), env = environment(formula))
  }

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

  # composite-analysis formula: add R.base, and optionally strata(trial)
  form_fit <- add_Rbase(formula)
  if (isTRUE(stratify)) {
    form_fit <- add_strata(form_fit, trial_id)
  }

  # trial-specific formula: add R.base, but do not include strata(trial)
  form_trial <- add_Rbase(formula)

  ## ---------------------- id & trial meta ---------------------- ##

  id_vec     <- dat[, id, drop = TRUE]
  ids_unique <- unique(id_vec)
  n_id       <- length(ids_unique)

  trial_vec   <- dat[, trial_id, drop = TRUE]
  trial_vals2 <- sort(unique(trial_vec))
  n_trials2   <- length(trial_vals2)

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
    fit_by_trial <- vector("list", n_trials2)
    names(fit_by_trial) <- paste0("fit_et", seq_len(n_trials2))

    for (k in seq_len(n_trials2)) {
      tr_val <- trial_vals2[k]
      dat_k  <- dat[trial_vec == tr_val, , drop = FALSE]

      if (nrow(dat_k) < 1L) {
        fit_by_trial[[k]] <- NULL
        next
      }

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
    call           = cl,
    formula        = form_fit,
    trtformula     = trtformula,
    tvtrt          = tvtrt,
    iv             = iv,
    coef           = coef_hat,
    var            = var_hat,
    se             = se_hat,
    z              = zval,
    p              = pval,
    ci             = ci_mat,
    n              = fit_used$n,
    nevent         = fit_used$nevent,
    fit            = fit_used,
    se_type        = se_type,
    B              = if (se_type == "bootstrap") B else NULL,
    jackknife      = jackknife_est,
    bootstrap      = boot_est,
    ipcw           = use_ipcw,
    ipcw_trunc     = ipcw_trunc,
    # the actual stacked analysis dataset with R.base and weights
    data_ste       = dat,
    wt_col         = wt_col,
    # first-stage info
    trtfit_by_trial = trtfit_by_trial,
    # STE meta
    n_ids          = n_id,
    n_trials       = n_trials2,
    trial_values   = trial_vals2,
    trial_id       = trial_id,
    stratify       = stratify,
    # by-trial second-stage fits
    by_trial       = by_trial,
    fit_by_trial   = fit_by_trial
  )

  class(res) <- "coxiv_seq"
  res
}


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


#' Print method for coxiv_seq objects
#'
#' @export
print.coxiv_seq <- function(x,
                            digits = max(3L, getOption("digits") - 3L),
                            ...) {
  if (!inherits(x, "coxiv_seq")) {
    return(NextMethod())
  }

  cat("Sequential 2SRI Cox IV analysis via sequential trials emulation\n")

  if (!is.null(x$tvtrt)) {
    cat("  Time-varying treatment: ", x$tvtrt, "\n", sep = "")
  }
  if (!is.null(x$iv)) {
    cat("  Baseline IV: ", x$iv, "\n", sep = "")
  }

  if (!is.null(x$n_trials) && !is.null(x$trial_id)) {
    cat(sprintf("Composite second-stage analysis across %d emulated trials (trial variable: '%s')\n",
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

  if (!is.null(x$data_ste)) {
    cat("Final stacked STE dataset with R.base and weights available in $data_ste",
        if (!is.null(x$wt_col)) sprintf(" (weight column: '%s')", x$wt_col) else "",
        "\n", sep = "")
  }

  if (!is.null(x$trtfit_by_trial)) {
    k_avail1 <- sum(!vapply(x$trtfit_by_trial, is.null, logical(1)))
    cat(sprintf("First-stage treatment models: %d fit(s) available in $trtfit_by_trial\n",
                k_avail1))
  }

  if (isTRUE(x$by_trial) && !is.null(x$fit_by_trial)) {
    k_avail2 <- sum(!vapply(x$fit_by_trial, is.null, logical(1)))
    cat(sprintf("Trial-specific second-stage fits: %d model(s) available in $fit_by_trial\n",
                k_avail2))
    cat("  (Access as fit_by_trial$fit_et1, ..., fit_et", length(x$fit_by_trial),
        "; mapping to actual trial values in $trial_values)\n", sep = "")
  }

  cat("\nCoefficients (composite second-stage fit):\n")

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
