#' Orthogonality method-of-moment IV estimation of Cox model
#'
#' Estimate marginal causal treatment effect (in terms of log-hazard ratio) in Cox model setting subject to unmeasured confounding \eqn{U} using instrumental variable based on the IV **O**rthogonality **M**ethod-**o**f-**M**oment estimator in \href{https://link.springer.com/article/10.1007/s10742-014-0117-x}{MacKenzie et al. 2014}.
#'
#'
#' A causal Cox counterfactual survival model for survival time \eqn{T} with additive unmeasured confounding is assumed for treatment \eqn{A}, observed confounder \eqn{X}, and unmeasured confounder \eqn{U}.
#' \deqn{P(T(a)=t|T(a)\ge t, U=u) = \lambda_0(t)\exp(\beta_a a + \beta'x) + h(u,t)}
#' OMOM IV estimator assumes the following conditions: 1) relevance \eqn{Z \not\amalg A | U,X}, 2) exclusion restriction \eqn{Z\amalg (T,D)|A, U,X}, and 3) randomization \eqn{Z\amalg U|X}.
#' The estimates are obtained as the numerical solution to zero-crossing of estimating equation:
#' \deqn{\psi(\beta) \approx \sum^n_{i=1}\int^\tau_0 \left\{ Z_i-\frac{\sum^n_{j=1} Z_j D_j(s)\exp\{\beta_a A_j+\beta'X\}}{\sum^n_{j=1}D_j(s)\exp\{\beta_a A_j+\beta'X\}}dN_i(s) \right\} = 0 }
#' which is identical to the partial Cox score function for instrument \eqn{Z}. \eqn{D(t)} is the event indicator at time t and \eqn{N(t)} is a counting process for survival.
#' Huber-White sandwich estimator is used to obtain variance for \eqn{\hat{\beta}}:
#' \deqn{\text{Cov}(g(\beta)) = \sum^n_{i=1}\int^\tau_0\left\{Z_i-\frac{\sum^n_{j=1} Z_j D_j(s)\exp\{\beta_a A_j+\beta'X\}}{\sum^n_{j=1}D_j(s)\exp\{\beta_a A_j+\beta'X\}}dN_i(s) \right\}dN_i(s)}
#' Suggestive diagnostic test for weak IV was provided by fitting a pseudo-first stage treatment linear model (since OMOM IV is an estimating equation based IV approach, not two-stage based)
#' regressing treatment `trt` on IV `iv` and confounders extracted from all predictors from `formula` input except for `trt`.
#' A nested \link[stats]{anova} F-test was performed to compare fitted treatment models with and without `iv` as predictor as suggested in \href{https://econpapers.repec.org/article/ecmemetrp/v_3a65_3ay_3a1997_3ai_3a3_3ap_3a557-586.htm}{Staiger and Stock, 1997} and \href{https://www.cambridge.org/core/books/abs/identification-and-inference-for-econometric-models/testing-for-weak-instruments-in-linear-iv-regression/8AD94FF2EFD214D05D75EE35015021E4}{Stock and Yogo, 1997}.
#' A popular rule of thumb of indication for weak IV in econometric literature is when F-statistics is less than 10 from comparing IV-exclusion model.
#'
#'
#' @param formula a regression \link[stats]{formula} object with \link[survival]{Surv} object as outcome including survival time and censoring indicator.
#' @param trt a string indicates the name of the treatment variable of interest. Only time-invariant __binary__ or __continuous__ treatment is supported. It is strongly recommended that binary treatment is converted to numeric type (0,1) instead of using factor type. Treatment cannot be multi-categorical.
#' @param iv a string indicates the name of the instrumental variable. Only one time-variant __binary__ or __continuous__ instrumental variable is supported. IV cannot be multi-categorical.
#' @param data a \link[base]{data.frame} object of input data containing all variables specified in `formula`, including observed confounders.
#'
#' @return Returns `coxivomom` and `surviv` object for OMOM IV estimation of Cox model with following components:
#' * `call`: the match call of the OMOM IV model.
#' * `input`: collected input of the OMOM IV model.
#' * `formula`: the outcome regression formula for Cox model
#' * `n`: the effective sample size used for model estimation after missingness removal.
#' * `events`: the effective number of events of survival outcome after missingness removal.
#' * `est_coef` : estimated log-hazard ratios, standard errors, and p-values.
#' * `vcov_mat` : estimated covariance matrix of log-HRs.
#' * `iv_diag`: list of IV diagnosis for weak IV. Include F-statistics and `anova` fit of nested treatment models comparison.
#' * `surv_curve`: a \link[survival]{survfit} object for fitted survival curve using complete data.
#'
#' @docType package
#'
#' @import MASS
#'
#' @name coxiv_omom
#'
#' @examples
#' surv_data = simdat_coxiv(n = 1000, 0.5)
#' coxiv_omom(Surv(time, case) ~ trt + x1, trt = "trt",
#'                            iv = "IV", data = surv_data)
#'
#' @references
#' MacKenzie, T.A., Tosteson, T.D., Morden, N.E. et al. Using instrumental variables to estimate a Cox’s proportional hazards regression subject to additive confounding. *Health Serv Outcomes Res Method* 14, 54–68 (2014). \url{https://doi.org/10.1007/s10742-014-0117-x}
#'
#' @export

coxiv_omom = function(formula, trt, iv, data) {

  # Save input arguments
  call <- match.call()
  input <- as.list(environment())

  # Validate inputs
  if (missing(formula) || missing(trt) || missing(iv) || missing(data)) {
    stop("All arguments 'formula', 'trt', 'iv', and 'data' must be provided.")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  if (iv %in% attr(terms(formula), "term.labels")) {
    stop(paste("Instrumental variable '", iv, "' should not be included in outcome model as predictor.", sep = ""))
  }

  # Extract model frame
  data_outcome_model <- model.frame(formula, data = data)

  # Validate survival data
  if (!inherits(data_outcome_model[, 1], "Surv")) {
    stop("The left-hand side of the formula must be a 'Surv' object.")
  }

  # Extract time, censoring, and covariates
  survtimes <- data_outcome_model[, 1][, 1]
  eventind <- data_outcome_model[, 1][, 2]
  covariates <- data_outcome_model[, setdiff(colnames(data_outcome_model[,-1]), trt),
                                   drop = FALSE]

  # Validate treatment and instrumental variable
  if (!(trt %in% colnames(data))) {
    stop(paste("Treatment variable", trt, "is not found in the data."))
  }
  if (!(iv %in% colnames(data))) {
    stop(paste("Instrumental variable", iv, "is not found in the data."))
  }

  trt_values <- data[[trt]]


  iv_values <- data[[iv]]


  # Handle factors that are binary (two levels only)
  # Process passable treatment/IV cases with factor/string type
  if (is.factor(trt_values) && length(levels(trt_values)) == 2) {
    message("Attention! Binary factor data type detected in treatment.
            Converting to 0/1 based on factor levels.")
    trt_values <- as.numeric(trt_values) - 1
  } else if (is.character(trt_values) && length(unique(trt_values)) == 2) {
    message("Attention! Binary string data type detected in treatment.
            Converting to 0/1 based on factor levels")
    trt_values <- as.numeric(as.factor(trt_values)) -1
  }
  if (is.factor(iv_values) && length(levels(iv_values)) == 2) {
    message("Attention! Binary factor data type detected in instrumental variable.
            Converting to 0/1 based on factor levels.")
    iv_values <- as.numeric(iv_values) - 1
  } else if (is.character(iv_values) && length(unique(iv_values)) == 2) {
    message("Attention! Binary string data type detected in instrumental variable.
            Converting to 0/1 based on factor levels")
    iv_values <- as.numeric(as.factor(iv_values)) -1
  }

  # Ensure treatment and IV at thie point are numeric
  if (!is.numeric(trt_values)) {
    stop("The treatment variable must be continuous or binary (0/1).")
  }
  if (!is.numeric(iv_values)) {
    stop("The instrumental variable must be continuous or binary (0/1).")
  }

  # Check treatment and IV types
  if (all(trt_values %in% c(0, 1))) {
    message("~ Treatment variable is binary.")
  } else if (is.numeric(trt_values)) {
    message("~ Treatment variable is continuous.")
  } else {
    stop("The treatment variable must be binary or continuous.")
  }

  if (all(iv_values %in% c(0, 1))) {
    message("~ IV is binary.")
  } else if (is.numeric(iv_values)) {
    message("~ IV is continuous.")
  } else {
    stop("The treatment variable must be binary or continuous.")
  }

  # Process covariates
  for (col in colnames(covariates)) {
    if (is.character(covariates[[col]])) {
      # Convert character columns to factors
      covariates[[col]] <- as.factor(covariates[[col]])
    } else if (is.factor(covariates[[col]])) {
      # Keep factors as-is
      next
    } else if (!is.numeric(covariates[[col]])) {
      stop(paste("Unsupported data type for covariate:", col))
    }
  }

  # Convert covariates into a design matrix
  # model.matrix automatically handles factors (multi-level or binary)
  cov_matrix <- model.matrix(~ ., data = covariates)[, -1, drop = FALSE]

  # Handle missing data
  complete_data <- data.frame(
    survtimes = survtimes,
    eventind = eventind,
    trt_values = trt_values,
    iv_values = iv_values,
    cov_matrix
  )
  complete_data <- na.omit(complete_data)

  # Reassign variables after removing missing data
  survtimes <- complete_data$survtimes
  eventind <- complete_data$eventind
  trt_values <- complete_data$trt_values
  iv_values <- complete_data$iv_values
  covariates <- complete_data[, -(1:4), drop = FALSE]  # Remaining covariates

  # Reorder data by survival times
  ord <- order(survtimes, decreasing = TRUE)
  survtimes <- survtimes[ord]
  eventind <- eventind[ord]
  trt_values <- trt_values[ord]
  iv_values <- iv_values[ord]
  covariates <- covariates[ord, , drop = FALSE]

  # Combine treatment and covariates
  trt_cov <- as.matrix(cbind(trt_values, covariates))
  W <- as.matrix(cbind(iv_values, covariates))  # Combine IV and covariates

  # Initialize output objects
  S.W1 <- matrix(0, nrow = length(survtimes), ncol = ncol(W))

  # Estimation equation
  Est.Equat <- function(beta) {
    HR <- exp(trt_cov %*% beta)
    S.0 <- cumsum(HR)
    for (i in 1:ncol(W)) {
      S.W1[, i] <- cumsum(W[, i] * HR)
    }
    colSums(eventind * (W - S.W1 / S.0))
  }

  # Solve for beta using root-finding
  result <- tryCatch({
    rootSolve::multiroot(Est.Equat, start = rep(0, ncol(W)))
  }, error = function(e) {
    stop("Root finding algorithm failed: ", e$message)
  })

  # Handle cases where no root is found
  if (any(is.na(result$root))) {
    stop("Root finding algorithm failed: no solution found.")
  }

  beta_hat <- result$root

  # Variance estimation  (MacKenzie et al. OMOM sandwich)
  HR <- exp(trt_cov %*% beta_hat)
  S.0 <- cumsum(HR)

  # Recompute these at beta_hat (do NOT rely on S.W1 left over from multiroot)
  S.W1 <- matrix(0, nrow = length(survtimes), ncol = ncol(W))
  S.X1 <- matrix(0, nrow = length(survtimes), ncol = ncol(trt_cov))

  for (i in 1:ncol(W)) {
    S.W1[, i] <- cumsum(W[, i] * HR)
  }
  for (j in 1:ncol(trt_cov)) {
    S.X1[, j] <- cumsum(trt_cov[, j] * HR)
  }

  # Cross-cumsums needed for the Jacobian ("bread")
  S.W1X1 <- array(0, dim = c(length(survtimes), ncol(W), ncol(trt_cov)))
  for (i in 1:ncol(W)) {
    for (j in 1:ncol(trt_cov)) {
      S.W1X1[, i, j] <- cumsum(W[, i] * trt_cov[, j] * HR)
    }
  }

  # "Meat" and "bread"
  Var.E.E <- Deriv <- matrix(0, nrow = ncol(W), ncol = ncol(W))
  for (i in 1:ncol(W)) {
    for (j in 1:ncol(W)) {
      # A = sum dN * (W - E[W|risk]) (W - E[W|risk])'
      Var.E.E[i, j] <- sum(eventind *
                             (W[, i] - S.W1[, i] / S.0) *
                             (W[, j] - S.W1[, j] / S.0))
      Deriv[i, j] <- -sum(eventind *
                            (S.W1X1[, i, j] / S.0 -
                               (S.W1[, i] * S.X1[, j]) / (S.0^2)))
    }
  }

  Inv.Deriv <- tryCatch(
    solve(Deriv),
    error = function(e) MASS::ginv(Deriv)  # generalized inverse if singular
  )

  Sandwich <- Inv.Deriv %*% Var.E.E %*% t(Inv.Deriv)
  se <- sqrt(diag(Sandwich))

  # Add covariance matrix to output
  rownames(Sandwich) <- colnames(Sandwich) <- c(trt, colnames(cov_matrix))

  # Output results
  result_mat <- data.frame(
    Estimate = beta_hat,
    `Std Error` = se,
    `Z value` = beta_hat / se,
    Pval = 2 * (1 - pnorm(abs(beta_hat / se)))
  )
  colnames(result_mat)[4] = "P(>|z|)"
  rownames(result_mat) <- c(trt, colnames(cov_matrix))

  # 95% Confidence Interval for treatment effect
  upper <- result_mat$Estimate + qnorm(0.975) * result_mat$Std.Error
  lower <- result_mat$Estimate - qnorm(0.975) * result_mat$Std.Error
  conf_int <- data.frame(lower, upper)
  rownames(conf_int) <- rownames(result_mat)
  names(conf_int) = c("2.5%", "97.5%")

  # Instrumental variable relevance diagnostics
  cders <- setdiff(attr(terms(formula), "term.labels"), trt)
  iv_diag <- iv_strgth(covs = cders, trt = trt, iv = iv, data = na.omit(data))
  if (iv_diag$f_stat < 10) message("Attention! Possible weak IV: nested anova has F < 10.0 for IV-exclusion treatment model.\n ")

  # Estimate survival function
  surv_curve <- survfit(Surv(survtimes, eventind) ~ 1, data = complete_data)

  # Organize output results
  result = list(
    call = call,
    input = input,
    formula = formula,
    n = nrow(complete_data),
    events = sum(eventind),
    est_coef = round(result_mat, 12),
    vcov_mat = Sandwich,
    conf_int = conf_int,
    iv_diag = iv_diag,
    surv_curve = surv_curve
  )

  # Define classes
  class(result) <- c("coxivomom", "survivmod")

  result
}

#' @noRd

iv_strgth <- function(covs, trt, iv, data) {
  # Construct saturated and unsaturated treatment model formulas
  formula_sat <- as.formula(paste(trt, "~", paste(c(covs, iv), collapse = " + ")))
  if (is.null(covs)) {
    formula_unsat <- as.formula(paste(trt, "~", 1))
  } else {
    formula_unsat <- as.formula(paste(trt, "~", paste(covs, collapse = " + ")))
  }

  # Fit models
  trtmod_sat <- lm(formula_sat, data = data)
  trtmod_unsat <- lm(formula_unsat, data = data)

  # Perform ANOVA to compute F-statistic
  anova_res <- anova(trtmod_unsat, trtmod_sat)

  # Output IV diagnostics result
  res <- list(f_stat = anova_res$`F`[2],
              nested_anova = anova_res)

  return(res)
}

