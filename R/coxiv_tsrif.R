#' Two-stage residual inclusion-frailty (TSRI-F) instrumental variable analysis of Cox model
#'
#' Flexible instrumental analysis of Cox model in a novel **T**wo-**S**tage **R**esidual **I**nclusion with **F**railty (TSRI-F) framework when treatment \eqn{A} and mortality \eqn{D} are subject to unmeasured confounding.
#' Allow estimation of time constant treatment effect (\href{https://doi.org/10.1093/biostatistics/kxx062}{Martinez-Camblor et al. 2019}) and time-varying treatment effect (\href{https://doi.org/10.1111/rssc.12341}{Martinez-Camblor et al. 2019})
#' in terms of log hazard ratio (log-HR). A set of instrumental variables \eqn{Z} and covariates for adjustment \eqn{X} (i.e., measured confounders) are required to estimate treatment effects \eqn{\beta^{(1)}_a} and \eqn{\beta^{(2)}_a} before and after the pre-specified time `tchange`.
#'
#' This function performs two-stage residual inclusion IV analysis of the Cox model with individual frailty
#' when the treatment effect \eqn{\beta_a} is constant over time (`tvareff = FALSE`) or \eqn{\beta_a (t)} changes value at some time _t_ (`tvareff = TRUE`).
#' The IV assumptions are that 1) relevance \eqn{Z \not\amalg A | U,X}, 2) exclusion restriction \eqn{Z\amalg (T,D)|A, U,X}, and 3) randomization \eqn{Z\amalg U|X} with \eqn{U} the unmeasured confounder.
#' TSRI-F posits a first-stage linear model that generates the treatment: \deqn{A = \gamma_0+ \gamma_1Z+\gamma_2U+\gamma'X+\epsilon}
#' and a second-stage outcome model follows Cox proportional hazards form: \deqn{\lambda(t|X,U)=\lambda_0(t)\exp\{\beta_aA+\beta_uU+\beta'X\}}
#' where \eqn{\lambda_0} is the baseline hazard function. The first stage of TSRIF procedure for constant treatment effect estimates _control function_ \eqn{\widehat{R}} from residual of treatment model:
#' \eqn{\widehat{R} = A - \hat{\gamma}_0 -\hat{\gamma}_1Z-\hat{\gamma}'X}, and fit a second-stage Cox model with user-specified individual frailty \eqn{\phi}:
#' \eqn{\widehat{\lambda}(t|A,Z,\widehat{R},X) = \phi\cdot \widehat{\lambda}_0(t)\exp\{\hat{\beta_a}A+\hat{\beta}_r\widehat{R}+\hat{\beta}'X\}}.
#'
#' For more intricate case of time-varying effect, the first stage proceeds identically to compute \eqn{\widehat{R}^{(1)} = A - \hat{\gamma}_0 -\hat{\gamma}_1Z-\hat{\gamma}'X}.
#' Second stage involves fitting two Cox models to estimate \eqn{\beta_a^{(1)} = \beta_a(t), t<} `tchange` and \eqn{\beta_a^{(2)} = \beta_a(t), t\ge} `tchange`.
#' The first-period Cox model is fitted as \eqn{\widehat{\lambda}^1(t|A,Z,\widehat{R}^{(1)},X) = \phi^1 \cdot \widehat{\lambda}_0(t)\exp\{\hat{\beta}^{(1)}_aA+\hat{\beta}^{(1)}_r\widehat{R}^{(1)}+\hat{\beta}'X\}}.
#' Then, update the control function \eqn{\widehat{R}^{(2)} = \widehat{R}^{(1)} + \hat{\beta}^{(1), -1}_r \hat{\phi_1}} for the second period, which are then included as a covariate in the second-period Cox model:
#' \eqn{\widehat{\lambda}^2(t|A,Z,\widehat{R}^{(2)},X) = \phi^2\cdot \widehat{\lambda}_0(t)\exp\{\hat{\beta}_a^{(2)}A+\hat{\beta}^{(2)}_r\widehat{R}^{(2)}+\hat{\beta}'X\}}
#'
#' This approach estimates conditional treatment effect of `trt`. Coefficient estimates of covariates other than `trt` are also estimated and readily available, but these estimates are
#' suggestive of their associational relations with survival outcome, and causal conclusions regarding these covariates with outcome should be cautious.
#'
#' @usage coxiv_tsrif(surv, cens, trt, iv, covs, data, tchange = NULL, tvareff = FALSE,
#'             fdist = "gamma", bootvar = FALSE, B = 50)
#'
#' @param surv a string indicating variable name of survival time \eqn{T}.
#' @param cens a string indicating variable name of event indicator \eqn{D}: 1 for terminal event occurrence and 0 for right censor.
#' @param trt a string indicating the treatment variable name \eqn{A}. Binary or continuous treatment supported.
#' @param iv a string or vector of strings indicating instrumental variable(s) \eqn{Z}. Binary or continuous IV supported.
#' @param covs a string or vector of strings indicating set of measured covariates \eqn{X}.
#' @param data a \link[base]{data.frame} containing all variables needed.
#' @param tchange a positive numeric value specifying the time point at which the treatment effect changes value.
#' @param fdist the frailty distribution, default is `gamma`. Other options include `gaussian` or `t` distribution. Read more about frailty distribution specification at \link[survival]{frailty}.
#' @param bootvar logical, if `TRUE`, bootstrap variance estimation is performed. Otherwise, analytical variance is provided from frailty Cox model.
#' @param B number of bootstrap iterations if `bootvar` is `TRUE`.
#' @param tvareff logical, by default `FALSE` and estimate a time-constant treatment effect model. If `TRUE`, a time-varying treatment effect model is applied.
#'
#' @return Returns `coxivtsrif` and `survivmod` object with following components:
#' * `est_coef` : estimated log-HR of treatment effect. If `tvareff = TRUE`, returns both \eqn{\beta_a^1} and \eqn{\beta_a^2}.
#' * `variances` : estiamted variance of treatment effect estimate. If `bootvar = TRUE`, bootstrap variance is computed.
#' * `ctrl_func` : estimated control function \eqn{\widehat{R}}, or `data.frame` of \eqn{\widehat{R}^1} and \eqn{\widehat{R}^2} if `tvareff = TRUE`.
#' * `trt_model` : a \link[stats]{lm} object fit of the first-stage treatment model.
#' * `out_model` : a \link[survival]{coxph} object fit of the second-stage outcome Cox model. If `tvareff = TRUE`, Cox models for first period `out_model1` and second period `out_model2` are both returned.
#' * `est_boot` :  a vector of length `B` of bootstrap estimates from each iteration. If `bootvar = FALSE`, returns `NULL`.
#' * `data_long` :  a `data.frame` of data for model fit transformed into start-stop form for two-period Cox estimation, if `tvareff = TRUE`.
#' * `surv_curve`: a \link[survival]{survfit} object for fitted survival curve using complete data.
#' * `tvareff` :  logical; mirror the input argument `tvareff`.
#'
#' @examples
#' # load example data
#' data(surviv::PractData)
#'
#' # estimate TSRI-F model with constant effect
#' mod_fit = coxiv_tsrif(formula = Surv(V1, V2) ~ V3 + V4 + V5,
#'                       trtformula = V3 ~ V4 + V5 + V6 + V7,
#'                       tchange = NULL, data = PractData, trt = "V3", iv = c("V6", "V7"),
#'                       fdist = "gaussian", bootvar = F, B = NULL, tvareff = F)
#'
#' # estimate TSRI-F model with time-varying effect and bootstrap variance
#' mod_fit = coxiv_tsrif(formula = Surv(V1, V2) ~ V3 + V4 + V5,
#'                       trtformula = V3 ~ V4 + V5 + V6 + V7,
#'                       tchange = 2, data = PractData, trt = "V3", iv = c("V6", "V7"),
#'                       fdist = "gamma", bootvar = T, B = 100, tvareff = T)
#'
#'
#' @import survival stats
#'
#' @references
#' 1. Pablo Martínez-Camblor, Todd Mackenzie, Douglas O Staiger, Philip P Goodney, A James O’Malley, Adjusting for bias introduced by instrumental variable estimation in the Cox proportional hazards model, *Biostatistics*, Volume 20, Issue 1, January 2019, Pages 80–96,
#'
#' 2. Pablo Martínez-Camblor, Todd A. MacKenzie, Douglas O. Staiger, Phillip P. Goodney, A. James O’Malley, An Instrumental Variable Procedure for Estimating Cox Models with Non-Proportional Hazards in the Presence Of Unmeasured Confounding, *Journal of the Royal Statistical Society Series C: Applied Statistics*, Volume 68, Issue 4, August 2019, Pages 985–1005
#' @export
#'
#'

coxiv_tsrif <- function(surv, cens, trt, iv, covs, data, tchange = NULL, tvareff = FALSE,
                        fdist = "gamma", bootvar = FALSE, B = 50) {
  # Input validation
  if (missing(surv) || missing(cens) || missing(covs) || missing(trt) || missing(iv) || missing(data)) {
    stop("All arguments 'surv', 'covs', 'trt', 'iv', and 'data' must be provided.")
  }
  if (!is.character(surv) || length(surv) != 1) stop("'surv' must be a string for survival time variable")
  if (!is.character(surv) || length(surv) != 1) stop("'cens' must be a string for censoring/mortality variable.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!is.character(trt) || length(trt) != 1) stop("'trt' must be a string for treatment variable.")
  if (!is.null(covs) && !is.character(covs)) stop("'covs' must be a string or a vector of strings for confounder(s).")
  if (!(is.character(iv) && (length(iv) >= 1))) {
    stop("'iv' must be a string or a vector of strings for instrumental variable(s).")
  }
  if (!(surv %in% colnames(data))) stop(paste("Survival time", surv, "not found in the data."))
  if (!(cens %in% colnames(data))) stop(paste("Censoring indicator", cens, "not found in the data."))
  if (!(trt %in% colnames(data))) stop(paste("Treatment variable", trt, "not found in the data."))
  if (!all(iv %in% colnames(data))) stop(paste("Instrumental variable(s)", iv, "not found in the data."))
  if (!all(covs %in% colnames(data))) stop("Some confounders are not found in the data.")
  if (tvareff & (is.null(tchange) || tchange < 0 || !is.numeric(tchange))) stop("A valid change time 'tchange' is needed for time-varying treatment effect estimation.")
  if (!tvareff & !is.null(tchange)) message("~ The supplied 'tchange' value ignored since 'tvareff' is FALSE.")
  if (!bootvar & !is.null(B)) message("~ Bootstrap iterations 'B' ignored since 'bootvar was set to FALSE.")
  if (bootvar && (is.null(B) || !is.numeric(B) || B <= 0)) stop("If bootvar = TRUE, you must specify a positive numeric value for B.")

  # Construct two-stage formulae
  formula <- as.formula(paste("Surv(", surv, ",", cens, ") ~",
                              paste(c(trt, if (!is.null(covs)) covs else NULL), collapse = " + ")))

  trtformula <- as.formula(paste(trt, "~",
                                 paste(c(iv, if (!is.null(covs)) covs else NULL), collapse = " + ")))

  input <- as.list(environment())

  # Handle missing data
  complete_data <- na.omit(data)
  if (nrow(complete_data) < nrow(data)) {
    message("~ Missing data found. Complete case analysis is performed.")
  }

  # Delegate to the appropriate function based on tvareff
  result <- if (tvareff) {
    message("~ A TSRI-F Cox model with *time-varying* treatment effect is estimated.")
    est_tsrif_tvar(formula = formula, trtformula = trtformula, data = complete_data,
                   trt = trt, iv = iv, t = tchange, fdist = fdist, bootvar = bootvar, B = B)
  } else {
    message("~ A TSRI-F Cox model with *constant* treatment effect is estimated.")
    est_tsrif_tfix(formula = formula, trtformula = trtformula, data = complete_data,
                   trt = trt, iv = iv, fdist = fdist, bootvar = bootvar, B = B)
  }

  # Instrumental variable relevance diagnostics
  iv_diag <- iv_strgth(covs = covs, trt = trt, iv = iv, data = complete_data)
  if (iv_diag$f_stat < 10) message("Attention! Possible weak IV: nested anova has F < 10.0 for IV-exclusion treatment model.")

  # Append additional output elements
  result$input <- input
  result$iv_diag <- iv_diag

  # Organize results
  class(result) <- c("coxivtsrif", "survivmod")

  return(result)
}

# TSRIF Cox model for time-constant effect

#' @noRd
est_tsrif_tfix <- function(formula, trtformula, data, trt, iv, fdist = "gamma",
                           bootvar = FALSE, B = NULL) {
  # Ensure the necessary libraries are loaded
  requireNamespace("survival")
  requireNamespace("stats")

  # Extract survival data
  data_outcome_model <- model.frame(formula, data = data)
  if (!inherits(data_outcome_model[,1], "Surv")) stop("The left-hand side of the outcome formula must be a 'Surv' object.")
  survtimes <- data_outcome_model[, 1][, 1]  # Extract survival times
  eventind <- data_outcome_model[, 1][, 2]  # Extract event indicators

  # Estimate survival function
  surv_curve <- survfit(Surv(survtimes, eventind) ~ 1)

  # Extract variables from trtformula
  trt_model <- lm(trtformula, data = data)
  Resid <- residuals(trt_model)

  # Add residuals to data for the second stage
  data$Resid <- Resid

  # Fit the second-stage Cox model with frailty
  cox_model <- coxph(
    formula = as.formula(paste0(
      deparse(formula),
      " + Resid + frailty(as.factor(1:",nrow(data),"), dist = '", fdist, "')"
    )),
    data = data
  )

  # Extract coefficients from outcome model
  est <- coef(cox_model)[trt]
  #names(est)[length(est)] <- "CtrlFunc"

  # Calculate variance
  if (!bootvar) {
    # Use variance from the coxph model
    variance <- vcov(cox_model)[trt, trt]
  } else {
    # Bootstrap variance estimation
    N <- nrow(data)
    boot_coefs <- numeric(B)
    message("... Bootstrap simulation starts ...")

    for (b in 1:B) {
      # Bootstrap resampling
      boot_indices <- sample(1:N, replace = TRUE)
      boot_data <- data[boot_indices, ]

      # First-stage model on bootstrap sample
      boot_trt_model <- lm(trtformula, data = boot_data)
      boot_data$Resid <- residuals(boot_trt_model)

      # Second-stage model on bootstrap sample
      boot_cox_model <- coxph(
        formula = as.formula(paste0(
          deparse(formula),
          " + Resid + frailty(as.factor(1:",nrow(data),"), dist = '", fdist, "')"
        )),
        data = boot_data
      )

      # Store treatment effect estimate
      boot_coefs[b] <- coef(boot_cox_model)[trt]
    }
    message("... Bootstrap simulation ends ...")

    # Variance of bootstrap estimates
    variance <- var(boot_coefs)
  }

  # 95% Confidence Interval for treatment effect
  upper <- est + qnorm(0.975) * sqrt(variance)
  lower <- est - qnorm(0.975) * sqrt(variance)
  conf_int <- c(lower, upper)
  names(conf_int) = c("2.5%", "97.5%")

  # Return results
  return(list(
    est_coef = est,
    variance = variance,
    conf_int = conf_int,
    trt_model = trt_model,
    ctrl_func = Resid,
    out_model = cox_model,
    est_boot = if (bootvar) boot_coefs else NULL,
    surv_curve = surv_curve
  ))
}

# TSRIF Cox model for time-varying effect

#' @noRd
est_tsrif_tvar <- function(formula, trtformula, data, trt, iv, fdist = "gamma",
                           t, bootvar = FALSE, B = NULL) {
  # Ensure the necessary libraries are loaded
  requireNamespace("survival")
  requireNamespace("stats")

  # Extract survival data
  data_outcome_model <- model.frame(formula, data = data)
  if (!inherits(data_outcome_model[,1], "Surv")) stop("The left-hand side of the outcome formula must be a 'Surv' object.")
  survtimes <- data_outcome_model[, 1][, 1]  # Extract survival times
  eventind <- data_outcome_model[, 1][, 2]  # Extract event indicators

  # Estimate survival function
  surv_curve <- survfit(Surv(survtimes, eventind) ~ 1)

  # First-stage model (using the original data)
  trt_model <- lm(trtformula, data = data)
  Resid <- as.vector(residuals(trt_model))
  data$Resid <- Resid

  # Split data for time-varying analysis
  data_long <- survSplit(formula = update(formula, . ~ . + Resid), data = data, cut = t, episode = "timegroup",
                         start = "tstart", end = "tstop", event = "event", id = "new_id")

  # Indices for each time period
  i1 <- which(data_long$timegroup == 1)
  i2 <- which(data_long$timegroup == 2)

  # Fit second-stage Cox models for first time period
  out_model1 <- coxph(
    #formula = update(formula, Surv(tstart, tstop, event) ~ . + Resid + frailty(1:length(i1), dist = fdist)),
    formula = as.formula(paste0(
      "Surv(tstart, tstop, event) ~ ",
      as.character(formula[3]),
      " + Resid + frailty(1:", length(i1), ", dist = '", fdist, "')"
    )),
    data = data_long[i1, ]
  )

  frail1 <- out_model1$frail
  Resid2 <- out_model1$coef["Resid"]^(-1) * frail1 + Resid
  data2 = as.data.frame(cbind(data, Resid2, frail1))
  data2_long= survSplit(formula = update(formula, . ~ . + Resid2 + frail1), data = data2, cut=t, episode ="timegroup",
                        start = "tstart", end = "tstop", event = "event", id = "new_id")

  # Fit the second Cox model for the second period, using updated residual
  out_model2 <- coxph(
    #formula = update(formula, Surv(tstart, tstop, event) ~ . + Resid2 + frailty(1:length(i2), dist = fdist)),
    formula = as.formula(paste0(
      "Surv(tstart, tstop, event) ~ ",
      as.character(formula[3]),
      " + Resid2 + frailty(1:", length(i2), ", dist = '", fdist, "')"
    )),
    data = data2_long[i2, ]
  )

  # Extract coefficients and variances
  est <- c(coef(out_model1)[trt], coef(out_model2)[trt])
  names(est) <- c(paste(trt,".t1",sep = ""), paste(trt,".t2",sep = ""))
  variance <- c(vcov(out_model1)[trt, trt], vcov(out_model2)[trt, trt])
  names(variance) <- c(paste(trt,".t1",sep = ""), paste(trt,".t2",sep = ""))

  if (bootvar) {
    # Bootstrap variance estimation
    N <- nrow(data)
    boot_coefs <- matrix(NA, nrow = B, ncol = 2)
    message("... Bootstrap simulation starts ...")

    for (b in 1:B) {
      # Bootstrap resampling
      boot_indices <- sample(1:N, replace = TRUE)
      boot_data <- data[boot_indices, ]

      # First-stage model on bootstrap sample
      boot_trt_model <- lm(trtformula, data = boot_data)
      boot_Resid <- residuals(boot_trt_model)
      boot_data$Resid <- boot_Resid

      # Split data into time-varying format
      boot_data_long <- survSplit(
        formula = update(formula, . ~ . + Resid), data = boot_data,
        cut = t, episode = "timegroup", start = "tstart", end = "tstop",
        event = "event", id = "new_id")

      # Indices for each time period
      i1 <- which(boot_data_long$timegroup == 1)
      i2 <- which(boot_data_long$timegroup == 2)

      # Fit the first period Cox model
      boot_model1 <- coxph(
        formula = as.formula(paste0(
          "Surv(tstart, tstop, event) ~ ",
          as.character(formula[3]),
          " + Resid + frailty(1:", length(i1), ", dist = '", fdist, "')"
        )),
        data = boot_data_long[i1, ]
      )
      boot_frail1 <- boot_model1$frail
      boot_Resid2 <- boot_model1$coef["Resid"]^(-1) * boot_frail1 + boot_Resid
      boot_data2 <- as.data.frame(cbind(boot_data, Resid2 = boot_Resid2, frail1 = boot_frail1))

      # Update the split data for the second period
      boot_data2_long <- survSplit(
        formula = update(formula, . ~ . + Resid2 + frail1), data = boot_data2,
        cut = t, episode = "timegroup", start = "tstart", end = "tstop",
        event = "event", id = "new_id")

      # Fit the second period Cox model
      boot_model2 <- coxph(
        formula = as.formula(paste0(
          "Surv(tstart, tstop, event) ~ ",
          as.character(formula[3]),
          " + Resid2 + frailty(1:", length(i2), ", dist = '", fdist, "')"
        )),
        data = boot_data2_long[i2, ]
      )

      # Store treatment effect estimates for each period
      boot_coefs[b, ] <- c(coef(boot_model1)[trt], coef(boot_model2)[trt])
    }
    message("... Bootstrap simulation ends ...")

    colnames(boot_coefs) <- c(paste(trt,".t1",sep = ""), paste(trt,".t2",sep = ""))

    # Calculate variances of bootstrap estimates
    variance <- apply(boot_coefs, 2, var)
  }

  # Organize some outputs
  cntrlfunc = as.data.frame(cbind(Resid, Resid2))

  # 95% Confidence Interval for treatment effect
  upper1 <- est[1] + qnorm(0.975) * sqrt(variance[1])
  lower1 <- est[1] - qnorm(0.975) * sqrt(variance[1])
  conf_int1 <- c(lower1, upper1)
  upper2 <- est[2] + qnorm(0.975) * sqrt(variance[2])
  lower2 <- est[2] - qnorm(0.975) * sqrt(variance[2])
  conf_int2 <- c(lower2, upper2)
  conf_int <- matrix(c(conf_int1, conf_int2), nrow = 2, byrow = T)
  colnames(conf_int) = c("2.5%", "97.5%")

  # Return results
  return(list(
    est_coef = est,
    variance = variance,
    conf_int = conf_int,
    trt_model = trt_model,
    ctrl_func = cntrlfunc,
    out_model1 = out_model1,
    out_model2 = out_model2,
    est_boot = if (bootvar) boot_coefs else NULL,
    data_long = data2_long,
    surv_curve = surv_curve
  ))
}


