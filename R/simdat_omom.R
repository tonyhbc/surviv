#' Simulate toy data for MOM IV estimation of Cox model
#'
#' Simulate survival data with treatment ("trt"), three covariates ("x1", "x2", "x3), and an instrumental variable ("IV").
#' The outcome is survival time ("Tiempo") and censorship ("censor"). Only one covariate is confounder.
#'
#' @param n sample size specified by user.
#' @param effect user-specified causal effect in terms of log(HR) or beta.
#'
#' @examples
#' simdat_omom(0.5)
#' @return a data frame contains survival time, censor indicator, treament, IV, and covariates/
#'
#' @export


simdat_omom = function (n, effect) {
  log.odds.IV.on.trt <- log(5)
  log.odds.Omit.on.trt <- log(3)
  log.HR.trt.on.Time <- 0.5 # that which we wish to estimate
  log.HR.Omitted.on.Time <- log(0.5)
  n <- n
  Censoring.rate <- 0.75
  IV <- rnorm(n)
  Omitted <- rnorm(n)
  Cov.1 <- effect * Omitted + rnorm(n) # Confounder, affects trt and Time,
  Cov.2 <- effect * Omitted + rnorm(n) # Affects trt
  Cov.3 <- effect * Omitted + rnorm(n) # Affects Time only
  trt <- 1/(1+exp(-(log.odds.IV.on.trt*IV + log.odds.Omit.on.trt*Omitted +
                      0.6*Cov.1 + 0.3 * Cov.2 +   0.0*Cov.3))) < runif(n)
  Y <- -log(runif(n)) / exp(log.HR.trt.on.Time*trt + log.HR.Omitted.on.Time*Omitted -
                              0.5*Cov.1 - 0.0*Cov.2   + 0.5*Cov.3)
  C <- runif(n)
  c.rescale <- sort(Y/C)[n*(1-Censoring.rate)]
  C <- C * c.rescale
  Tiempo <- pmin(Y, C)
  censor  <- ifelse(Y <= C, 1, 0)
  Cov <- cbind(Cov.1, Cov.2, Cov.3)
  # generate data for simulation
  data = as.data.frame(cbind(Tiempo, censor, trt, Cov, IV))
  data
}
