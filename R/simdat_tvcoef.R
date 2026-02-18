library(survival)

## Pre-processing steps to generate data for iv analysis
##  In practice, a user could already have a data set ready that they wish to analyze,
##  in which case this step would not be performed.

#' Generate a data set for test purposes or to conduct a simulation study
#'
#' @param n Number of data values
#' @param p Number of exogeneous covariates
#' @param q Number of IVs
#' @param efftrtsurv Treatment effect on survival time
#' @param effcovsurv Observed confounders effect on survival time
#' @param effusurv Unmeasured confounder effect on survival time
#' @param effcovtrt Observed covariates effect on treatment selection
#' @param effivtrt Instrumental variables effect on treatment selection
#' @param effvtrt Unmeasured confounder effect on treatment selection
#' @return Generated data and full info as list
#' @export

simdat_tvcoef <- function(n,p,q,efftrtsurv,effcovsurv,effusurv,effcovtrt,effivtrt,effvtrt)	{
  cov=rnorm(n*p,0,1)    # measured confounding variable (can generalize to a matrix containing values for p covariates)
  if (p>1) {
    cov=matrix(cov,n,p) #establish multiple covariates
  }
  u<- rgamma(n,1,1)-1 # unmeasured confounding variable within the survival model
  v=u                 # unmeasured confounding variable within the assignment treatment model (equal to u in this case)
  iv=rnorm(n*q,0,1)     # instrumental variable (can generalize to a matrix containing values for q covariates)
  if (q>1) {
    iv=matrix(iv,n,q) #establish multiple covariates
  }

  #Generate treatment values
  coveff=cov
  if (p>1) {
    coveff = cov %*% effcovtrt
  }
  iveff=iv
  if (q>1) {
    iveff = iv %*% effivtrt
  }
  x= iveff + coveff + effvtrt*v + rnorm(n,0,1) # model assignment (continuous in this case.
  # Include an additional unmeasured confounding: rnorm(n,0,1)
  #Generate survival times
  coveff=cov
  if (p>1) {
    coveff = cov %*% effcovsurv # can generalize to effects than (1,1)
  }
  lt=exp(coveff + efftrtsurv*(x-mean(x))+ effusurv*u)     # risk model. Covariates centered in order to obtain realistic times
  y=(rexp(n)/lt)^(0.5)             # Times

  #Generate censoring times:
  censtimes=(rexp(n)/0.02)^(0.5)
  eventind=ifelse(y<=censtimes,1,0)
  survtimes=y*eventind+censtimes*(1-eventind)

  #Form data including unmeasured covariate
  data=as.data.frame(cbind(survtimes,eventind,x,cov,iv,u)) # data
  FCM<- coxph(Surv(survtimes,eventind)~ x + coveff + u) # Real model

  #Form data frame: this is what a user would instead input
  data=as.data.frame(cbind(survtimes,eventind,x,cov,iv))
  return(list(data=data,FullInfo=FCM))
}
