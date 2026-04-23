#' vascular: longitudinal vascular surgery follow-up data
#'
#' A synthetic longitudinal dataset designed to mimic repeated follow-up in a
#' vascular surgery study. The data are stored in long start-stop format for
#' survival analyses with a time-varying exposure and a time-varying binary
#' prognostic factor. In particular, `reint` indicates reintervention surgery
#' status over time and `endoleak` indicates time-varying endoleak status over
#' time. Both are coded so that once the status becomes 1 it remains 1 at later
#' visits.
#'
#' @format A data frame with variables:
#' * `id` : Unique patient identifier
#' * `surv.time` : Observed overall follow-up time to death or censoring
#' * `death` : Terminal event indicator
#' * `time.start` : Interval start time
#' * `time.stop` : Interval stop time
#' * `event` : Counting-process event indicator for the interval
#' * `reint` : Time-varying reintervention status at the current visit
#' * `endoleak` : Time-varying endoleak status at the current visit
#' * `iv` : Baseline instrumental variable for reintervention propensity
#' * `smoke` : Baseline smoking status indicator
#' * `age` : Baseline age
#' * `diameter` : Baseline aneurysm diameter
#' * `gender` : Baseline gender indicator
#' * `hyper` : Baseline hypertension indicator
#' * `obese` : Baseline obesity indicator
#' * `reintlag1` : Lagged reintervention status from the previous interval
#' * `visit` : Visit number
#' @source Synthetic data generated to mimic the distributional features of a
#'   vascular surgery longitudinal follow-up study.
#' @usage data(vascular)
"vascular"
