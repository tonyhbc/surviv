#' VitD: cohort data on vitamin D and mortality (mirrored)
#'
#' This dataset is mirrored from the \pkg{ivtools} package (version 2.3.0).
#' The upstream documentation notes the data originate from a real cohort study
#' and were modified for public availability.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{`age`}{Age at baseline}
#'   \item{`filaggrin`}{Indicator of filaggrin mutation}
#'   \item{`vitd`}{Serum 25-OH-D (nmol/L) at baseline}
#'   \item{`time`}{Follow-up time}
#'   \item{`death`}{Death indicator during follow-up}
#' }
#' @source Mirrored from \pkg{ivtools} (LGPL (>= 3)).
#' @references Sjolander, Arvid and Martinussen, Torben. "Instrumental Variable Estimation with the R Package ivtools" Epidemiologic Methods, vol. 8, no. 1, 2019, pp. 20180024. https://doi.org/10.1515/em-2018-0024
#' @usage data(VitD)
"VitD"
