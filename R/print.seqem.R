#' Print method for seqem objects
#'
#' @param x An object of class `seqem` returned by [seqem()].
#' @param ... Ignored.
#'
#' @export
print.seqem <- function(x, ...) {
  if (!inherits(x, "seqem")) {
    return(NextMethod())
  }
  if (!is.null(x$data_seqem) && is.data.frame(x$data_seqem)) {
    cat(.seqem_summary_text(x$data_seqem))
  } else if (!is.null(x$summary) && is.character(x$summary)) {
    cat(x$summary)
  } else {
    cat("Object of class 'seqem' with no summary available.\n")
  }
  invisible(x)
}
