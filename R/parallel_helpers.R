#' Temporarily configure the requested multisession worker count
#'
#' The caller must restore the returned plan with on.exit(). The helper
#' also restores the old plan if backend creation itself fails.
#' @importFrom future plan multisession sequential
#' @noRd
.spa_multisession_plan <- function(n_cores) {
  if (!is.numeric(n_cores) || length(n_cores) != 1L ||
      !is.finite(n_cores) || n_cores < 1 || n_cores != floor(n_cores)) {
    stop("n_cores must be a finite positive integer.", call. = FALSE)
  }
  previous <- future::plan("list")
  tryCatch({
    if (n_cores == 1) {
      future::plan(future::sequential)
    } else {
      future::plan(future::multisession, workers = n_cores)
    }
  }, error = function(e) {
    future::plan(previous)
    stop(e)
  })
  previous
}
