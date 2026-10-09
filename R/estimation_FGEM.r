#' Estimate the parameters for simplified Fisher-Gaussian kernel mixture model via EM algorithm
#'
#' Fits a Fisher-Gaussian kernel mixture model to the spatial coordinates
#' (per domain) using the EM algorithm, and stores results in the \code{spaDesign} object.
#' 
#' @param spaDesign A \code{spaDesign} object containing spatial spatial metadata in \code{refcolData}
#'   with columns \code{"x"}, \code{"y"}, and \code{"domain"} (one row per spot).
#' @param iter_max Maximum number of EM iterations (default = 1000)
#' @param M_candidates Integer vector of candidate mixture sizes to try (default = 2:5)
#' @param tol Positive convergence tolerance for EM algorithm (default = 1e-1)
#' @param n_cores Number of CPU cores for parallization (default = 1). Uses multisession workers on Windows, macOS, and Linux.
#' @param verbose Logical; if \code{TRUE}, print progress messages (default = TRUE).
#' 
#' @return Updated \code{spaDesign} object with Fisher-Gaussian parameter estimates and M list:
#'   \itemize{
#'     \item \code{paramsFG}: list of per-domain fit results (from \code{select_best_M}).
#'     \item \code{selected_M_list_BIC}: named integer vector of best M by BIC (NA if failed).
#'     \item \code{selected_M_list_AIC}: named integer vector of best M by AIC (NA if failed).
#'   }
#'
#' @details
#' Coordinates are rescaled to \([0,1]\times[0,1]\) with \code{igraph::norm_coords} before fitting.
#' For any domain with fewer spots than some candidate M, those candidates are dropped for that domain.
#'
#' @importFrom igraph norm_coords
#' @import dplyr
#' @import future.apply
#' @export
#' @examples
#' \dontrun{
#' # Assuming spaDesign is a valid spaDesign object
#' # result <- estimation_FGEM(spaDesign, 
#'                             iter_max = 1000, 
#'                             M_candidates = 2:5, 
#'                             tol = 1e-1, 
#'                             n_cores = 2, 
#'                             verbose = FALSE)
#' }
#' 
estimation_FGEM <- function(spaDesign, 
                            iter_max = 1000, 
                            M_candidates = 2:5, 
                            tol = 1e-1, 
                            n_cores = 4, 
                            verbose = FALSE){
  
  # Use the requested worker count and restore the caller's plan on exit.
  .spa_previous_plan <- .spa_multisession_plan(n_cores)
  on.exit(future::plan(.spa_previous_plan), add = TRUE)

  # input validation
	if (!is.numeric(iter_max) || iter_max <= 0 || iter_max != round(iter_max)) {
		stop("iter_max must be a positive integer")
	}
	if (!is.numeric(M_candidates) || any(M_candidates <= 0) || any(M_candidates != round(M_candidates))) {
		stop("M_candidates must be positive integers")
	}
	if (!is.numeric(tol) || tol <= 0) {
		stop("tol must be a positive number")
	}
  
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    stop("n_cores must be a positive integer.")
  }
	
  loc_file <- refcolData(spaDesign)[, c('x','y','domain')]
    
  ## scale the coordinates to [0,1] range
  coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]), 
									   xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_file$domain

  DOMAIN <- sort(unique(loc_file$domain))
  
  if(verbose){
    message('Fitting models for ', length(DOMAIN), 'domain(s)')
  }
  
  # fit models for each domain
	RST <- future.apply::future_lapply(seq_along(DOMAIN), function(d){
		          coords_sub <- coords_norm %>% dplyr:: filter(domain == DOMAIN[d])
		          
		          FIT <- tryCatch({
		            select_best_M(
		              x = as.matrix(coords_sub[, c('x', 'y')]), 
		              M_candidates = M_candidates, 
		              iter_max = iter_max, 
		              tol = tol
		              )
		            }, error = function(e){
		               warning('Failed to fit model for domain ', DOMAIN[d], ':', e$messsage)
		              return(NULL)
		              })
		          return(FIT)	
	}, future.seed = TRUE)
	
	if (verbose){
  message('Completed fitting Fisher-Gaussian mixture models for all domains')
	}
	
  names(RST) <- DOMAIN
  spaDesign@paramsFG <- RST
	spaDesign@selected_M_list_BIC <- sapply(RST, function(fit) {
	  if (is.null(fit)) return(NA_integer_)
	  fit$best_M_BIC
	  })
	spaDesign@selected_M_list_AIC <- sapply(RST, function(fit) {
	  if (is.null(fit)) return(NA_integer_)
	  fit$best_M_AIC
	  })
  return(spaDesign)
}




