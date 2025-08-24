#' Estimate the parameters for simplified Fisher-Gaussian kernel mixture model via EM algorithm
#'
#' Fits a Fisher-Gaussian kernel mixture model to the spatial coordinates
#' (per domain) using the EM algorithm, and stores results in the \code{shinyDesign} object.
#' 
#' @param shinyDesign A \code{shinyDesign} object containing spatial spatial metadata in \code{refcolData}
#'   with columns \code{"x"}, \code{"y"}, and \code{"domain"} (one row per spot).
#' @param iter_max Maximum number of EM iterations (default = 1000)
#' @param M_candiate Integer vector of candidate mixture sizes to try (default = 2:5)
#' @param tol Positive convergence tolerance for EM algorithm (default = 1e-1)
#' @param n_cores Number of CPU cores for parallization (default = 1). Uses \code{pbmclapply}; 
#' on windows this is treated as sequential.
#' @param verbose Logical; if \code{TRUE}, print progress messages (default = TRUE).
#' 
#' return Updated \code{shinyDesign} object with Fisher-Gaussian parameter estimates and M list:
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
#' @import invgamma
#' @import pbmclapply
#' @import MASS
#' @import movMF
#' @import Rfast
#' @import truncnorm
#' @export
#' @examples
#' # Assuming shinyDesign is a valid shinyDesign object
#' # result <- estimation_FGEM(shinyDesign, iter_max = 1000, M_candidates = 2:5, tol = 1e-1, n_cores = 2, verbose = FALSE)

estimation_FGEM <- function(shinyDesign, iter_max = 1000, M_candidates = 2:5, tol = 1e-1, n_cores = 2, verbose = FALSE){
	message("DEBUG: Entering estimation_FGEM with iter_max=", iter_max)
	
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
	
  loc_file <- refcolData(shinyDesign)[, c('x','y','domain')]
    
  ## scale the coordinates to [0,1] range
  coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]), 
									   xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_file$domain

  DOMAIN <- sort(unique(loc_file$domain))
	RST <- pbmclapply(seq_along(DOMAIN), function(d){
		          coords_sub <- coords_norm %>% dplyr:: filter(domain == DOMAIN[d])
		          FIT <- select_best_M(x = as.matrix(coords_sub[, c('x', 'y')]), M_candidates = M_candidates, iter_max = iter_max, tol = tol)
		          return(FIT)	
	}, mc.cores = n_cores)
	   
  message('Completed fitting Fisher-Gaussian mixture models for all domains')
  names(RST) <- DOMAIN
  shinyDesign@paramsFG <- RST
	shinyDesign@selected_M_list_BIC <- sapply(RST, function(fit) fit$best_M_BIC)
	shinyDesign@selected_M_list_AIC <- sapply(RST, function(fit) fit$best_M_AIC)
  return(shinyDesign)
}




