#' Estimate the parameters for simplified Fisher-Gaussian kernel mixture model using EM algorithm
#'
#' Fits a Fisher-Gaussian kernel mixture model to spatial coordinates using the EM algorithm. This models the spatial distribution of spots within each tissue domain
#'
#' @param shinyDesign A \code{shinyDesign} object containing spatial data
#' @param iter_max Maximum number of EM iterations (default = 1000)
#' @param M Number of mixture components(clusters) (default = 5)
#' @param tol Convergence tolerance for EM algorithm (default = 1e-1)
#' return Updated \code{shinyDesign} object with Fisher-Gaussian parameter estimates
#' @import igraph
#' @import dplyr
#' @import invgamma
#' @import pbapply
#' @import MASS
#' @import movMF
#' @import Rfast
#' @import truncnorm
#' @export
#' @examples
#' # Assuming shinyDesign is a valid shinyDesign object
#' # result <- estimation_FGEM(shinyDesign, iter_max = 1000, M =5, tol = 1e-1)

estimation_FGEM <- function(shinyDesign, iter_max = 1000, M_candidates = 2:5, tol = 1e-1){
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
	
    loc_file <- refcolData(shinyDesign)[, c('x','y','domain')]
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]), 
									   xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain

    DOMAIN <- sort(unique(loc_file$domain))
	RST <- lapply(seq_along(DOMAIN), function(d){
		coords_sub <- coords_norm %>% filter(domain == DOMAIN[d])
		FIT <- select_best_M(x = as.matrix(coords_sub[, c('x', 'y')]), M_candidates = M_candidates, iter_max = iter_max, tol = tol)
		return(FIT)	
	})
	   
    message('Completed fitting Fisher-Gaussian mixture models for all domains')
    names(RST) <- DOMAIN
    shinyDesign@paramsFG <- RST
	shinyDesign@selected_M_list_BIC <- sapply(RST, function(fit) fit$best_M_BIC)
	shinyDesign@selected_M_list_AIC <- sapply(RST, function(fit) fit$best_M_AIC)
    return(shinyDesign)
}




