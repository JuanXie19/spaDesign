#' Estimate the parameters for Fisher-Gaussian kernel mixture model
#' @param spadesign A \code{spaDesign} object
#' @param n_iter number of MCMC iterations
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
#' # Assuming spadesign is a valid spaDesign object
#' # result <- estimateParamsFG(spadesign, n_iter = 2000)

estimateParamsFG <- function(spadesign, n_iter = 5000){
	
    loc_file <- refcolData(spadesign)[, c('x','y','domain')]
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain

    DOMAIN <- sort(unique(loc_file$domain))
    message(paste0('Fitting Fisher-Gaussian mixture models for ', length(DOMAIN), ' domains'))
    
    results <- pblapply(DOMAIN, function(d){
        message(sprintf('Fitting Fisher-Gaussian mixture model for domain: %s', d))
        coords_domain <- coords_norm %>% filter (domain == d)
        FG_fit <- fgFit_single(coords_domain[, c("x", "y")], n_iter = n_iter)
        message(sprintf('Completed Fisher-Gaussian mixture model for domain: %s', d))
        return(FG_fit)
    })
    
    message('Completed fitting Fisher-Gaussian mixture models for all domains')
    names(results) <- DOMAIN
    spadesign@paramsFG <- results
    return(spadesign)
}

#' Fit Fisher-Gaussian kernel mixture for a single tissue domain
#' @title Fit Fisher-Gaussian kernel mixture model for the spatial coordinates of a given tissue domain
#' @description The Fisher-Gaussian kernel mixture model was proposed by M Mukhopadhyay et al., J R Stat Soc Series B Stat Methodol. 2020. 
#' The model was originally developed for density estimation. Here, it is used to model the spatial distribution of spots from a tissue domain.
#' @param loc data frame. The normalized 2-d spatial coordinates of a spatial domain
#' @param n_iter number of MCMC iterations
#' @return The fitted Fisher-Gaussian kernel mixture model
#' @examples
#' # Assuming loc is a data frame with normalized coordinates
#' result <- fgFit_single(loc, n_iter = 2000)
#' @export
#'
 fgFit_single <- function(loc, n_iter){
    
    if (!is.data.frame(loc)) stop("loc must be a data frame.")
    if (!all(c("x", "y") %in% names(loc))) stop("loc must contain 'x' and 'y' columns for coordinates.")
    if (!is.numeric(n_iter) || n_iter <= 0) stop("n_iter must be a positive integer.")

    message(sprintf('Fitting Fisher-Gaussian mixture model for %d spots', nrow(loc)))
    message('This process may take a while; please be patient.')
	
    FG_fit <- FG_mixture(loc, Iter = n_iter)
	
    message('Completed fitting Fisher-Gaussian mixture model!')
    return(FG_fit)    
}
