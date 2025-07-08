#' Estimate Parameters for Gaussian Process Model
#'
#' Fits a Nearest-Neighbor Gaussian Process (NNGP) model to estimate spatial parameters
#' for each domain-informative gene.
#' @param shinyDesign A \code{shinyDesign} object containing spatial coordiantes and expression values for selected domain-informative genes
#' @param n_neighbors Number of nearest neighbors for NNGP fitting (default = 10)
#' @param order Ordering scheme for coordinates('AMMD' or 'Sum_coords', default = 'AMMD')
#' @param verbose Whether to display progress messages (default = FALSE)
#' @return Updated \code{shinyDesign} object with parameter estimates
#' @import igraph
#' @import dplyr
#' @import BRISC
#' @import Matrix
#' @export
#' @examples
#' ## Example usage
#' toyDATA <- createshinyDesignObject(count_matrix = toyData$toyCount, loc = toyData$loc)
#' toyDATA <- estimation_NNGP(toyDATA, n_neighbors = 10, order = 'AMMD', verbose = FALSE)
#'

estimation_NNGP <- function(shinyDesign, n_neighbors = 10, ORDER = 'AMMD',X = NULL, verbose = FALSE){

    loc_file <- refcolData(shinyDesign)[, c('x','y','domain')]
    count_matrix <- refCounts(shinyDesign)
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]),xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain

    DOMAIN <- sort(unique(loc_file$domain))
    topGenes <- topGenes(shinyDesign)

    if (is.null(topGenes) || length(topGenes) == 0) {
        stop("topGenes is NULL or empty. Check the shinyDesign object.")
    }
    
    RST <- list()
    for (i in seq_along(topGenes)){
        d <- names(topGenes)[i] # the domain name
        idx <- which(coords_norm$domain == d) # fit the model for within domain expression
		coords_norm_sub <- coords_norm[idx, ]  # domain coordinates
		
		GENES <- rownames(topGenes[[i]])  # the informative genes for that domain
		
		if (length(GENES) == 0) {
          next
        }
		
		FIT <- list()
		for(g in 1:length(GENES)){
			gene <- GENES[g]
			counts.gene <- count_matrix[rownames(count_matrix) == gene, , drop = FALSE]
			counts.sub <- counts.gene[idx] # informative gene within-domain gene expression
			test.data <- data.frame(gene = counts.sub,
								x = coords_norm_sub$x,
								y = coords_norm_sub$y)
		
			coords.train <- as.matrix(test.data[, c('x','y')])
			## fit the NNGP on the log transformed data using BRISC
			GPrst <- nnGPfit(spatial_coords = coords.train, 
							logCount = log(test.data$gene + 1),
							X = X,
							n_neighbors = n_neighbors,
							order = ORDER,
							verbose = verbose)
			FIT[[g]] <- GPrst
		}
		names(FIT) <- GENES
		RST[[i]] <- FIT
    }
    names(RST) <- names(topGenes)
    
    shinyDesign@paramsGP <- RST
    return(shinyDesign)    
}


#' Fit a Nearest-Neighbor Gaussian Process(NNGP) model
#' The internal function that fits an NNGP model using the BRISC package
#'
#' @param spatial_coords Matrix of spatial coordinates
#' @param logCount vector of log-transformed gene expression counts
#' @param X Optional design matrix of covariates
#' @param n_neighbors Number of nearest neighbors for fitting NNGP model with BRISC. Default = 10.
#' @param order: Ordering scheme to use for ordering coordinates with BRISC.('AMMD' or 'Sum_coords', default = 'AMMD'). See BRISC documentation for details.
#' @param verbose Whether to display process messages (default = FALSE)
#' @return BRISC fitted model objet
#' @import BRISC
#' @keywords internal

nnGPfit <- function(spatial_coords, logCount, X = NULL, 
					n_neighbors = 10, order = 'AMMD', verbose = FALSE){					
		order_brisc <- BRISC::BRISC_order(coords = spatial_coords, 
										  order = order,
										  verbose = verbose)										  
		nn_brisc <- BRISC_neighbor(coords = spatial_coords,
								   n.neighbors = n_neighbors,
								   n_omp = 1, 
								   order = order,
								   search.type = "tree",
								   ordering = order_brisc, 
								   verbose = verbose)
		BRISCfit <- BRISC_estimation(coords = spatial_coords, 
									 y = logCount,
									 x = X,
									 n.neighbors = n_neighbors,
									 neighbor = nn_brisc,
									 order = order,
									 ordering = order_brisc)
	   return(BRISCfit)
	   }