#' Estimate the NNGP Parameters for domain-informative genes
#'
#' Fits a Nearest-Neighbor Gaussian Process (NNGP) model for each domain-informative gene
#' using the BRISC package.
#'
#' @param spaDesign A \code{spaDesign} object containing spatial coordinates in \code{refcolData}
#' and expression values in \code{refCounts}.
#' @param n_neighbors Number of nearest neighbors for NNGP fitting (default = 10)
#' @param order Ordering scheme for coordinates('AMMD' or 'Sum_coords', default = 'AMMD'). 
#'  use 'Sum_coords' for non-visium data
#' @param X Optional design matrix of covariates
#' @param verbose Logical, whether to display progress messages (default = FALSE)
#' 
#' @return Updated \code{spaDesign} object. The slot \code{paramsGP} is a nested list:
#' \code{paramsGP[[domain]][[gene]]} contains the BRISC NNGP fit for that gene.
#' 
#' @details
#' Coordinates are normalized to [0,1] range before fitting. For each domain, NNGP models
#' are fit to log-transformed expression (log(count + 1)) for all selected genes within that domain.
#' 
#' @return Updated \code{spaDesign} object with \code{paramsGP} slot populated. 
#'   This is a nested list: \code{paramsGP[[domain]][[gene]]} contains the BRISC fit object.
#' 
#' @import igraph
#' @import dplyr
#' @import BRISC
#' @import Matrix
#' @export
#' 
#' @examples
#' \dontrun{
#' toyDATA <- createDesignObject(count_matrix = toyData$toyCount, 
#'                                loc = toyData$loc)
#' toyDATA <- featureSelection(toyDATA, logfc_cutoff = 0.7, 
#'                             mean_in_cutoff = 2, max_num_gene = 50, n_cores = 4)
#' toyDATA <- estimation_NNGP(toyDATA, n_neighbors = 10, order = 'AMMD')
#' }

estimation_NNGP <- function(spaDesign, 
                            n_neighbors = 10, 
                            order = 'AMMD',
                            X = NULL,
                            verbose = FALSE){
    # Input validation
    if (!inherits(spaDesign, "spaDesign")) {
      stop("spaDesign must be a spaDesign object")
    }
  
    if (!is.numeric(n_neighbors) || n_neighbors <= 0 || n_neighbors != round(n_neighbors)) {
      stop("n_neighbors must be a positive integer")
    }
  
    if (!order %in% c('AMMD', 'Sum_coords')) {
      stop("order must be either 'AMMD' or 'Sum_coords'")
   }
  
    if (!is.logical(verbose)) {
      stop("verbose must be logical (TRUE or FALSE)")
    }
    
    # extract data
    loc_file <- refcolData(spaDesign)[, c('x','y','domain')]
    count_matrix <- refCounts(spaDesign)
    
    # Check for required columns
    if (!all(c('x', 'y', 'domain') %in% colnames(spaDesign@refcolData))) {
      stop("refcolData must contain columns: x, y, domain")
    }
    
    # scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(
      as.matrix(loc_file[,c('x','y')]),
      xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain

    DOMAIN <- sort(unique(loc_file$domain))
    topGenes <- topGenes(spaDesign)

    if (is.null(topGenes) || length(topGenes) == 0) {
        stop("topGenes is NULL or empty. Please run featureSelection first.")
    }
    
    if (verbose) {
      message("Fitting NNGP models for ", length(topGenes), " domain(s)")
    }
    
    RST <- lapply(seq_along(topGenes), function(i){
		    d <- names(topGenes)[i] # the domain name
        idx <- which(coords_norm$domain == d) # fit the model for within domain expression
		    coords_norm_sub <- coords_norm[idx, ]  # domain coordinates
		
		    GENES <- rownames(topGenes[[i]])  # the informative genes for that domain
		
		    if (length(GENES) == 0) {
          next
        }
		
		    FIT <- lapply(GENES, function(gene){
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
							               order = order,
							               verbose = verbose)
			  return(GPrst)
		    })		
		    names(FIT) <- GENES
		    return(FIT)
    })
	
    names(RST) <- names(topGenes)
    
    spaDesign@paramsGP <- RST
    return(spaDesign)    
}


#' Fit a Nearest-Neighbor Gaussian Process(NNGP) model
#' 
#' Internal function that fits an NNGP model using the BRISC package
#'
#' @param spatial_coords Matrix of spatial coordinates (n × 2).
#' @param logCount Numeric vector of log-transformed gene expression.
#' @param X Optional design matrix of covariates.
#' @param n_neighbors Integer. Number of nearest neighbors (default: 10).
#' @param order Character. Coordinate ordering: 'AMMD' (default) or 'Sum_coords'.
#' @param verbose Logical. Display progress messages (default: FALSE).
#' 
#' @return BRISC fitted model object.
#' 
#' @import BRISC
#' @noRd

nnGPfit <- function(spatial_coords, logCount, X = NULL, 
					n_neighbors = 10, order = 'AMMD', verbose = FALSE){
  
		      order_brisc <- BRISC::BRISC_order(coords = spatial_coords, 
										                        order = order,
										                        verbose = verbose)
		      
		      nn_brisc <- BRISC::BRISC_neighbor(coords = spatial_coords,
								                      n.neighbors = n_neighbors,
								                      n_omp = 1, 
								                      order = order,
								                      search.type = "tree",
								                      ordering = order_brisc, 
								                      verbose = verbose)
		      
		      BRISCfit <- BRISC::BRISC_estimation(coords = spatial_coords, 
									                      y = logCount,
									                      x = X,
									                      n.neighbors = n_neighbors,
									                      neighbor = nn_brisc,
									                      order = order,
									                      ordering = order_brisc)
		      
	        return(BRISCfit)
}
