#' Estimate Parameters for Gaussian Process Model (Parallel Version)
#'
#' Fits a Nearest-Neighbor Gaussian Process (NNGP) model to estimate spatial parameters
#' for each domain-informative gene using parallel processing.
#' @param spadesign A \code{spaDesign} object containing spatial coordinates and expression values
#' @param n_neighbors Number of nearest neighbors for NNGP fitting (default = 10)
#' @param order Ordering scheme for coordinates ('AMMD' or 'Sum_coords', default = 'AMMD')
#' @param verbose Whether to display progress messages (default = FALSE)
#' @param n_cores Number of cores to use for parallel processing (default = parallel::detectCores() - 1)
#' @return Updated \code{spaDesign} object with parameter estimates
#' @import future_apply
#' @import parallel
#' @export

library(future.apply)
library(BRISC)
library(progressr)

estimation_NNGP_parallel <- function(shinyDesign, n_neighbors = 10, ORDER = 'AMMD', X = NULL, verbose = FALSE, workers = 4){
  
  # Prepare data outside loop
  loc_file <- refcolData(shinyDesign)[, c('x','y','domain')]
  count_matrix <- refCounts(shinyDesign)
  
  coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]),
                                     xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_file$domain
  
  DOMAIN <- sort(unique(loc_file$domain))
  topGenes_list <- topGenes(shinyDesign)
  
  if (is.null(topGenes_list) || length(topGenes_list) == 0) {
    stop("topGenes is NULL or empty. Check the shinyDesign object.")
  }
  
  # Setup future
  plan(multisession, workers = workers)
  
  # Parallel over domains
  handlers(global = TRUE)  # allow progress bars
  with_progress({
    p <- progressor(along = DOMAIN)
    
    RST <- future_lapply(seq_along(DOMAIN), function(i) {
      d <- DOMAIN[i]
      idx <- which(coords_norm$domain == d)
      coords_norm_sub <- coords_norm[idx, ]
      
      GENES <- rownames(topGenes_list[[i]])
      if (length(GENES) == 0) return(NULL)
      
      FIT <- lapply(GENES, function(gene) {
        counts.gene <- count_matrix[rownames(count_matrix) == gene, , drop = FALSE]
        counts.sub <- counts.gene[idx]
        test.data <- data.frame(gene = counts.sub,
                                x = coords_norm_sub$x,
                                y = coords_norm_sub$y)
        
        coords.train <- as.matrix(test.data[, c('x','y')])
        
        # Fit GP
        GPrst <- nnGPfit(spatial_coords = coords.train, 
                         logCount = log(test.data$gene + 1),
                         X = X,
                         n_neighbors = n_neighbors,
                         order = ORDER,
                         verbose = verbose)
        return(GPrst)
      })
      names(FIT) <- GENES
      p(sprintf("Domain %s done", d))
      return(FIT)
    })
  })
  
  names(RST) <- DOMAIN
  
  shinyDesign@paramsGP <- RST
  return(shinyDesign)
}


### use foreach instead of future.apply
estimation_NNGP_parallel <- function(spadesign, n_neighbors = 10, order = 'AMMD',verbose = FALSE, n_cores = parallel::detectCores() - 2) {
    
    # Input validation
    if (!order %in% c('AMMD', 'Sum_coords')) {
        stop("order must be either 'AMMD' or 'Sum_coords'")
    }
    if (!is.numeric(n_neighbors) || n_neighbors < 1) {
        stop("n_neighbors must be a positive integer")
    }
    
    loc_file <- refcolData(spadesign)[, c('x','y','domain')]
    count_matrix <- refCounts(spadesign)
    
    # Scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]), 
                                      xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain

    DOMAIN <- sort(unique(loc_file$domain))
    topGenes <- topGenes(spadesign)

    if (is.null(topGenes) || length(topGenes) == 0) {
        stop("topGenes is NULL or empty. Check the spaDesign object.")
    }
    
    # Set up parallel backend
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))  # Ensure cluster is stopped on exit
    
    RST <- foreach::foreach(i = seq_along(topGenes), .packages = c("BRISC", "Matrix"), .export = c('nnGPfit')) %dopar% {
        d <- names(topGenes)[i]  # the domain name
        idx <- which(coords_norm$domain == d)  # indices for current domain
        coords_norm_sub <- coords_norm[idx, ]  # domain coordinates
        
        GENES <- rownames(topGenes[[i]])  # informative genes for that domain
        
        if (length(GENES) == 0) {
            return(NULL)
        }
        
        # Process genes in serial within each parallel domain
        FIT <- lapply(GENES, function(gene){
			counts.gene <- count_matrix[rownames(count_matrix) == gene, , drop = FALSE]
            counts.sub <- counts.gene[idx]  # gene expression within domain
            
            test.data <- data.frame(gene = counts.sub,
                                   x = coords_norm_sub$x,
                                   y = coords_norm_sub$y)
            
            coords.train <- as.matrix(test.data[, c('x','y')])
			nnGPfit(spatial_coords = coords.train, 
                       logCount = log(test.data$gene + 1),
                       X = NULL,  # Replace with your X if needed
                       n_neighbors = n_neighbors,
                       order = order,
                       verbose = verbose)
		})
        names(FIT) <- GENES
        FIT
    }
    
    names(RST) <- names(topGenes)
    spadesign@paramsGP <- RST
    return(spadesign)
}