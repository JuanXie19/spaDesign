#' Simualte spatial transcriptomics data with modified effect size and sequencing depth
#'
#' @param shinyDesign A \code{shinyDesign} object containing spatial coordiantes, gene expression and estimated parameters
#' @param seq_depth_factor Numeric scaling factor for sequencing depth (must be > 0 )
#' @param effect_size_factor Numeric scaling factor for effect size (must be > 0)
#' @param SEED random seed for reproducibility (integer)
#' @return A \code{spatialdesign} object with simulated count matrix in the \code{simCounts} slot
#' @import future.apply
#' @export
#'
#' @examples
#' \dontrun{
#' simulated_data <- simulation_EffectSize(spadesign_obj, 
#'                                        seq_depth_factor = 1.5,
#'                                        effect_size_factor = 2.0,
#'                                        SEED = 123)
#' }

library(future.apply)
library(pbapply)  # optional for progress bar if you want

simulation_EffectSize <- function(shinyDesign, seq_depth_factor, effect_size_factor, SEED, workers = 4) {
  
  # Validate inputs as before
  if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) stop("seq_depth_factor must be positive numeric.")
  if (!is.numeric(effect_size_factor) || effect_size_factor <= 0) stop("effect_size_factor must be positive numeric.")
  
  # Extract only the data you need outside future_lapply
  count_matrix <- refCounts(shinyDesign)
  loc_data <- refcolData(shinyDesign)[, c('x', 'y', 'domain')]
  par_GP <- paramsGP(shinyDesign)
  
  # Normalize coords once
  coords_norm <- igraph::norm_coords(as.matrix(loc_data[, c('x', 'y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_data$domain
  
  # Setup parallel plan and increase globals max size
  options(future.globals.maxSize = 2 * 1024^3) # 6 GiB just in case
  plan(multisession, workers = workers)
  
  # Helper function for each domain; pass only what is needed
  simulate_domain <- function(domain_name, par_GP_domain, coords_norm_df, count_mat, seqDepth_factor, ES_factor, SEED) {
    
    idx <- which(coords_norm_df$domain == domain_name)
    coords_norm_sub <- coords_norm_df[idx, ]
    coords_norm_out <- coords_norm_df[-idx, ]
    
    gene_counts_list <- lapply(names(par_GP_domain), function(gene) {
      nnGP_fit <- par_GP_domain[[gene]]
      
      # Call the lower-level simulation for a single gene
      simulate_genecount_ES(
        SEED = SEED,
        coords_norm_sub = coords_norm_sub,
        coords_norm_out = coords_norm_out,
        nnGP_fit = nnGP_fit,
        seqDepth_factor = seqDepth_factor,
        ES_factor = ES_factor,
        counts = count_mat,
        gene = gene,
        spot_idx = idx
      )
    })
    
    # Combine all gene counts for this domain into matrix/data.frame
    gene_counts_mat <- do.call(rbind, gene_counts_list)
    colnames(gene_counts_mat) <- colnames(count_mat)
    rownames(gene_counts_mat) <- paste0(domain_name, "-", names(par_GP_domain))
    
    return(gene_counts_mat)
  }
  
  # Run parallel over domains explicitly passing minimal globals
  all_counts <- future_lapply(
    names(par_GP),
    function(domain_name) {
      simulate_domain(
        domain_name = domain_name,
        par_GP_domain = par_GP[[domain_name]],
        coords_norm_df = coords_norm,
        count_mat = count_matrix,
        seqDepth_factor = seq_depth_factor,
        ES_factor = effect_size_factor,
        SEED = SEED
      )
    }
  )
  
  # Combine all domain results
  all_counts <- do.call(rbind, all_counts)
  
  if (is.matrix(all_counts) || is.data.frame(all_counts)) {
    shinyDesign@simCounts <- all_counts
    shinyDesign@simcolData <- refcolData(shinyDesign)
    message("Simulation complete.\n")
    return(shinyDesign)
  } else {
    stop("Combined all_count is not a matrix or data frame.")
  }
}


#' Simulate expression count for a single gene
#'
#' @param SEED Random seed for reproducibility
#' @param coords_norm_sub Normalized coordinates for spots within the domain
#' @param coords_norm_out Normalized coordinates for spots outside the domain
#' @param nnGP_fit ouutput from parameter estimation step
#' @param seq_depth_factor Sequencing depth factor
#' @param es_factor Effect size factor
#' @param counts Count matrix
#' @param gene Gene name
#' @param spot_idx Indices of spots within the domain
#' @return Simulated counts for the gene
#' @export


simulate_genecount_ES <- function(SEED, coords_norm_sub, coords_norm_out, nnGP_fit, seqDepth_factor, ES_factor, counts, gene, spot_idx){
    
    b.condition <- mean_in(SEED, coords_norm_sub, nnGP_fit)
    
    # the inside mean
    mu1 <- sum(b.condition)
    n1 <- nrow(coords_norm_sub)
    
    # the outside mean
    mu.out <- mean_out(counts, gene, spot_idx)
    n2 <- ncol(counts) - n1
    
    ## calculate lambda_in and lambda_out
    A <- mu1 + n2 * mu.out
    lambda.in <- A * ES_factor * b.condition / (ES_factor * mu1 + n2 * mu.out)
    lambda.out <- mu.out * A / (ES_factor * mu1 + n2 * mu.out)
    
    set.seed(SEED)    
    y.post <- rpois(n = nrow(coords_norm_sub), lambda = seqDepth_factor * lambda.in)
	
    sim.in <- data.frame(x = coords_norm_sub$x, y = coords_norm_sub$y, sim = y.post)
    rownames(sim.in) <- rownames(coords_norm_sub)
    
    set.seed(SEED)
    y.out <- rpois(n = n2, lambda = seqDepth_factor * lambda.out)
    sim.out <- data.frame(x = coords_norm_out$x, y = coords_norm_out$y, sim = y.out)
    rownames(sim.out) <- rownames(coords_norm_out)
    
    count.sim <- rbind(sim.in, sim.out)
    count.sim <- count.sim[order(match(rownames(count.sim),colnames(counts))),]
    identical(rownames(count.sim),colnames(counts))
    tt <- count.sim$sim
    return(tt) 
}

#' This function simulates data based on fitted NNGP model
#' 
#' @param SEED Random seed for reproducibility
#' @param coords_norm_sub Normalized coordinates for spots within the domain
#' @param nnGP_fit ouutput from parameter estimation step 
#' @return A numeric vector representing the sampled expression levels within the domain
#' @export
mean_in <- function(SEED, coords_norm_sub, nnGP_fit) {
    if (missing(SEED) || !is.numeric(SEED)) stop("SEED must be a numeric value.")
    if (!is.data.frame(coords_norm_sub) || !all(c("x", "y") %in% names(coords_norm_sub))) stop("coords_norm_sub must be a data frame with 'x' and 'y' columns.")
	
	coords_spatial <- as.matrix(coords_norm_sub[, c('x', 'y')])
	brisc_pred <- BRISC::BRISC_prediction(nnGP_fit, coords_spatial)
    sigma2 <- nnGP_fit$Theta['sigma.sq']
	tau2 <- nnGP_fit$Theta['tau.sq']   
    b.condition <- round(exp(brisc_pred$prediction + tau2/2))-1
    return(b.condition)
}

#' Calculate the Mean for Outside Domain Spots
#'
#' This function calculates the mean of the lower half of counts outside the domain for a given gene.
#' 
#' @param counts A matrix where rows represent genes and columns represent spots.
#' @param gene A character string representing the gene of interest.
#' @param spot_idx Numeric vector indicating the indices of spots within the domain.
#' @return A numeric value representing the mean of the lower half of counts outside the domain for the given gene.
#' @export
mean_out <- function(counts, gene, spot_idx) {
    if (!inherits(counts, "matrix") && !inherits(counts, "CsparseMatrix")) {
      stop("counts must be a matrix or dgCMatrix.")
    }
    if (!is.character(gene) || length(gene) != 1) stop("gene must be a single character string.")
    if (!is.numeric(spot_idx)) stop("spot_idx must be numeric.")
    
    gene.count <- counts[which(rownames(counts) == gene), ]
    OUT <- gene.count[-spot_idx]
    OUT.lower <- OUT[OUT <= median(OUT)]
    mean.out.low <- mean(OUT.lower)
    return(mean.out.low)
}


