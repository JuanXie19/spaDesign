#' simualte data with modified effect size and sequencing depth
#'
#' @param spadesign The \code{spadesign} object
#' @param seq_depth_factor A scaling factor for sequencing depth
#' @param effect_size_factor A scaling factor for effect size
#' @param SEED random seed for reproducibility
#' @return A spatialdesign object with simulated count matrix
#' @export
#'
simulateDataEffectSize <- function(spadesign, seq_depth_factor, effect_size_factor, SEED){

    if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) {
        stop("seq_depth_factor must be a positive numeric value.")
    }
    if (!is.numeric(effect_size_factor) || effect_size_factor <= 0) {
        stop("effect_size_factor must be a positive numeric value.")
    }

    count_matrix <- refCounts(spadesign)    
    loc_file <- refcolData(spadesign)[, c('x', 'y', 'domain')]
    par_GP <- paramsPoissonGP(spadesign)
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x', 'y')]),xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain
    
    seqDepth_factor <- seq_depth_factor
    ES_factor <- effect_size_factor

    message("Starting simulation...\n")
    all_count <- pblapply(seq_along(par_GP), function(d){
        domain <- names(par_GP)[d]
        GP.par <- par_GP[[d]]
        
        idx <- which(coords_norm$domain == domain)
        coords_norm_sub <- coords_norm[idx, ]
        coords_norm_out <- coords_norm[-idx, ]
                
        message(sprintf("Simulating genes for domain: %s\n", domain))
                
        gene_count <- lapply (seq_along(GP.par), function(g) {
            tryCatch({
                    gene <- names(GP.par)[g]
                    message(sprintf("Simulating gene: %s in domain: %s\n", gene, domain))
                    alpha.est <- GP.par[[g]][['alpha.est']]
                    rho.est <- GP.par[[g]][['rho.est']]
                    mu.est <- GP.par[[g]][['mu.est']]
                    s.est <- GP.par[[g]][['s.est']]
                    spot_idx <- idx
                
                    tt <- simulate_genecount_ES(SEED, coords_norm_sub, coords_norm_out, alpha.est, rho.est, mu.est, s.est, seqDepth_factor, ES_factor, count_matrix, gene, spot_idx)
                    return(tt)
                    }, error = function(e) {
                    warning(paste('Error in gene', gene, 'of domain', domain, ':', e$message))
                    return(NULL)            
                    })
        })
        
        # Filter out NULL values
        gene_count <- Filter(Negate(is.null), gene_count)
    
        if (length(gene_count) > 0) {
          if (length(gene_count) == 1) {
            gene_count <- t(data.frame(gene_count[[1]]))  # Ensure it's a data frame if only one row
          } else {
            gene_count <- do.call(rbind, gene_count)
          }
        colnames(gene_count) <- colnames(count_matrix)
        rownames(gene_count) <- paste0(domain, "-", names(GP.par))
        } else {
          warning(paste("All gene simulations for domain", domain, "returned NULL."))
          return(NULL)
        }
    
        return(gene_count)
    })
    
    # Filter out NULL values from all_count
    all_count <- Filter(Negate(is.null), all_count)

    if (length(all_count) > 0) {
        all_count <- do.call(rbind, all_count)
        if (is.matrix(all_count) || is.data.frame(all_count)) {
            spadesign@simCounts <- all_count
            spadesign@simcolData <- refcolData(spadesign)
            message("Simulation complete.\n")
            return(spadesign)
        } else {
            stop("Combined all_count is not a matrix or data frame.")
        }
    } else {
        stop("All domain simulations returned NULL.")
    }    
}



#' Simulate expression count for a single gene
#'
#' @param SEED Random seed for reproducibility
#' @param coords_norm_sub Normalized coordinates for spots within the domain
#' @param coords_norm_out Normalized coordinates for spots outside the domain
#' @param alpha_est Estimated alpha parameter
#' @param rho_est Estimated rho parameter
#' @param mu_est Estimated mean parameter
#' @param s_est Estimated s parameter
#' @param seq_depth_factor Sequencing depth factor
#' @param es_factor Effect size factor
#' @param counts Count matrix
#' @param gene Gene name
#' @param spot_idx Indices of spots within the domain
#' @return Simulated counts for the gene
#' @export


simulate_genecount_ES <- function(SEED, coords_norm_sub, coords_norm_out, alpha.est, rho.est, mu.est, s.est, seqDepth_factor, ES_factor, counts, gene, spot_idx){
    
    b.condition <- mean_in(SEED, coords_norm_sub, alpha.est, rho.est, mu.est, s.est)
    
    # the inside mean
    mu1 <- sum(exp(b.condition))
    n1 <- nrow(coords_norm_sub)
    
    # the outside mean
    mu.out <- mean_out(counts, gene, spot_idx)
    n2 <- ncol(counts) - n1
    
    ## calculate lambda_in and lambda_out
    A <- mu1 + n2 * mu.out
    lambda.in <- A * ES_factor * exp(b.condition) / (ES_factor * mu1 + n2 * mu.out)
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


