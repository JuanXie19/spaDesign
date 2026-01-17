#' Simulate spatial counts with scaled effect size and sequencing depth
#'
#' Simulates spatial gene expression counts using fitted NNGP models from a \code{spaDesign} object.
#' For each spatial domain, the function generates Poisson counts for domain-informative genes by:
#' (i) predicting within-domain mean expression from the fitted NNGP model, and
#' (ii) estimating outside-domain baseline expression from observed counts.
#'
#' Two global scaling factors are supported:
#' \itemize{
#'   \item \code{seq_depth_factor}: controls overall sequencing depth (global count scaling)
#'   \item \code{effect_size_factor}: controls within-domain enrichment strength
#' }
#' 
#' Simulated counts are stored in the \code{simCounts} slot and spot metadata is copied
#' into \code{simcolData}.
#'
#' @param spaDesign A \code{spaDesign} object containing:
#'   \itemize{
#'     \item Reference count matrix (accessible via \code{refCounts})
#'     \item Spatial coordinates and domain labels (accessible via \code{refcolData})
#'     \item Fitted NNGP parameters for domain-informative genes (accessible via \code{paramsGP})
#'   }
#' @param seq_depth_factor Numeric scalar > 0. Scaling factor for sequencing depth.
#'   Values greater than 1 increase the expected counts; values less than 1 decrease counts.
#' @param effect_size_factor Numeric scalar > 0. Scaling factor for within-domain effect size.
#'   Values greater than 1 strengthen domain enrichment; values less than 1 weaken it.
#' @param SEED Integer random seed for reproducibility.
#' @return A \code{spaDesign} object with:
#'   \itemize{
#'     \item \code{simCounts}: simulated gene-by-spot count matrix
#'     \item \code{simcolData}: spot metadata copied from \code{refcolData(spaDesign)}
#'   }
#'
#' @import igraph
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate counts with increased depth and stronger domain effect
#' sim_obj <- simulation_EffectSize(
#'   spaDesign = my_spadesign,
#'   seq_depth_factor = 1.5,
#'   effect_size_factor = 2.0,
#'   SEED = 123
#' )
#'
#' sim_counts <- sim_obj@simCounts
#' dim(sim_counts)
#' }

simulation_EffectSize <- function(spaDesign, seq_depth_factor, effect_size_factor, SEED){

    if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) {
        stop("seq_depth_factor must be a positive numeric value.")
    }
    if (!is.numeric(effect_size_factor) || effect_size_factor <= 0) {
        stop("effect_size_factor must be a positive numeric value.")
    }
	
	# extract reference data
    count_matrix <- refCounts(spaDesign)    
    loc_data <- refcolData(spaDesign)[, c('x', 'y', 'domain')]
    par_GP <- paramsGP(spaDesign)
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_data[,c('x', 'y')]),xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_data$domain
    
    seqDepth_factor <- seq_depth_factor
    ES_factor <- effect_size_factor

    message("Starting simulation...\n")
    all_count <- lapply(seq_along(par_GP), function(d){
        domain <- names(par_GP)[d]
        GP.par <- par_GP[[d]]
        
        idx <- which(coords_norm$domain == domain)
        coords_norm_sub <- coords_norm[idx, ]
        coords_norm_out <- coords_norm[-idx, ]
                
        message(sprintf("Simulating genes for domain: %s\n", domain))
                
        gene_count <- lapply(seq_along(GP.par), function(g) {
            tryCatch({
                    gene <- names(GP.par)[g]
                    message(sprintf("Simulating gene: %s in domain: %s\n", gene, domain))
                    nnGP_fit <- GP.par[[g]]
                    
                    spot_idx <- idx
                
                    tt <- simulate_genecount_ES(SEED, coords_norm_sub, coords_norm_out, nnGP_fit, seqDepth_factor, ES_factor, count_matrix, gene, spot_idx)
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
            spaDesign@simCounts <- all_count
            spaDesign@simcolData <- refcolData(spaDesign)
            message("Simulation complete.\n")
            return(spaDesign)
        } else {
            stop("Combined all_count is not a matrix or data frame.")
        }
    } else {
        stop("All domain simulations returned NULL.")
    }    
}


#' (Refactored)Simulate spatial transcriptomics data with modified effect size and sequencing depth
#'
#' This function generates simulated spatial gene expression counts based on a fitted NNGP model.
#' It allows scaling of sequencing depth and effect size for domain-informative genes.
#'
#' @param spaDesign A \code{spaDesign} object containing spatial coordinates, gene expression, 
#'   and estimated NNGP parameters.
#' @param seq_depth_factor Numeric scaling factor for sequencing depth (must be > 0).
#' @param effect_size_factor Numeric scaling factor for effect size (must be > 0).
#' @param SEED Integer random seed for reproducibility.
#' @return A \code{spaDesign} object with simulated count matrix stored in \code{simCounts}, 
#'   and updated \code{simcolData} with spot metadata.
#'
#'
#'
#'
#'

#' Simulate Spatial Counts with Scaled Sequencing Depth and Effect Size (Refactored)
#'
#' A refactored version of \code{\link{simulation_EffectSize}} with reduced memory usage.
#' Instead of passing the full reference count matrix into each gene simulation step,
#' this function extracts one gene vector at a time and passes only the required data.
#'
#' The simulation results and output format are consistent with \code{\link{simulation_EffectSize}}.
#'
#' @param spaDesign A \code{spaDesign} object containing:
#'   \itemize{
#'     \item Reference count matrix (accessible via \code{refCounts})
#'     \item Spatial coordinates and domain labels (accessible via \code{refcolData})
#'     \item Fitted NNGP parameters for domain-informative genes (accessible via \code{paramsGP})
#'   }
#' @param seq_depth_factor Numeric scalar > 0. Scaling factor for sequencing depth.
#' @param effect_size_factor Numeric scalar > 0. Scaling factor for within-domain effect size.
#' @param SEED Integer random seed for reproducibility.
#'
#' @return A \code{spaDesign} object with simulated counts stored in \code{simCounts}
#' and spot metadata stored in \code{simcolData}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim_obj <- simulation_EffectSize_refactored(
#'   spaDesign = my_spadesign,
#'   seq_depth_factor = 1.2,
#'   effect_size_factor = 1.5,
#'   SEED = 1
#' )
#' }

simulation_EffectSize_refactored <- function(spaDesign, seq_depth_factor, effect_size_factor, SEED){
  
  if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) {
    stop("seq_depth_factor must be a positive numeric value.")
  }
  if (!is.numeric(effect_size_factor) || effect_size_factor <= 0) {
    stop("effect_size_factor must be a positive numeric value.")
  }
  
  # extract reference data
  count_matrix <- refCounts(spaDesign)
  all_spot_names <- colnames(count_matrix) # <-- Extract spot names once
  loc_data <- refcolData(spaDesign)[, c('x', 'y', 'domain')]
  par_GP <- paramsGP(spaDesign)
  
  ## scale the coordinates ... (this part is unchanged) ...
  coords_norm <- igraph::norm_coords(as.matrix(loc_data[,c('x', 'y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_data$domain
  
  seqDepth_factor <- seq_depth_factor
  ES_factor <- effect_size_factor
  
  message("Starting simulation...\n")
  all_count <- lapply(seq_along(par_GP), function(d){
    domain <- names(par_GP)[d]
    GP.par <- par_GP[[d]]
    
    idx <- which(coords_norm$domain == domain)
    coords_norm_sub <- coords_norm[idx, ]
    coords_norm_out <- coords_norm[-idx, ]
    
    message(sprintf("Simulating genes for domain: %s\n", domain))
    
    gene_count <- lapply(seq_along(GP.par), function(g) {
      tryCatch({
        gene <- names(GP.par)[g]
        # message(sprintf("Simulating gene: %s in domain: %s\n", gene, domain)) # Can be noisy
        nnGP_fit <- GP.par[[g]]
        
        # --- KEY CHANGE IS HERE ---
        # Check if gene exists and extract its single row of counts BEFORE calling the function
        if (!gene %in% rownames(count_matrix)) {
          warning(paste('Gene', gene, 'not found in count matrix for domain', domain))
          return(NULL)
        }
        gene_count_row <- count_matrix[gene, ]
        
        # Call the refactored function with only the necessary data
        tt <- simulate_genecount_ES_refactored(
          SEED, coords_norm_sub, coords_norm_out, nnGP_fit, 
          seqDepth_factor, ES_factor, 
          gene_count_row = gene_count_row,      
          all_spot_names = all_spot_names,      
          spot_idx = idx
        )
        return(tt)
      }, error = function(e) {
        warning(paste('Error in gene', gene, 'of domain', domain, ':', e$message))
        return(NULL)
      })
    })
    
    gene_count <- Filter(Negate(is.null), gene_count)
    
    if (length(gene_count) > 0) {
      if (length(gene_count) == 1) {
        gene_count <- t(data.frame(gene_count[[1]]))  # Ensure it's a data frame if only one row
      } else {
        gene_count <- do.call(rbind, gene_count)
      }
      colnames(gene_count) <- all_spot_names
      rownames(gene_count) <- paste0(domain, "-", names(GP.par))
    } else {
      warning(paste("All gene simulations for domain", domain, "returned NULL."))
      return(NULL)
    }
    
    return(gene_count)
  })
  
  
  all_count <- Filter(Negate(is.null), all_count)
  if (length(all_count) > 0) {
    all_count <- do.call(rbind, all_count)
    spaDesign@simCounts <- all_count
    spaDesign@simcolData <- refcolData(spaDesign)
    message("Simulation complete.\n")
    return(spaDesign)
  } else {
    stop("All domain simulations returned NULL.")
  }
}


#' Simulate Expression Counts for a Single Gene
#'
#' Simulates Poisson gene expression counts for one gene across all spatial spots.
#' Counts are generated for spots inside a domain using NNGP-predicted mean expression,
#' and for spots outside the domain using an empirical baseline mean estimated from the
#' observed data.
#'
#' The function supports global scaling of sequencing depth (\code{seqDepth_factor})
#' and within-domain effect size (\code{ES_factor}).
#'
#' @param SEED Integer random seed for reproducibility.
#' @param coords_norm_sub Data frame of normalized coordinates for spots within the domain.
#'   Must contain columns \code{x} and \code{y}.
#' @param coords_norm_out Data frame of normalized coordinates for spots outside the domain.
#'   Must contain columns \code{x} and \code{y}.
#' @param nnGP_fit Fitted NNGP model object for the gene.
#' @param seqDepth_factor Numeric scalar > 0. Scaling factor for sequencing depth.
#' @param ES_factor Numeric scalar > 0. Scaling factor for within-domain effect size.
#' @param counts Reference count matrix (genes x spots).
#' @param gene Character string specifying the gene name.
#' @param spot_idx Integer vector indicating indices of spots within the domain
#'   (relative to \code{colnames(counts)}).
#'
#' @return Numeric vector of simulated counts for all spots, ordered to match
#' \code{colnames(counts)}.
#'
#' @noRd

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
    
    #set.seed(SEED)
    y.out <- rpois(n = n2, lambda = seqDepth_factor * lambda.out)
    sim.out <- data.frame(x = coords_norm_out$x, y = coords_norm_out$y, sim = y.out)
    rownames(sim.out) <- rownames(coords_norm_out)
    
    count.sim <- rbind(sim.in, sim.out)
    count.sim <- count.sim[order(match(rownames(count.sim),colnames(counts))),]
    identical(rownames(count.sim),colnames(counts))
    tt <- count.sim$sim
    return(tt) 
}


#' Simulate Expression Counts for a Single Gene (Refactored)
#'
#' A refactored version of \code{\link{simulate_genecount_ES}} that reduces memory usage
#' by using a single gene count vector instead of the full count matrix.
#'
#' Outside-domain baseline expression is computed from \code{gene_count_row}, while the
#' ordering of spots is enforced using \code{all_spot_names}.
#'
#' @param SEED Integer random seed for reproducibility.
#' @param coords_norm_sub Data frame of normalized coordinates for spots within the domain.
#' @param coords_norm_out Data frame of normalized coordinates for spots outside the domain.
#' @param nnGP_fit Fitted NNGP model object for the gene.
#' @param seqDepth_factor Numeric scalar > 0. Scaling factor for sequencing depth.
#' @param ES_factor Numeric scalar > 0. Scaling factor for within-domain effect size.
#' @param gene_count_row Numeric vector of observed counts for the target gene across all spots.
#' @param all_spot_names Character vector of all spot names (canonical ordering).
#' @param spot_idx Integer vector indicating indices of spots within the domain
#'   (relative to \code{all_spot_names}).
#'
#' @return Numeric vector of simulated counts for all spots, ordered to match
#' \code{all_spot_names}.
#'
#' @noRd

simulate_genecount_ES_refactored <- function(SEED, coords_norm_sub, coords_norm_out, nnGP_fit, seqDepth_factor, ES_factor, gene_count_row, all_spot_names, spot_idx){
  
  b.condition <- mean_in(SEED, coords_norm_sub, nnGP_fit)
  
  # the inside mean
  mu1 <- sum(b.condition)
  n1 <- nrow(coords_norm_sub)
  
  # the outside mean
  mu.out <- mean_out_refactored(gene_count_row, spot_idx)
  n2 <- length(all_spot_names) - n1
  
  # calculate lambda_in and lambda_out
  A <- mu1 + n2 * mu.out
  lambda.in <- A * ES_factor * b.condition / (ES_factor * mu1 + n2 * mu.out)
  lambda.out <- mu.out * A / (ES_factor * mu1 + n2 * mu.out)
  
  set.seed(SEED)    
  y.post <- rpois(n = nrow(coords_norm_sub), lambda = seqDepth_factor * lambda.in)
  
  sim.in <- data.frame(x = coords_norm_sub$x, y = coords_norm_sub$y, sim = y.post)
  rownames(sim.in) <- rownames(coords_norm_sub)
  
  y.out <- rpois(n = n2, lambda = seqDepth_factor * lambda.out)
  sim.out <- data.frame(x = coords_norm_out$x, y = coords_norm_out$y , sim = y.out)
  rownames(sim.out) <- rownames(coords_norm_out)
  
  count.sim <- rbind(sim.in, sim.out)
  # use the passed-in spot names for ordering
  count.sim <- count.sim[order(match(rownames(count.sim), all_spot_names)),] 
  
  tt <- count.sim$sim
  return(tt)  
  
}
  


#' Predict Within-Domain Mean Expression from a Fitted NNGP Model
#'
#' Uses a fitted NNGP model to predict expected expression levels for each spot inside a domain.
#' The predicted values are transformed to the count scale and returned as a numeric vector,
#' which is used as the within-domain mean intensity in simulation.
#'
#' @param SEED Integer random seed for reproducibility.
#' @param coords_norm_sub Data frame of normalized coordinates for spots within the domain.
#'   Must contain columns \code{x} and \code{y}.
#' @param nnGP_fit Fitted NNGP model object for the gene (compatible with
#'   \code{BRISC::BRISC_prediction}).
#'
#' @return Numeric vector of predicted within-domain mean expression values.
#'
#' @noRd

mean_in <- function(SEED, coords_norm_sub, nnGP_fit) {
    if (missing(SEED) || !is.numeric(SEED)) stop("SEED must be a numeric value.")
    if (!is.data.frame(coords_norm_sub) || !all(c("x", "y") %in% names(coords_norm_sub))) stop("coords_norm_sub must be a data frame with 'x' and 'y' columns.")
	
	coords_spatial <- as.matrix(coords_norm_sub[, c('x', 'y')])
	temp_file <- tempfile()
	sink(temp_file)
	brisc_pred <- BRISC::BRISC_prediction(nnGP_fit, coords_spatial)
	sink()
	unlink(temp_file)
    sigma2 <- nnGP_fit$Theta['sigma.sq']
	tau2 <- nnGP_fit$Theta['tau.sq']   
    b.condition <- round(exp(brisc_pred$prediction + tau2/2))-1
    return(b.condition)
}

#' Compute Outside-Domain Baseline Mean Expression
#'
#' Calculates the outside-domain baseline mean for a given gene by taking the mean of the
#' lower half (<= median) of counts outside the domain. This provides a robust baseline
#' against high-expression outliers.
#'
#' @param counts Matrix of gene expression counts (genes x spots).
#' @param gene Character string specifying the gene name.
#' @param spot_idx Integer vector indicating indices of spots within the domain.
#'
#' @return Numeric scalar representing the outside-domain baseline mean expression.
#' @noRd
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


#' Compute Outside-Domain Baseline Mean Expression (Refactored)
#'
#' Refactored version of \code{mean_out} that works on a single gene count vector rather than a
#' full gene-by-spot matrix.
#'
#' @param gene_count_row Numeric vector of observed counts for a single gene across all spots.
#' @param spot_idx Integer vector indicating indices of spots within the domain.
#'
#' @return Numeric scalar representing the outside-domain baseline mean expression.
#' @noRd
mean_out_refactored <- function(gene_count_row, spot_idx){
  if(!is.numeric(spot_idx)) stop('spot_idx must be numeric.')
  
  OUT <- gene_count_row[-spot_idx]
  OUT.lower <- OUT[OUT <= median(OUT)]
  mean.out.low <- mean(OUT.lower)
  return(mean.out.low)
}



