#' Simulate data with modified sequencing depth and dispersion
#'
#' This function generates simulated spatial gene expression counts by combining baseline domain-informative
#' expression patterns with perturbed spot locations. i.e.,creates a mixture where a proportion of genes retain their original spatial
#' patterns while others are simulated with disturbed locations.
#' It allows scaling of sequencing depth and the \code{sigma} parameter from the FG model, and can keep a proportion of genes undisturbed.
#' 
#' @param spaDesign A \code{spaDesign} object containing:
#'   \itemize{
#'     \item Spatial coordinates and gene expression data
#'     \item Fitted GP parameters (from \code{paramsGP} slot)
#'     \item Fitted FG parameters (from \code{paramsFG} slot)
#'     \item Domain assignments
#'   }
#' @param selected_M_list Optional integer vector specifying which FG model (number of clusters M)
#'   to use for each domain. If \code{NULL}, uses \code{spaDesign@selected_M_list_BIC}.
#'   Length must match number of domains.
#' @param seq_depth_factor Numeric scalar > 0. Multiplicative scaling factor for sequencing depth.
#'   Values > 1 increase depth; < 1 decrease depth.
#' @param SIGMA Numeric scalar > 0. Multiplicative scaling factor for the \code{sigma_sq} parameter
#'   from the FG model. Larger values increase spatial noise/dispersion.
#' @param SEED Integer. Random seed for reproducibility of GP and FG sampling.
#' @param prop Numeric in [0,1]. Proportion of genes to keep with original (undisturbed) 
#'   spatial patterns. The remaining (1-prop) genes will have perturbed locations.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing via \code{mclapply}.
#'
#'
#' @details
#' The simulation proceeds in three steps:
#' \enumerate{
#'   \item Generate baseline counts for all genes using fitted GP models
#'   \item Generate perturbed counts using FG-model-based location disturbance
#'   \item Create final counts by randomly mixing baseline and perturbed patterns,
#'         with averaging over 1000 random selections
#' }
#' 
#' The FG model parameters are modified as follows:
#' \itemize{
#'   \item \code{sigma_sq} is multiplied by \code{SIGMA}
#'   \item \code{tau} (concentration) is multiplied by 2 to increase spatial perturbation
#' }
#'
#' @return A \code{spaDesign} object with updated slots:
#'   \describe{
#'     \item{simCounts}{Numeric matrix (genes × spots) of simulated expression counts,
#'                      averaged over 1000 replications of random gene selection}
#'     \item{simcolData}{Data frame of spot metadata, copied from \code{refcolData(spaDesign)}}
#'   }
#'
#' @import pdist
#' @import clue
#' @import RANN
#' @import invgamma
#' @import Rfast
#' @import movMF
#' @import MASS
#' @import dplyr
#' @import pbapply
#' @import pbmcapply
#' @import future.apply
#' @importFrom igraph norm_coords
#' @importFrom parallel mclapply
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming spadesign_obj has fitted GP and FG parameters
#' simulated_data <- simulation_Spatial(
#'   spaDesign = spadesign_obj,
#'   selected_M_list = c(3, 4, 5),  # Use M=3,4,5 for domains 1,2,3
#'   seq_depth_factor = 1.5,         # 50% increase in sequencing depth
#'   SIGMA = 1.2,                    # 20% increase in spatial noise
#'   SEED = 123,
#'   prop = 0.7,                     # Keep 70% genes undisturbed
#'   n_cores = 4
#' )
#' 
#' # Access simulated counts
#' sim_counts <- simulated_data@simCounts
#' }



simulation_Spatial <- function(spaDesign, selected_M_list = NULL, seq_depth_factor, SIGMA, SEED, prop, n_cores){
  
  if (is.null(selected_M_list)) {
    if (!is.null(spaDesign@selected_M_list_BIC)) {
      selected_M_list <- spaDesign@selected_M_list_BIC
    } else {
      stop("No selected_M_list provided and spaDesign does not contain selected_M_list_BIC.")
    }
  }
  
  # Validate inputs
  if (length(selected_M_list) != length(spaDesign@paramsFG)) {
    stop("Length of selected_M_list must match number of domains in spaDesign")
  }
  
  if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) {
    stop("seq_depth_factor must be a positive numeric value.")
  }
  
  count_matrix <- refCounts(spaDesign)    
  loc_file <- refcolData(spaDesign)[, c('x','y','domain')]
  par_GP <- paramsGP(spaDesign)
  par_FG <- paramsFG(spaDesign)
  
  # Check if par_GP and par_FG exist and are not NULL
  if (is.null(par_GP) || is.null(par_FG)) {
    stop("par_GP and/or par_FG do not exist or are NULL.Please run parameter estimation first")
  }
  
  ## scale the coordinates to [0,1] range
  coords_norm <- igraph::norm_coords(as.matrix(loc_file[, c('x','y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_file$domain
  
  ## extract FG models
  FG_selected_model <- list()
  for (d in 1:length(spaDesign@paramsFG)){
    domain <- names(spaDesign@paramsFG)[d]
    M <- selected_M_list[d]
    model_name <- paste0("M_", M)
    domain_models <- spaDesign@paramsFG[[domain]]$all_models
    if (!model_name %in% names(domain_models)) {
      available_M <- gsub("M_", "", names(domain_models))
      stop(paste0("M =", M, " not found for domain ", domain, 
                  ". Available M values: ", paste(available_M, collapse = ", ")))
    }
    
    FG_selected_model[[d]] <- domain_models[[model_name]]
  }
  names(FG_selected_model) <- names(spaDesign@paramsFG)
  
  
  seqDepth_factor <- seq_depth_factor
  
  
  message('Simualting baseline count matrix for all domain-informative genes...')
  ### simulate baseline count matrix for all domain-informative genes
  base_count <- pbapply::pblapply(seq_along(par_GP), function(d){
    
    domain <- names(par_GP)[d]
    GP.par <- par_GP[[d]]
    
    simulate_base_count(SEED, seqDepth_factor, domain, GP.par, count_matrix, coords_norm)        
  })
  base_count <- do.call(rbind, base_count)
  
  ## simulate count matrix where the spots location are disturbed
  message('Simulating count matrix for disturbed spots location...')
  
  worse_count <- mclapply(seq_along(par_GP), function(d){
    domain <- names(par_GP)[d]
    GP.par <- par_GP[[d]]
    FG.par <- FG_selected_model[[d]]
    
    simulate_worse_count(SEED, seqDepth_factor, domain, GP.par, FG.par, count_matrix, coords_norm, SIGMA)
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  # Check for errors
  if (any(sapply(worse_count, inherits, "try-error"))) {
    stop("Error in parallel processing of worse_count. Try reducing n_cores or check error messages.")
  }
  
  worse_count <- do.call('rbind', worse_count)
  
  generateCount_fast <- function(){
    count_combined <- base_count
    disturbed_genes <- sample(rownames(base_count), round((1 - prop) * nrow(base_count)))
    count_combined[disturbed_genes, ] <- worse_count[disturbed_genes, ]
    count_combined
  }
  
  REPEAT <- 1000
  COUNT.SIM_LIST <- future.apply::future_lapply(seq_len(REPEAT), function(i) {
    generateCount_fast()
  }, future.seed = TRUE)
  COUNT.SIM <- Reduce(`+`, COUNT.SIM_LIST) / REPEAT
  
  message("Simulation complete.\n")
  spaDesign@simCounts <- COUNT.SIM
  spaDesign@simcolData <- refcolData(spaDesign)
  return(spaDesign)
}


#' Simulate gene expression for spots inside a domain
#'
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param coords_norm_sub Data frame of normalized coordinates for spots inside the domain.
#' @param nnGP_fit Fitted NNGP model object for the gene.
#' @return Numeric vector of simulated counts for the spots inside the domain.
#' @noRd
simulate_geneCounts_in <- function(SEED, seqDepth_factor, coords_norm_sub, nnGP_fit){
  
  set.seed(SEED)
  b.condition <- mean_in(SEED, coords_norm_sub, nnGP_fit)
  
  ## inside domain, seq depth R
  y.post <- rpois(n = nrow(coords_norm_sub), lambda = seqDepth_factor * b.condition)
  return(y.post)
}

#' Simulate gene expression for spots outside a domain
#'
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param counts Matrix of gene expression counts (genes x spots).
#' @param gene Character string specifying the gene name.
#' @param spot_idx Numeric vector of indices of spots inside the domain.
#' @param COORDS.OUT Data frame of coordinates for spots outside the domain.
#' @return Numeric vector of simulated counts for the spots outside the domain.
#' @noRd
simulate_geneCounts_out <- function(SEED, seqDepth_factor, counts, gene, spot_idx, COORDS.OUT){
   
  mean.out.low <- mean_out(counts, gene, spot_idx)
  set.seed(SEED)
  
  ## modify the outside gene expression by seqDepth factor R
  y.out <- rpois(n = nrow(COORDS.OUT),lambda = seqDepth_factor * mean.out.low)
  return(y.out)
}


#' (Refactored) Simulate gene expression for spots outside a domain
#'
#' Uses a single gene vector instead of full count matrix to reduce memory usage.
#'
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param gene_count_row Numeric vector of counts for the gene.
#' @param spot_idx Numeric vector of indices of spots inside the domain.
#' @param COORDS.OUT Data frame of coordinates for spots outside the domain.
#' @return Numeric vector of simulated counts for the spots outside the domain.
#' @noRd
simulate_geneCounts_out_refactored <- function(SEED, seqDepth_factor, gene_count_row, spot_idx, COORDS.OUT){
  
  mean.out.low <- mean_out_refactored(gene_count_row, spot_idx)
  set.seed(SEED)
  
  ## modify the outside gene expression by seqDepth factor R
  y.out <- rpois(n = nrow(COORDS.OUT),lambda = seqDepth_factor * mean.out.low)
  return(y.out)
}



#' Combine simulated inside and outside domain counts
#'
#' @param COORDS.IN Data frame of coordinates for spots inside the domain.
#' @param COORDS.OUT Data frame of coordinates for spots outside the domain.
#' @param y.post Numeric vector of simulated counts for inside-domain spots.
#' @param y.out Numeric vector of simulated counts for outside-domain spots.
#' @param counts Original count matrix (for column ordering).
#' @return Data frame of simulated counts for all spots, ordered to match the original count matrix.
#' @noRd
combine_in_out <- function(COORDS.IN, COORDS.OUT, y.post, y.out, counts){
  sim.in <- data.frame(x = COORDS.IN$x,y = COORDS.IN$y, sim = y.post)
  rownames(sim.in) <- rownames(COORDS.IN)
  
  sim.out <- data.frame(x = COORDS.OUT$x, y = COORDS.OUT$y,sim = y.out)
  rownames(sim.out) <- rownames(COORDS.OUT)
  
  count.sim <- rbind(sim.in, sim.out)
  count.sim <- count.sim[order(match(rownames(count.sim), colnames(counts))), ]
  return(count.sim)
}

#' Process simulated counts and remove NULL entries
#'
#' @param GENE.COUNT List of simulated gene counts (each element a numeric vector or NULL).
#' @param counts Original count matrix (genes x spots).
#' @param domain Character string of domain name.
#' @param GP.par List of GP parameters for the genes in the domain.
#' @return Matrix of simulated counts for all genes in the domain, with appropriate row and column names.
#' @noRd
process_count <- function(GENE.COUNT, counts, domain, GP.par){
  if (any(sapply(GENE.COUNT, is.null)) == 'FALSE') {
    GENE.COUNT <- Reduce('rbind', GENE.COUNT)     
    rownames(GENE.COUNT) <- paste0(domain, "-", names(GP.par))
  } else {
    null.idx <- which(sapply(GENE.COUNT, is.null) == 'TRUE')
    GENE.COUNT <- Reduce('rbind', GENE.COUNT)
    rownames(GENE.COUNT) <- paste0(domain, "-", names(GP.par)[-null.idx])
  }
  colnames(GENE.COUNT) <- colnames(counts)
  return(GENE.COUNT)
}
 
#' Simulate expression counts for a set of domain-informative genes without modifying their spatial patterns.
#'
#' Generates baseline (undisturbed) expression counts for genes in a domain using fitted GP models.
#' Inside-domain expression follows the GP model; outside-domain expression uses the lower 50% quantile.
#'
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param domain Character string specifying the domain.
#' @param GP.par List of GP parameters for the domain genes.
#' @param counts Original count matrix (genes x spots).
#' @param coords_norm Data frame of normalized coordinates for all spots.
#' @return Matrix of simulated expression counts for the domain genes.
#' @noRd
simulate_base_count <- function(SEED, seqDepth_factor, domain, GP.par, counts, coords_norm){
    
    spot_idx <- which(coords_norm$domain == domain)        
    coords_norm_sub <- coords_norm[spot_idx, ]
    coords_norm_out <- coords_norm[-spot_idx, ]

    GENE.COUNT <- vector('list',length(GP.par))
  
    for (g in seq_along(GP.par)){
        tryCatch({
            gene <- names(GP.par)[g]
            nnGP_fit <- GP.par[[g]]
            
            gene_count_row <- counts[gene, ]
            ## simulate the inside domain gene expression
            y.post <- simulate_geneCounts_in(SEED, seqDepth_factor, coords_norm_sub, nnGP_fit)
            
            ## simulate the outside gene expression: based on the lower 50% of the values
            y.out <- simulate_geneCounts_out_refactored(SEED, seqDepth_factor, gene_count_row, spot_idx, coords_norm_out)
            count.sim <- combine_in_out(coords_norm_sub, coords_norm_out, y.post, y.out, counts)
            identical(rownames(count.sim),colnames(counts))
            tt <- count.sim$sim
            GENE.COUNT[[g]] <- tt
            },error = function(e) {
            GENE.COUNT[[g]] <- NULL
            message(paste('Error in simulate_base_count for gene',gene,'in domain',domain,":",e$message))
        })
    }
  
    GENE.COUNT <- process_count(GENE.COUNT,counts,domain,GP.par)
    return(GENE.COUNT)
}

#' Simulate perturbed expression counts for a set of domain-informative genes.
#'
#' Generates expression counts with spatially perturbed spot locations using the FG model.
#' New locations are sampled from a modified FG distribution, then matched to original spots
#' via Hungarian algorithm.
#' 
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param domain Character string specifying the domain.
#' @param GP.par List of GP parameters for the domain genes.
#' @param counts Original count matrix (genes x spots).
#' @param coords_norm Data frame of normalized coordinates for all spots.
#' @param FG.par FG model parameters for the domain.
#' @param SIGMA Numeric scaling factor for FG sigma.
#' @return Matrix of simulated counts for the domain genes.
#' @noRd
simulate_worse_count <- function(SEED, seqDepth_factor, domain, GP.par, FG.par, counts, coords_norm, SIGMA){

  ## modulate the spatial pattern via changing the sigma in the FG fit
    FG.sim <- FG.par$model
    FG.sim$sigma_sq <- SIGMA * FG.par$model$sigma_sq
    FG.sim$tau <- FG.par$model$tau * 2
  
    spot_idx <- which(coords_norm$domain == domain)
    coords_norm_sub <- coords_norm[spot_idx,]
    coords_norm_sub <- as.data.frame(coords_norm_sub)
    #colnames(coords_norm_sub) <- c('x','y')
  
 # simulate new locations
  
    n.loc <- nrow(coords_norm_sub)
    N.loc <- round(n.loc * 2)
  
  
    set.seed(SEED)
    pred.loc1 <- plot_density_FG_EM(N.loc, FG.sim)[[2]]
    set.seed(SEED + 1)
    pred.loc2 <- plot_density_FG_EM(N.loc, FG.sim)[[2]]
    set.seed(SEED + 2)
    pred.loc3 <- plot_density_FG_EM(N.loc, FG.sim)[[2]]
    pred.loc <- rbind(pred.loc1,pred.loc2,pred.loc3)
    colnames(pred.loc) <- c('x', 'y')
  
  ## assign simulated locs to its NN points in the original coordinates
  
    source_matrix <- pred.loc
    target_matrix <- coords_norm[, c('x', 'y')]
  
    nearest_neighbors <- Nearest_RANN_vectorized(source_matrix, target_matrix)
  
    coords_sim_nn <- as.data.frame(cbind(pred.loc, nearest_neighbors))
    colnames(coords_sim_nn) <- c('x', 'y', 'nearest_neighbors')
    coords_sim_nn <- coords_sim_nn %>% filter(nearest_neighbors >0)
  
    NN <- coords_norm[coords_sim_nn$nearest_neighbors,]
   ## note that some points are assigned to the same point, just keep the uniq
    NN.names <- sapply(strsplit(rownames(NN),'[.]'),'[',1)
    NN.unique <- NN[!duplicated(NN.names),]  # this is the final coordinates for the simulated spots
    nrow(NN.unique)
  
    if(nrow(NN.unique) > n.loc){
      NN.unique <- NN.unique[sample(1:nrow(NN.unique),n.loc),]
    }
  
    NN.INDX <- which(rownames(coords_norm) %in% rownames(NN.unique))
    REST <- as.data.frame(coords_norm[-NN.INDX,])
  
  #### find the best matching between simulated points and original domain points using Hungarian algorithm
  
    distances <- as.matrix(pdist::pdist(coords_norm_sub[, c('x', 'y')], NN.unique[,c('x', 'y')]))  # distance between existing loc and simulated loc
    sol <- clue::solve_LSAP(distances)
  
    reorder.syn <- as.vector(sol)
  
    sim_newloc <- as.data.frame(NN.unique[reorder.syn, ])
  
    GENE.COUNT <- vector('list',length(GP.par))
  
  for (g in seq_along(GP.par)){
    tryCatch({
      gene <- names(GP.par)[g]
      nnGP_fit <- GP.par[[g]]
      
      gene_count_row <- counts[gene, ]
      # inside domain gene expression
      y.post <- simulate_geneCounts_in(SEED, seqDepth_factor, coords_norm_sub, nnGP_fit)
      
      ## simulate outside domain gene expression
      y.out <- simulate_geneCounts_out_refactored(SEED, seqDepth_factor, gene_count_row, spot_idx, REST)
      
      count.sim <- combine_in_out(sim_newloc, REST, y.post, y.out, counts)                
      
      identical(rownames(count.sim),colnames(counts))
      tt <- count.sim$sim
      GENE.COUNT[[g]] <- tt 
    },error = function(e){
        GENE.COUNT[[g]] <- NULL
        message(paste('Error in simulate_worse_count for gene',gene, 'in domain', domain,":",e$message))
    } )
  }
  
  GENE.COUNT <- process_count(GENE.COUNT,counts,domain,GP.par)
  return(GENE.COUNT)
}

    
#' Get nearest neighbor indices between two sets of points
#'
#' Finds the nearest neighbor in the target matrix for each point in the source matrix.
#' Points whose nearest neighbor is farther than the minimum pairwise distance in the
#' target set are marked as outliers with index -1.
#'
#' @param source_matrix Matrix of coordinates for source points (rows = points, columns = x,y).
#' @param target_matrix Matrix of coordinates for target points (rows = points, columns = x,y).
#' @return Integer vector of nearest neighbor indices for each source point in the target matrix.
#'   If no neighbor is within minimal distance, the index is -1.
#' @noRd
Nearest_RANN <- function(source_matrix, target_matrix) {
    nn_result <- RANN::nn2(target_matrix, source_matrix, k = 1)
    nearest_indices <- nn_result$nn.idx
    distances <- nn_result$nn.dists
  
    dist_min <- min(dist(target_matrix))
  
    temp <- sapply(1:nrow(source_matrix), function(row) {
        if (min(distances[row, ]) > dist_min) {
            temp <- -1
        } else {
            temp <- nearest_indices[row]
        }
        temp
    })
  
    return(temp)
}


#' Get nearest neighbor indices (vectorized version)
#'
#' Vectorized implementation of nearest neighbor search. Finds the nearest neighbor
#' in the target matrix for each point in the source matrix.
#'
#' @param source_matrix Numeric matrix of coordinates for source points (rows = points, columns = x,y).
#' @param target_matrix Numeric matrix of coordinates for target points (rows = points, columns = x,y).
#' @return Integer vector of nearest neighbor indices for each source point in the target matrix.
#'   Returns -1 for source points whose nearest neighbor distance exceeds the minimum
#'   pairwise distance within the target set (i.e., isolated outliers).
#' @noRd
Nearest_RANN_vectorized <- function(source_matrix, target_matrix) {
  nn_result <- RANN::nn2(target_matrix, source_matrix, k = 1)
  nearest_indices <- as.vector(nn_result$nn.idx)
  distances <- as.vector(nn_result$nn.dists)
  
  dist_min <- min(dist(target_matrix))
  
  result <- nearest_indices
  result[distances > dist_min] <- -1
  
  return(result)
}

