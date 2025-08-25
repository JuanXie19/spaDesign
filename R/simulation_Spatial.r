#' Simulate data with modified sequencing depth and dispersion
#'
#' This function generates simulated spatial gene expression counts by combining baseline domain-informative
#' expression patterns with perturbed spot locations. It allows scaling of sequencing depth and the 
#' \code{sigma} parameter from the FG model, and can keep a proportion of genes undisturbed.
#' 
#' @param shinyDesign A \code{shinyDesign} object containing spatial coordinates, gene expression,
#'   and fitted GP/FG parameters.
#' @param selected_M_list Optional list of selected FG model indices for each domain. If \code{NULL},
#'   the function will attempt to use \code{shinyDesign@selected_M_list_BIC}.
#' @param seq_depth_factor Numeric scaling factor for sequencing depth
#' @param SIGMA Numeric scaling factor for the \code{sigma} parameter from the FG model.
#' @param SEED Integer random seed for reproducibility.
#' @param prop Numeric, proportion of genes to keep undisturbed (between 0 and 1).
#' @return A \code{shinyDesign} object with simulated count matrix stored in \code{simCounts}, and
#'   updated spot metadata in \code{simcolData}.
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
#' @export
#'
#' @examples
#' \dontrun{
#' simulated_data <- simulation_Spatial(spadesign_obj, 
#'                                    seq_depth_factor = 1.5,
#'                                    SIGMA = 1.2,
#'                                    SEED = 123,
#'                                    prop = 0.7)
#' }

simulation_Spatial <- function(shinyDesign, selected_M_list = NULL, seq_depth_factor, SIGMA, SEED, prop){
  
  if (is.null(selected_M_list)) {
    if (!is.null(shinyDesign@selected_M_list_AIC)) {
      selected_M_list <- shinyDesign@selected_M_list_BIC
    } else {
      stop("No selected_M_list provided and shinyDesign does not contain selected_M_list_AIC.")
    }
  }
  
  # Validate inputs
  if (length(selected_M_list) != length(shinyDesign@paramsFG)) {
    stop("Length of selected_M_list must match number of domains in shinyDesign")
  }
  
  if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) {
    stop("seq_depth_factor must be a positive numeric value.")
  }
  
  count_matrix <- refCounts(shinyDesign)    
  loc_file <- refcolData(shinyDesign)[, c('x','y','domain')]
  par_GP <- paramsGP(shinyDesign)
  par_FG <- paramsFG(shinyDesign)
  
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
  for (d in 1:length(shinyDesign@paramsFG)){
    domain <- names(shinyDesign@paramsFG)[d]
    M <- selected_M_list[d]
    model_name <- paste0("M_", M)
    domain_models <- shinyDesign@paramsFG[[domain]]$all_models
    if (!model_name %in% names(domain_models)) {
      available_M <- gsub("M_", "", names(domain_models))
      stop(paste0("M =", M, " not found for domain ", domain, 
                  ". Available M values: ", paste(available_M, collapse = ", ")))
    }
    
    FG_selected_model[[d]] <- domain_models[[model_name]]
  }
  names(FG_selected_model) <- names(shinyDesign@paramsFG)
  
  
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
  

  worse_count <- pbmcapply::pbmclapply(seq_along(par_GP), function(d){
    domain <- names(par_GP)[d]
    GP.par <- par_GP[[d]]
    FG.par <- FG_selected_model[[d]]
    
    simulate_worse_count(SEED, seqDepth_factor, domain, GP.par, FG.par, count_matrix, coords_norm, SIGMA)
  }, mc.cores = 4)
  worse_count <- do.call('rbind', worse_count)
  
  generateCount <- function(){
    COUNT <- lapply(seq_along(par_GP), function(d){
      domain <- names(par_GP)[d]
      count.base <- base_count[grep(domain, rownames(base_count)), ]   
      count.worst <- worse_count[grep(domain, rownames(worse_count)), ]
      n <- nrow(count.base)
      
      set.seed(SEED)
      keep.idx <- sample(seq_len(n),round(prop * n))
      count_combined <- rbind(count.base[keep.idx, ], count.worst[-keep.idx, ])
      rownames(count_combined) <- c(rownames(count.base)[keep.idx],rownames(count.worst)[-keep.idx])
      count_combined <- count_combined[order(match(rownames(count_combined),rownames(count.base))),]
      return(count_combined)
    })
    COUNT <- do.call('rbind', COUNT)
    return(COUNT)
  }
  
  generateCount_fast <- function(){
    count_combined <- base_count
    disturbed_genes <- sample(rownames(base_count), round((1 - prop) * nrow(base_count)))
    count_combined[disturbed_genes, ] <- worse_count[disturbed_genes, ]
    count_combined
  }
  
  REPEAT <- 1000
  COUNT.SIM_LIST <- future_lapply(seq_len(REPEAT), function(i) {
    generateCount_fast()
  }, future.seed = TRUE)
  COUNT.SIM <- Reduce(`+`, COUNT.SIM_LIST) / REPEAT
  
  if(FALSE){
  REPEAT <- 1000
  #COUNT.SIM <- lapply(seq_len(REPEAT), function(x) generateCount())%>% Reduce('+',.)/REPEAT
  
  COUNT.SIM <- replicate(REPEAT, generateCount(), simplify = F)
  COUNT.SIM <- Reduce('+', COUNT.SIM)/REPEAT
  }
  message("Simulation complete.\n")
  shinyDesign@simCounts <- COUNT.SIM
  shinyDesign@simcolData <- refcolData(shinyDesign)
  return(shinyDesign)
}


#' Simulate gene expression for spots inside a domain
#'
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param coords_norm_sub Data frame of normalized coordinates for spots inside the domain.
#' @param nnGP_fit Fitted NNGP model object for the gene.
#' @return Numeric vector of simulated counts for the spots inside the domain.
#' @export
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
#' @export
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
#' @export
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
#' @export

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
#' @export
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
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param domain Character string specifying the domain.
#' @param GP.par List of GP parameters for the domain genes.
#' @param counts Original count matrix (genes x spots).
#' @param coords_norm Data frame of normalized coordinates for all spots.
#' @return Matrix of simulated expression counts for the domain genes.
#' @export
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
#' @param SEED Integer random seed for reproducibility.
#' @param seqDepth_factor Numeric scaling factor for sequencing depth.
#' @param domain Character string specifying the domain.
#' @param GP.par List of GP parameters for the domain genes.
#' @param counts Original count matrix (genes x spots).
#' @param coords_norm Data frame of normalized coordinates for all spots.
#' @param FG.par FG model parameters for the domain.
#' @param SIGMA Numeric scaling factor for FG sigma.
#' @return Matrix of simulated counts for the domain genes.
#' @export
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
  
  
    set.seed(123456)
    pred.loc1 <- plot_density_FG_EM(N.loc, FG.sim)[[2]]
    pred.loc2 <- plot_density_FG_EM(N.loc, FG.sim)[[2]]
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
  
    distances <- as.matrix(pdist(coords_norm_sub[, c('x', 'y')], NN.unique[,c('x', 'y')]))  # distance between existing loc and simulated loc
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
#' @param source_matrix Matrix of coordinates for source points (rows = points, columns = x,y).
#' @param target_matrix Matrix of coordinates for target points (rows = points, columns = x,y).
#' @return Integer vector of nearest neighbor indices for each source point in the target matrix.
#'   If no neighbor is within minimal distance, the index is -1.
#' @export
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


# Vectorized function
Nearest_RANN_vectorized <- function(source_matrix, target_matrix) {
  nn_result <- RANN::nn2(target_matrix, source_matrix, k = 1)
  nearest_indices <- as.vector(nn_result$nn.idx)
  distances <- as.vector(nn_result$nn.dists)
  
  dist_min <- min(dist(target_matrix))
  
  result <- nearest_indices
  result[distances > dist_min] <- -1
  
  return(result)
}

