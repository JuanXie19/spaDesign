#' simulate data with modified sequencing depth and dispersed pattern
#'
#' @param spadesign The \code{spadesign} object
#' @param seq_depth_factor A scaling factor for sequencing depth
#' @param SIGMA scaling factor for the \code{sigma} parameter from the FG model
#' @param SEED random seed for reproducibility
#' @param prop Proportion of genes with undisturbed expression pattern
#' @return A \code{spadesign} object with simulated count matrix
#' @import pdist
#' @import clue
#' @import RANN
#' @import invgamma
#' @import Rfast
#' @import movMF
#' @import MASS
#' @export
#'
simulateDataSpatial <- function(spadesign, seq_depth_factor, SIGMA, SEED, prop){
 
    if (!is.numeric(seq_depth_factor) || seq_depth_factor <= 0) {
        stop("seq_depth_factor must be a positive numeric value.")
    }
    
    count_matrix <- refCounts(spadesign)    
    loc_file <- refcolData(spadesign)[, c('x','y','domain')]
    par_GP <- paramsPoissonGP(spadesign)
    par_FG <- paramsFG(spadesign)
    
    # Check if par_GP and par_FG exist and are not NULL
    if (is.null(par_GP) || is.null(par_FG)) {
        stop("par_GP and/or par_FG do not exist or are NULL.Please run parameter estimation first")
    }
    
    #message("Starting simulation...")
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[, c('x','y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain
    
    seqDepth_factor <- seq_depth_factor
    
    
    message('Simualting baseline count matrix for all domain-informative genes...')
    ### simulate baseline count matrix for all domain-informative genes
    base_count <- pbapply::pblapply(seq_along(par_GP), function(d){
        
        domain <- names(par_GP)[d]
        GP.par <- par_GP[[d]]

        simulate_base_count(SEED, seqDepth_factor, domain, GP.par, count_matrix, coords_norm)        
        })%>% Reduce('rbind', .)
    
    ## simulate count matrix where the spots location are disturbed
    message('Simulating count matrix for disturbed spots location...')
    worse_count <- pbapply::pblapply(seq_along(par_GP), function(d){
        domain <- names(par_GP)[d]
        GP.par <- par_GP[[d]]
        FG.par <- par_FG[[d]]

        simulate_worse_count(SEED, seqDepth_factor, domain, GP.par, FG.par, count_matrix, coords_norm, SIGMA)
        }) %>% Reduce('rbind', .)

    generateCount <- function(){
        COUNT <- lapply(seq_along(par_GP), function(d){
            domain <- names(par_GP)[d]
            count.base <- base_count[grep(domain,rownames(base_count)), ]   
            count.worst <- worse_count[grep(domain,rownames(worse_count)), ]
            n <- nrow(count.base)
         
            set.seed(SEED)
            keep.idx <- sample(seq_len(n),round(prop * n))
            count_combined <- rbind(count.base[keep.idx, ], count.worst[-keep.idx, ])
            rownames(count_combined) <- c(rownames(count.base)[keep.idx],rownames(count.worst)[-keep.idx])
            count_combined <- count_combined[order(match(rownames(count_combined),rownames(count.base))),]
            return(count_combined)
        })%>% Reduce('rbind', .)
	return(COUNT)
    }
	
	REPEAT <- 1000
	COUNT.SIM <- lapply(seq_len(REPEAT), function(x) generateCount())%>% Reduce('+',.)/REPEAT
	message("Simulation complete.\n")
	spadesign@simCounts <- COUNT.SIM
	spadesign@simcolData <- refcolData(spadesign)
	return(spadesign)
}


## function to simulate gene expression for within domain spots
simulate_geneCounts_in <- function(SEED, seqDepth_factor, coords_norm_sub, alpha.est, rho.est, mu.est, s.est){
  
  set.seed(SEED)
  b.condition <- mean_in(SEED, coords_norm_sub, alpha.est, rho.est, mu.est, s.est)
  
  ## inside domain, seq depth R
  set.seed(SEED)
  y.post <- rpois(n = nrow(coords_norm_sub), lambda = seqDepth_factor * exp(b.condition))
  return(y.post)
}

## function to simulate gene expression for outside domain spots
simulate_geneCounts_out <- function(SEED, seqDepth_factor, counts, gene, spot_idx, COORDS.OUT){
   
  mean.out.low <- mean_out(counts, gene, spot_idx)
  set.seed(SEED)
  
  ## modify the outside gene expression by seqDepth factor R
  y.out <- rpois(n = nrow(COORDS.OUT),lambda = seqDepth_factor * mean.out.low)
  return(y.out)
}

## function to combine the inside and outside domain gene expression
combine_in_out <- function(COORDS.IN, COORDS.OUT, y.post, y.out, counts){
  sim.in <- data.frame(x = COORDS.IN$x,y = COORDS.IN$y, sim = y.post)
  rownames(sim.in) <- rownames(COORDS.IN)
  
  sim.out <- data.frame(x = COORDS.OUT$x, y = COORDS.OUT$y,sim = y.out)
  rownames(sim.out) <- rownames(COORDS.OUT)
  
  count.sim <- rbind(sim.in, sim.out)
  count.sim <- count.sim[order(match(rownames(count.sim), colnames(counts))), ]
  return(count.sim)
}

## process the gene expression: if GP conditional sampling has an error, null value will be assigned to the gene, here to remove null expression generated 
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
 
# simulate base count for one set of domain-informative genes. Specifically, simulate count for multiple genes, including their within-domain and outside domain expression

simulate_base_count <- function(SEED, seqDepth_factor, domain, GP.par, counts, coords_norm){
    
    spot_idx <- which(coords_norm$domain == domain)        
    coords_norm_sub <- coords_norm[spot_idx, ]
    coords_norm_out <- coords_norm[-spot_idx, ]

    GENE.COUNT <- vector('list',length(GP.par))
  
    for (g in seq_along(GP.par)){
        tryCatch({
            gene <- names(GP.par)[g]
            alpha.est <- GP.par[[g]][['alpha.est']]
            rho.est <- GP.par[[g]][['rho.est']]
            mu.est <- GP.par[[g]][['mu.est']]
            s.est <- GP.par[[g]][['s.est']]
            
            ## simulate the inside domain gene expression
            y.post <- simulate_geneCounts_in(SEED, seqDepth_factor, coords_norm_sub, alpha.est, rho.est, mu.est, s.est)
            
            ## simulate the outside gene expression: based on the lower 50% of the values
            y.out <- simulate_geneCounts_out(SEED, seqDepth_factor, counts, gene, spot_idx, coords_norm_out)
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

simulate_worse_count <- function(SEED, seqDepth_factor, domain, GP.par, FG.par, counts, coords_norm, SIGMA){

  ## modulate the spatial pattern via changing the sigma in the FG fit
    FG.sim <- FG.par
    FG.sim$trace_sigma <- SIGMA * FG.par$trace_sigma
    FG.sim$tau <- FG.par$tau * 2
  
    spot_idx <- which(coords_norm$domain == domain)
    coords_norm_sub <- coords_norm[spot_idx,]
    coords_norm_sub <- as.data.frame(coords_norm_sub)
    #colnames(coords_norm_sub) <- c('x','y')
  
 # simulate new locations
  
    n.loc <- nrow(coords_norm_sub)
    N.loc <- round(n.loc * 2)
  
  
    set.seed(123456)
    pred.loc1 <- plot_density_FGM(N.loc, FG.sim)[[2]]
    pred.loc2 <- plot_density_FGM(N.loc, FG.sim)[[2]]
    pred.loc3 <- plot_density_FGM(N.loc, FG.sim)[[2]]
    pred.loc <- rbind(pred.loc1,pred.loc2,pred.loc3)
    colnames(pred.loc) <- c('x', 'y')
  
  ## assign simulated locs to its NN points in the original coordinates
  
    source_matrix <- pred.loc
    target_matrix <- coords_norm[, c('x', 'y')]
  
    nearest_neighbors <- Nearest_RANN(source_matrix, target_matrix)
  
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
      alpha.est <- GP.par[[g]][['alpha.est']]
      rho.est <- GP.par[[g]][['rho.est']]
      mu.est <- GP.par[[g]][['mu.est']]
      s.est <- GP.par[[g]][['s.est']]
      
      # inside domain gene expression
      y.post <- simulate_geneCounts_in(SEED, seqDepth_factor, coords_norm_sub, alpha.est, rho.est, mu.est, s.est)
      
      ## simulate outside domain gene expression
      y.out <- simulate_geneCounts_out(SEED, seqDepth_factor, counts, gene, spot_idx, REST)
      
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

    
#' get the nearest neighbor index
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
