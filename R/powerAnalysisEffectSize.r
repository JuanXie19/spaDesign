#' Power analysis 
#' 
#' This function performs power analysis by simulating data with different effect sizes and sequencing depths.
#' 
#' @param shinyDesign A \code{shinyDesign} object with estimated model parameters
#' @param es_range A numeric vector indicating the range of effect size scaling factors
#' @param seq_depth_range A numeric vector indicating the range of sequencing depth scaling factors
#' @param n_rep An integer indicating the number of replications for each scenario
#' @param conda_env_path Specify the path to the conda environment for spaGCN
#' @param python_script Specify the path to the python script
#' @return A data frame containing the sequencing depth, effect size, total counts, and NMI
#' @name powerAnalysisEffectSize
#' @import parallel
#' @import pbapply
#' @import foreach
#' @import doParallel
#' @examples
#' # spadesign <- ...
#' # es_range <- c(1, 2, 3)
#' # seq_depth_range <- c(0.5, 1, 3, 5)
#' # n_rep <- 10
#' # results <- powerAnalysisEffectSize(shinyDesign, es_range, seq_depth_range, n_rep, conda_env_path) 
#' @export


powerAnalysisEffectSize <- function(shinyDesign, es_range, seq_depth_range, n_rep) {
    # Create a data frame of all combinations of es_range, seq_depth_range, and n_rep
    param_grid <- expand.grid(seq_depth = seq_depth_range, 
                              effect_size = es_range, 
                              SEED = seq_len(n_rep))

    # Define a function to process each row of the param_grid
    process_row <- function(row) {
        effect_size_factor <- as.numeric(row["effect_size"])
        seq_depth_factor <- as.numeric(row["seq_depth"])
        SEED <- as.integer(row["SEED"])
        
        # Print values for debugging
        message(sprintf("Simulating data with seq_depth_factor = %s, effect_size_factor = %s, SEED = %s", 
                        seq_depth_factor, effect_size_factor, SEED))
    
        shinyDesign <- simulation_EffectSize_refactored(shinyDesign, seq_depth_factor, effect_size_factor, SEED)
        
        message(sprintf("Evaluating power for simulated data with seq_depth_factor = %s, effect_size_factor = %s, SEED = %s", 
                        seq_depth_factor, effect_size_factor, SEED))
        
        #rst <- evaluatePower(shinyDesign, conda_env_path)
		    rst <- evaluatePowerSeurat(shinyDesign)
        NMI <- rst@NMI
        total_counts <- sum(simCounts(shinyDesign))
        
        data.frame(seq_depth = seq_depth_factor, effect_size = effect_size_factor, 
                   total_counts = total_counts, NMI = NMI)
    }
    
    
    
    results <- mclapply(1:nrow(param_grid), function(i){
      process_row(param_grid[i, ])
    }, mc.cores = 4)
    results <- do.call(rbind, results)
    
}


## another version. coords_norm_sub and coords_norm_out are shared by genes within the same domain, so no need to recalculate them 
# 
powerAnalysisEffectSize2 <- function(shinyDesign, es_range, seq_depth_range, n_rep) {
  
  # === Pre-calculation Phase 1 (Single Pass): Calculate static parameters ===
  
  message("Calculating static simulation parameters...")
  
  # Extract unchanging data once
  count_matrix <- refCounts(shinyDesign)
  loc_data <- refcolData(shinyDesign)[, c('x', 'y', 'domain')]
  par_GP <- paramsGP(shinyDesign)
  all_spot_names <- colnames(count_matrix)
  
  # Pre-calculate normalized coordinates and domain indices once per domain
  coords_norm <- igraph::norm_coords(as.matrix(loc_data[,c('x', 'y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  coords_norm <- as.data.frame(coords_norm)
  coords_norm$domain <- loc_data$domain
  
  domain_static_data <- list()
  for(domain in names(par_GP)){
    idx <- which(coords_norm$domain == domain)
    domain_static_data[[domain]] <- list(
      idx = idx,
      coords_norm_sub = coords_norm[idx, ],
      coords_norm_out = coords_norm[-idx, ]
    )
  }
  
  # Pre-calculate deterministic parameters (b.condition, mu.out, etc.)
  # This part is still done sequentially but is fast due to the small number of genes (~100)
  pre_calc_params <- list()
  for (domain in names(par_GP)) {
    static_data <- domain_static_data[[domain]]
    for (gene in names(par_GP[[domain]])) {
      gene_count_row <- count_matrix[gene, ]
      nnGP_fit <- par_GP[[domain]][[gene]]
      
      # Calculate b.condition once for each gene
      b.condition <- mean_in(SEED = 1, static_data$coords_norm_sub, nnGP_fit)
      mu1 <- sum(b.condition)
      
      # Calculate mu.out once for each gene
      mu.out <- mean_out_refactored(gene_count_row, static_data$idx)
      
      pre_calc_params[[domain]][[gene]] <- list(
        b.condition = b.condition,
        mu1 = mu1,
        mu.out = mu.out
      )
    }
  }
  message("Static parameter calculation complete.")
  
  
  # === Phase 2: Run all simulations in a single parallel loop ===
  
  # Create a grid for all simulation runs (replicates)
  param_grid_stochastic <- expand.grid(
    seq_depth = seq_depth_range,
    effect_size = es_range,
    SEED = seq_len(n_rep)
  )
  
  message("Running simulations and power analysis...")
  
  # Use mclapply to run the full simulation for each replicate
  results <- mclapply(1:nrow(param_grid_stochastic), function(i) {
    row <- param_grid_stochastic[i, ]
    seq_depth_factor <- as.numeric(row["seq_depth"])
    effect_size_factor <- as.numeric(row["effect_size"])
    SEED <- as.integer(row["SEED"])
    
    # Run the simulation for this single parameter combination
    simulated_counts_list <- list()
    for (domain in names(par_GP)) {
      static_data <- domain_static_data[[domain]]
      n1 <- nrow(static_data$coords_norm_sub)
      n2 <- length(all_spot_names) - n1
      
      for (gene in names(par_GP[[domain]])) {
        gene_params <- pre_calc_params[[domain]][[gene]]
        
        # Calculate lambda terms (the only part dependent on ES)
        A <- gene_params$mu1 + n2 * gene_params$mu.out
        lambda.in <- A * effect_size_factor * gene_params$b.condition / (effect_size_factor * gene_params$mu1 + n2 * gene_params$mu.out)
        lambda.out <- gene_params$mu.out * A / (effect_size_factor * gene_params$mu1 + n2 * gene_params$mu.out)
        
        # === Only the rpois call depends on SEED ===
        # The core simulation step is now as fast as it can be
        simulated_counts_list[[domain]][[gene]] <- simulate_genecount_ES_final(
          SEED = SEED,
          seqDepth_factor = seq_depth_factor,
          lambda.in = lambda.in,
          lambda.out = lambda.out,
          coords_norm_sub = static_data$coords_norm_sub,
          coords_norm_out = static_data$coords_norm_out,
          all_spot_names = all_spot_names
        )
      }
    }
    
    # Combine the simulated counts into a single matrix
    all_gene_sims <- do.call(rbind, unlist(simulated_counts_list, recursive = FALSE))
    rownames(all_gene_sims) <- unlist(lapply(names(par_GP), function(domain) {
      paste0(domain, "-", names(par_GP[[domain]]))
    }))
    colnames(all_gene_sims) <- all_spot_names
    
    # Update the shinyDesign object for evaluation
    shinyDesign@simCounts <- all_gene_sims
    shinyDesign@simcolData <- shinyDesign@refcolData
    
    # === The bottleneck is likely here ===
    rst <- evaluatePowerSeurat(shinyDesign)
    NMI <- rst@NMI
    total_counts <- sum(simCounts(shinyDesign))
    
    data.frame(seq_depth = seq_depth_factor, effect_size = effect_size_factor,
               total_counts = total_counts, NMI = NMI)
  }, mc.cores = 4)
  
  results <- do.call(rbind, results)
  return(results)
}


simulate_genecount_ES_final <- function(SEED, seqDepth_factor, lambda.in, lambda.out, coords_norm_sub, coords_norm_out, all_spot_names){
  
  set.seed(SEED)
  y.post <- rpois(n = nrow(coords_norm_sub), lambda = seqDepth_factor * lambda.in)
  
  sim.in <- data.frame(sim = y.post)
  rownames(sim.in) <- rownames(coords_norm_sub)
  
  y.out <- rpois(n = nrow(coords_norm_out), lambda = seqDepth_factor * lambda.out)
  sim.out <- data.frame(sim = y.out)
  rownames(sim.out) <- rownames(coords_norm_out)
  
  count.sim <- rbind(sim.in, sim.out)
  count.sim <- count.sim[order(match(rownames(count.sim), all_spot_names)),, drop = FALSE]
  
  tt <- as.numeric(count.sim$sim)
  return(tt)
}
