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
    
        shinyDesign <- simulation_EffectSize(shinyDesign, seq_depth_factor, effect_size_factor, SEED)
        
        message(sprintf("Evaluating power for simulated data with seq_depth_factor = %s, effect_size_factor = %s, SEED = %s", 
                        seq_depth_factor, effect_size_factor, SEED))
        
        #rst <- evaluatePower(shinyDesign, conda_env_path)
		rst <- evaluatePowerSeurat(shinyDesign)
        NMI <- rst@NMI
        total_counts <- sum(simCounts(shinyDesign))
        
        data.frame(seq_depth = seq_depth_factor, effect_size = effect_size_factor, 
                   total_counts = total_counts, NMI = NMI)
    }
    
     
  
    results_list <- lapply(1:nrow(param_grid), function(i){
        process_row(param_grid[i, ])
		})
    
    # Combine the results into a single data frame
    results <- do.call(rbind, results_list)
    results.t <- as.data.frame(t(results))
    return(results.t)
}
