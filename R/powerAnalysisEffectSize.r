#' Power analysis via effect size and sequencing depth simulations
#'
#' This function performs power analysis by simulating data under varying
#' effect sizes and sequencing depths, followed by clustering evaluation.
#'
#' @param shinyDesign A \code{shinyDesign} object with estimated model parameters.
#' @param es_range A numeric vector of effect size scaling factors
#' @param seq_depth_range A numeric vector of sequencing depth scaling factors.
#' @param n_rep Integer, number of replications for each scenario.
#' @param n_cores Integer, number of CPU cores for parallelization (default: 2)
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{seq_depth}: Sequencing depth scaling factor
#'     \item \code{effect_size}: Effect size scaling factor
#'     \item \code{total_counts}: Total counts in simulated dataset
#'     \item \code{NMI}: Normalized Mutual Information between predicted and true labels
#'   }
#'   
#' @import pbmcapply
#' @examples
#' # spadesign <- ...
#' # es_range <- c(1, 2, 3)
#' # seq_depth_range <- c(0.5, 1, 3, 5)
#' # n_rep <- 10
#' # results <- powerAnalysisEffectSize(spadesign, es_range, seq_depth_range, n_rep, n_cores)
#' @export


powerAnalysisEffectSize <- function(shinyDesign, es_range, seq_depth_range, n_rep, n_cores) {
    # Create a data frame of all combinations of es_range, seq_depth_range, and n_rep
    param_grid <- expand.grid(seq_depth = seq_depth_range, 
                              effect_size = es_range, 
                              SEED = seq_len(n_rep))

    # Function to process each row of the param_grid
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
        
		    rst <- evaluatePowerSeurat(shinyDesign)

        data.frame(
          seq_depth = seq_depth_factor, 
          effect_size = effect_size_factor, 
          total_counts = sum(shinyDesign@simCounts), 
          NMI = rst@NMI)
    }
    
    
    
    results <- mclapply(1:nrow(param_grid), function(i){
      process_row(param_grid[i, ])
    }, mc.cores = n_cores)
    results <- do.call(rbind, results)
    
}


