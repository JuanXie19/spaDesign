#' Power analysis for spatial data
#' 
#' This function performs power analysis by simulating spatial data under
#'  different proportions of undisturbed patterns and sequencing depths.
#' 
#' @param shinyDesign A \code{shinyDesign} object with estimated model parameters.
#' @param SIGMA A numeric value indicating the Gaussian noise level.
#' @param prop_range A numeric vector indicating the range of proportions 
#'  (proportion of genes with undisturbed pattern).
#' @param seq_depth_range A numeric vector indicating the range of sequencing depth scaling factors.
#' @param n_rep An integer indicating the number of replications for each scenario.
#' @param n_cores Number of COP cores to use for parallellization (default = 1).
#' 
#' @return A data frame with columns:
#' \describe{
#'   \item{seq_depth}{Sequencing depth scaling factor.}
#'   \item{prop}{Proportion of undisturbed patterns.}
#'   \item{total_counts}{Total counts in the simulated data.}
#'   \item{NMI}{Normalized Mutual Information from power evaluation.}
#' }
#' 
#' @import pbmcapply
#' @export
#' 
#' @examples
#' \dontrun{
#'   prop_range <- c(0.1, 0.5, 1)
#'   seq_depth_range <- c(0.5, 1, 3, 5)
#'   n_rep <- 10
#'   results <- powerAnalysisSpatial(shinyDesign, SIGMA = 1, prop_range, seq_depth_range, n_rep, n_cores = 2)
#' }

powerAnalysisSpatial <- function(shinyDesign, SIGMA, prop_range, seq_depth_range, n_rep, n_cores) {
    
    # Input checks
    stopifnot(is.numeric(prop_range), all(prop_range >= 0), all(prop_range <= 1))
    stopifnot(is.numeric(seq_depth_range), all(seq_depth_range > 0))
    stopifnot(is.numeric(n_rep), length(n_rep) == 1, n_rep > 0)
    
    # Create a data frame of all combinations of es_range, seq_depth_range, and n_rep
    param_grid <- expand.grid(seq_depth = seq_depth_range, 
                              prop = prop_range, 
                              SEED = seq_len(n_rep))

    # Define a function to process each row of the param_grid
    process_row <- function(row) {
        prop <- as.numeric(row['prop'])
        seq_depth_factor <- as.numeric(row["seq_depth"])
        SEED <- as.numeric(row["SEED"])
		

        message(sprintf("Simulating data with seq_depth_factor = %s, prop = %s, SEED = %s", 
                        seq_depth_factor, prop, SEED))
        DATA <- simulation_Spatial(shinyDesign,selected_M_list = NULL, seq_depth_factor, SIGMA, SEED, prop, n_cores = n_cores)
        
        message(sprintf("Evaluating power for simulated data with seq_depth_factor = %s, prop = %s, SEED = %s", 
                        seq_depth_factor, prop, SEED))
                        
		    rst <- evaluatePowerSeurat(DATA)
        NMI <- rst@NMI
        total_counts <- sum(simCounts(DATA))
        
        data.frame(seq_depth = seq_depth_factor, 
                   prop = prop, 
                   total_counts = total_counts, 
                   NMI = NMI)
    }
  
    results <- mclapply(1:nrow(param_grid), function(i){
      process_row(param_grid[i, ])
    }, mc.cores = n_cores)
    results <- do.call(rbind, results)
    return(results)
}
