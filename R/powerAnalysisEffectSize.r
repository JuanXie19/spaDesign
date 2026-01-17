#' Power analysis via effect size and sequencing depth simulations
#'
#' Performs power analysis by simulating data under varying effect sizes and sequencing depths,
#' then evaluating clustering performance.
#'
#' @param spaDesign A \code{spaDesign} object with estimated GP and FG model parameters.
#' @param es_range Numeric vector of effect size scaling factors.
#' @param seq_depth_range Numeric vector of sequencing depth scaling factors.
#' @param n_rep Integer. Number of replications per scenario.
#' @param n_cores Integer. Number of CPU cores for parallelization (default: 2).
#'
#' @details
#' For each combination of effect size and sequencing depth, \code{n_rep} datasets are simulated
#' and evaluated. Clustering performance is measured using Normalized Mutual Information (NMI)
#' between predicted and true domain labels.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{seq_depth}{Sequencing depth scaling factor}
#'     \item{effect_size}{Effect size scaling factor}
#'     \item{total_counts}{Total counts in simulated dataset}
#'     \item{NMI}{Normalized Mutual Information between predicted and true labels}
#'   }
#'
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#' \dontrun{
#' results <- powerAnalysisEffectSize(
#'   spaDesign = my_spadesign,
#'   es_range = c(1, 2, 3),
#'   seq_depth_range = c(0.5, 1, 3, 5),
#'   n_rep = 5,
#'   n_cores = 4
#' )
#' head(results)
#' }


powerAnalysisEffectSize <- function(spaDesign, es_range, seq_depth_range, n_rep, n_cores) {
    # Input validation
    if (!inherits(spaDesign, "spaDesign")) {
      stop("spaDesign must be a spaDesign object")
    }
  
    if (!is.numeric(es_range) || any(es_range <= 0)) {
      stop("es_range must contain positive numeric values")
    }
  
    if (!is.numeric(seq_depth_range) || any(seq_depth_range <= 0)) {
      stop("seq_depth_range must contain positive numeric values")
    }
  
    if (!is.numeric(n_rep) || n_rep <= 0 || n_rep != round(n_rep)) {
      stop("n_rep must be a positive integer")
    }
  
    if (!is.numeric(n_cores) || n_cores <= 0 || n_cores != round(n_cores)) {
      stop("n_cores must be a positive integer")
    }
  
  
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
    
        spaDesign <- simulation_EffectSize_refactored(spaDesign, seq_depth_factor, effect_size_factor, SEED)
        
        message(sprintf("Evaluating power for simulated data with seq_depth_factor = %s, effect_size_factor = %s, SEED = %s", 
                        seq_depth_factor, effect_size_factor, SEED))
        
		    rst <- evaluatePowerSeurat(spaDesign)

        data.frame(
          seq_depth = seq_depth_factor, 
          effect_size = effect_size_factor, 
          total_counts = sum(spaDesign@simCounts), 
          NMI = rst@NMI)
    }
  
    results <- mclapply(1:nrow(param_grid), function(i){
      process_row(param_grid[i, ])
    }, mc.cores = n_cores)
    
    # Check for errors
    if (any(sapply(results, inherits, "try-error"))) {
      stop("Error in parallel processing. Check individual scenario results.")
    }
    
    results <- do.call(rbind, results)
    
}


