#' Power analysis for spatial data
#' 
#' This function performs power analysis by simulating spatial data with different proportions of undisturbed patterns and sequencing depths.
#' 
#' @param shinyDesign A \code{shinyDesign} object with estimated model parameters
#' @param prop_range A numeric vector indicating the range of proportions (proportion of genes with undisturbed pattern)
#' @param seq_depth_range A numeric vector indicating the range of sequencing depth scaling factors
#' @param n_rep An integer indicating the number of replications for each scenario
#' @param conda_env_path Specify the path to the conda environment for spaGCN
#' @param python_script Specify the path to the python script
#' @return A data frame containing the sequencing depth, proportion, total counts, and NMI
#' @name powerAnalysisSpatial
#' @examples
#' # spadesign <- ...
#' # prop_range <- c(0.1, 0.5, 1)
#' # seq_depth_range <- c(0.5, 1, 3, 5)
#' # n_rep <- 10
#' # results <- powerAnalysisSpatial(shinyDesign, prop_range, seq_depth_range, n_rep, conda_env_path) 
#' @export

library(parallel)

powerAnalysisSpatial <- function(shinyDesign, prop_range, seq_depth_range, n_rep) {
    # Create a data frame of all combinations of es_range, seq_depth_range, and n_rep
    param_grid <- expand.grid(seq_depth = seq_depth_range, 
                              prop = prop_range, 
                              SEED = seq_len(n_rep))

    # Define a function to process each row of the param_grid
    process_row <- function(row) {
        prop <- as.numeric(row['prop'])
        seq_depth_factor <- as.numeric(row["seq_depth"])
        SEED <- as.numeric(row["SEED"])
		
		# current just allow SIGMA = 1
		SIGMA <- 1
        
        message(sprintf("Simulating data with seq_depth_factor = %s, prop = %s, SEED = %s", 
                        seq_depth_factor, prop, SEED))
        DATA <- simulation_Spatial(shinyDesign, selected_M_list = NULL, seq_depth_factor, SIGMA, SEED, prop)
        
        message(sprintf("Evaluating power for simulated data with seq_depth_factor = %s, prop = %s, SEED = %s", 
                        seq_depth_factor, prop, SEED))
                        
        #rst <- evaluatePower(DATA, conda_env_path)
		rst <- evaluatePowerSeurat(DATA)
        NMI <- rst@NMI
        total_counts <- sum(simCounts(DATA))
        
        data.frame(seq_depth = seq_depth_factor, prop = prop, 
                   total_counts = total_counts, NMI = NMI)
    }
    
    
    
    #num_cores <- min(detectCores()-2, 8)
    #cl <- makeCluster(num_cores)
    #on.exit(stopCluster(cl))
	#registerDoParallel(cl)
  
   # results_list <- foreach(i = 1:nrow(param_grid), .combine = rbind, .packages = 'shinyDesign2') %dopar% {
   #     process_row(param_grid[i, ])
   # }
	
	results_list <- lapply(1:nrow(param_grid), function(i){
		process_row(param_grid[i, ]
		)
	})
	
    
    # Combine the results into a single data frame
    results <- do.call(rbind, results_list)
    results.t <- as.data.frame(t(results))
    return(results.t)
}
