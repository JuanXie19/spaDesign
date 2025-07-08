#' Apply spaGCN to the simulated dataset to obtain predicted labels
#' @param shinyDesign A \code{shinyDesign} object with simulated count matrix
#' @return The updated \code{shinyDesign} object with predicted labels and NMI
#' @import pbapply
#' @import readr
#' @import aricode
#' @examples
#' # Assuming shinyDesign is a valid shinyDesign object
#' # result <- evaluatePower(shinyDesign)
#' @export
#'
evaluatePower <- function(shinyDesign, conda_env_path) {
    library(pbapply)
    library(readr)
    library(reticulate)
    python_executable <- file.path(conda_env_path, "bin/python")
    Sys.setenv(RETICULATE_PYTHON = python_executable)
    #reticulate::use_condaenv('SpaGCN', required = TRUE)
    reticulate::py_config()
    
    # locate the python script
    python_script <- system.file('python','runSpaGCN_singleSim.py', package = 'shinyDesign')
    if (!file.exists(python_script)) {
        stop("Python script not found: ", python_script)
    }
    
    
    # Check if simCounts and simcolData are available
    if (is.null(simCounts(shinyDesign))) {
        stop("simCounts is not available in the provided object. Please run data simulation first.")
    }
  
    if (is.null(simcolData(shinyDesign))) {
        stop("simcolData is not available in the provided object. Please double check!")
    }
    
    
    counts <- simCounts(shinyDesign)
    loc <- simcolData(shinyDesign)[, c('x', 'y', 'domain')]

    count.t <- t(counts)
    count.t <- as.data.frame(count.t)

    temp_file1 <- tempfile(fileext = '.tsv')
    write_tsv(count.t, temp_file1)

    spatial <- loc
    colnames(spatial) <- c('x', 'y', 'regions')
    temp_file2 <- tempfile(fileext = '.csv')
    write.csv(spatial, temp_file2, row.names = FALSE)

    output_file <- tempfile(fileext = '.csv')
    
    command <- sprintf("%s %s %s %s %s", python_executable, python_script, temp_file1, temp_file2, output_file)

    message('Running spaGCN...')
    #message('Running spaGCN with command:\n', command)
    python_output <- tryCatch({
        system(command, intern = TRUE)
    }, error = function(e) {
        stop("Failed to run Python script: ", e$message)
    })

    if (length(python_output) > 0) {
        warning("Python script output:\n", paste(python_output, collapse = "\n"))
    }

    results <- tryCatch({
        read.csv(output_file)
    }, error = function(e) {
        stop("Failed to read output file: ", e$message)
    })

    if (!all(c('cluster', 'pred') %in% colnames(results))) {
        stop("Results do not contain the required columns 'cluster' and 'pred'")
    }

    message('Computing NMI...')
    nmi <- NMI(results$cluster, results$pred, variant = "sqrt")
    shinyDesign@simcolData$label_pred <- results[,'pred']
    shinyDesign@NMI <- nmi

    file.remove(temp_file1)
    file.remove(temp_file2)
    file.remove(output_file)

    message('Performance evaluation complete. NMI: ', round(nmi,3))
    return(shinyDesign)
}
