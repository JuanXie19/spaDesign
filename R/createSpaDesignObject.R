#' Create spaDesign object
#' @param count_matrix A gene expression count \code{matrix}
#' @param loc A \code{data.frame} containing the spatial coordinates and domain information for spots, corresponding to columns named as 'x', 'y', and 'domain'
#' @return Returns a \code{spaDesign} object 
#' @importFrom methods new validObject
#' @export 
#' @examples
#'
#' ## Create a spaDesign object
#' toyDATA  <- createSpaDesignObject(count_matrix = toyData$toyCount, loc = toyData$loc)
#'
createSpaDesignObject <- function(count_matrix, loc) {
    ## Check dimensions
    if (ncol(count_matrix) != nrow(loc)) {
        stop("The number of spots in count_matrix and loc should be consistent!")
    }
    
    ## Check location dimension
    required_columns <- c("x", "y", "domain")
    missing_columns <- setdiff(required_columns, colnames(loc))
    if (length(missing_columns) > 0) {
        stop(paste("loc file must contain these three columns: x, y and domain. The following columns are/is missing:", paste(missing_columns, collapse = ", ")))
    }
    
    ## Check location index consistency
    if (!identical(colnames(count_matrix), rownames(loc))) {
        stop("colnames in the count_matrix file is different from the rownames in the loc file")
    }
    
    
    ## Convert loc to data.frame if not already
    if (!inherits(loc, "data.frame")) {
        loc <- as.data.frame(loc)
    }
    
    ## Create the spaDesign object
    object <- new(
        Class = "spaDesign",
        refCounts = count_matrix,
        refcolData = loc,
        refrowData = loc,  # Assuming similar structure for rowData; adjust if needed
        simCounts = matrix(numeric(), nrow = 0, ncol = 0),  # Empty matrix; adjust if needed
        simcolData = data.frame(),  # Initialize with empty data.frame
        simrowData = data.frame(),  # Initialize with empty data.frame
        paramsPoissonGP = NULL,
        paramsFG = NULL,
		topGenes = NULL,
		GP.post = NULL,
        NMI = NULL
    )
    
    ## Validate the object
    if (!validObject(object)) {
        stop("The created spaDesign object is not valid.")
    }
    
    return(object)
}
