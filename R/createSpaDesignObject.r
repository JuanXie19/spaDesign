#' Create a spaDesign object
#' 
#' @description
#' Creates a \code{spaDesign} object to store reference and simulated count data,
#' along with associated spatial metadata and optional model parameters.
#'
#' @param count_matrix A numeric matrix of gene expression count(genes x spots)
#' @param loc A \code{data.frame} containing the spatial coordinates and domain information for spots, 
#'            Required columns: 'x', 'y', and 'domain'.
#'          
#' @return A \code{spaDesign} object 
#' @importFrom methods new validObject
#' @export
#' 
#' @examples
#' \dontrun{
#' sd <- createDesignObject(count_matrix = toyData$toyCount, loc = toyData$loc)
#' }
createDesignObject <- function(count_matrix, loc) {
    # Check count_matrix
    if (!(is.matrix(count_matrix) || inherits(count_matrix, "Matrix"))) {
      stop("count_matrix must be a base matrix or a Matrix object")
    }
    if (!(is.numeric(count_matrix) || inherits(count_matrix, "Matrix"))) {
      stop("count_matrix must contain numeric values or be a sparse Matrix object")
    }
  
    # Check loc
    if (!inherits(loc, "data.frame")) loc <- as.data.frame(loc)
    required_cols <- c("x", "y", "domain")
    missing_cols <- setdiff(required_cols, colnames(loc))
    if (length(missing_cols) > 0) {
      stop("loc is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  
    # Check dimensions
    if (ncol(count_matrix) != nrow(loc)) {
      stop("Number of columns in count_matrix must equal number of rows in loc")
    }
  
    # Check names consistency
    if (!identical(colnames(count_matrix), rownames(loc))) {
      stop("colnames(count_matrix) and rownames(loc) must be identical")
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
        paramsGP = NULL,
        paramsFG = NULL,
		    topGenes = NULL,
        NMI = NULL
    )
    
    ## Validate the object
    if (!validObject(object)) {
        stop("The created spaDesign object is not valid.")
    }
    
    return(object)
}
