#' Access the reference count matrix from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return The reference count matrix
#' @export
setMethod(
    f = 'refCounts',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"refCounts" %in% slotNames(x)) stop("Slot 'refCounts' not found in 'shinyDesign' object")
        x@refCounts
    }
)

#' Access the reference column data from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return The reference column data
#' @export
setMethod(
    f = 'refcolData',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"refcolData" %in% slotNames(x)) stop("Slot 'refcolData' not found in 'shinyDesign' object")
        x@refcolData    
    }
)

#' Access the reference row data from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return The reference row data
#' @export
setMethod(
    f = 'refrowData',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"refrowData" %in% slotNames(x)) stop("Slot 'refrowData' not found in 'shinyDesign' object")
        x@refrowData    
    }
)

#' Access the simulated count matrix from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return The simulated count matrix
#' @export
setMethod(
    f = 'simCounts',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"simCounts" %in% slotNames(x)) stop("Slot 'simCounts' not found in 'shinyDesign' object")
        x@simCounts    
    }
)

#' Access the simulated column data from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return The simulated column data
#' @export
setMethod(
    f = 'simcolData',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"simcolData" %in% slotNames(x)) stop("Slot 'simcolData' not found in 'shinyDesign' object")
        x@simcolData    
    }
)

#' Access fitted Fisher-Gaussian kernel model parameters from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return A list of estimated parameters from fitted Fisher-Gaussian kernel mixture model
#' @export
setMethod(
    f = 'paramsFG',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"paramsFG" %in% slotNames(x)) stop("Slot 'paramsFG' not found in 'shinyDesign' object")
        x@paramsFG    
    }
)

#' Access fittedNNGP model parameters from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return A list of estimated parameters from fitted NNGP model
#' @export
setMethod(
    f = 'paramsGP',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"paramsGP" %in% slotNames(x)) stop("Slot 'paramsGP' not found in 'shinyDesign' object")
        x@paramsGP    
    }
)

#' Access selected top genes for each domain from a shinyDesign object
#' @param x A \code{shinyDesign} object
#' @return A list of selected top genes
#' @export
setMethod(
    f = 'topGenes',
    signature = 'shinyDesign',
    definition = function(x) {
        if (!"topGenes" %in% slotNames(x)) stop("Slot 'topGenes' not found in 'shinyDesign' object")
        x@topGenes    
    }
)

#' Show method for shinyDesign objects
#' @param object A \code{shinyDesign} object
#' @export
setMethod(
  f = 'show',
  signature = 'shinyDesign',
  definition = function(object) {
    cat("class:", class(object), "\n")
    
    if (!is.null(object@refCounts)) {
      cat("Reference dim:", dim(object@refCounts), "\n")
      
      rownames_ref <- ifelse(!is.null(rownames(object@refCounts)), paste(head(rownames(object@refCounts), 5), collapse = ", "), "NULL")
      cat("Reference rownames (first 5):", rownames_ref, ifelse(nrow(object@refCounts) > 5, ", ...", ""), "\n")
      
      colnames_ref <- ifelse(!is.null(colnames(object@refCounts)), paste(head(colnames(object@refCounts), 5), collapse = ", "), "NULL")
      cat("Reference colnames (first 5):", colnames_ref, ifelse(ncol(object@refCounts) > 5, ", ...", ""), "\n")
    } else {
      cat("Reference count matrix: NULL\n")
    }

    if (!is.null(object@refcolData)) {
      cat("refcolData names:", paste(names(object@refcolData), collapse = ", "), "\n")
    } else {
      cat("refcolData: NULL\n")
    }

    if (!is.null(object@simCounts)) {
      cat("Simulated Count dim:", dim(object@simCounts), "\n")
      
      rownames_sim <- ifelse(!is.null(rownames(object@simCounts)), paste(head(rownames(object@simCounts), 5), collapse = ", "), "NULL")
      cat("Simulated Count rownames (first 5):", rownames_sim, ifelse(nrow(object@simCounts) > 5, ", ...", ""), "\n")
      
      colnames_sim <- ifelse(!is.null(colnames(object@simCounts)), paste(head(colnames(object@simCounts), 5), collapse = ", "), "NULL")
      cat("Simulated Count colnames (first 5):", colnames_sim, ifelse(ncol(object@simCounts) > 5, ", ...", ""), "\n")
      
      cat("simcolData names:", paste(names(object@simcolData), collapse = ", "), "\n")
    } else {
      cat("Simulated Count matrix: NULL\n")
    }
  }
)


