#' @title Accessors for \code{shinyDesign} objects
#' 
#' @description Access slots of a \code{shinyDesign} object, including reference
#' and simulated data, model parameters, and results.
#' 
#' @param x A \code{shinyDesign} object.
#' 
#' @return Slot contents:
#' \itemize{
#'   \item \code{refCounts} – reference count matrix
#'   \item \code{refcolData} – reference column metadata
#'   \item \code{refrowData} – reference row metadata
#'   \item \code{simCounts} – simulated count matrix
#'   \item \code{simcolData} – simulated column metadata
#'   \item \code{simrowData} – simulated row metadata
#'   \item \code{paramsFG} – Fisher–Gaussian mixture model parameters
#'   \item \code{paramsGP} – Gaussian process/NNGP model parameters
#'   \item \code{topGenes} – list of top selected genes
#' }
#' 
#' @examples
#' \dontrun{
#' obj <- new("shinyDesign",
#'            refCounts = matrix(1:6, nrow = 2),
#'            refcolData = data.frame(sample = 1:3),
#'            refrowData = data.frame(gene = c("A", "B")),
#'            simCounts = matrix(7:12, nrow = 2),
#'            simcolData = data.frame(cell = 1:3),
#'            simrowData = data.frame(gene = c("A", "B")))
#' 
#' refCounts(obj)
#' simCounts(obj)
#' paramsFG(obj)
#' }
#' 
#' @name shinyDesign-accessors
NULL


# Reference data
#' @rdname shinyDesign-accessors
#' @export
setMethod("refCounts", "shinyDesign", function(x) x@refCounts)

#' @rdname shinyDesign-accessors
#' @export
setMethod("refcolData", "shinyDesign", function(x) x@refcolData)


#' @rdname shinyDesign-accessors
#' @export
setMethod("refrowData", "shinyDesign", function(x) x@refrowData)

# Simulated data
#' @rdname shinyDesign-accessors
#' @export
setMethod("simCounts", "shinyDesign", function(x) x@simCounts)

#' @rdname shinyDesign-accessors
#' @export
setMethod("simcolData", "shinyDesign", function(x) x@simcolData)

#' @rdname shinyDesign-accessors
#' @export
setMethod("simrowData", "shinyDesign", function(x) x@simrowData)

# Parameters
#' @rdname shinyDesign-accessors
#' @export
setMethod("paramsFG", "shinyDesign", function(x) x@paramsFG)

#' @rdname shinyDesign-accessors
#' @export
setMethod("paramsGP", "shinyDesign", function(x) x@paramsGP)

# Top genes
#' @rdname shinyDesign-accessors
#' @export
setMethod("topGenes", "shinyDesign", function(x) x@topGenes)

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


