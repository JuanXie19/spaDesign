#==================================
# Accessor Generics for shinyDesign
#==================================
#' @title Access slots of a shinyDesign object
#' @description Access reference counts, simulated counts, metadata, model parameters, and top genes
#' stored in a \code{shinyDesign} object.
#' @param x A \code{shinyDesign} object
#' @return Slot content (matrix, data.frame, or list) depending on the accessor
#' @examples
#' \dontrun{
#' sd <- new("shinyDesign",
#'           refCounts = matrix(1:6, nrow = 2),
#'           refcolData = data.frame(sample = 1:3),
#'           refrowData = data.frame(gene = c("A","B")),
#'           simCounts = matrix(7:12, nrow = 2),
#'           simcolData = data.frame(cell = 1:3),
#'           simrowData = data.frame(gene = c("A","B")),
#'           paramsFG = list(param1 = 1),
#'           paramsGP = list(param2 = 2),
#'           topGenes = list(c("gene1","gene2")))
#'
#' refCounts(sd)
#' refcolData(sd)
#' refrowData(sd)
#' simCounts(sd)
#' simcolData(sd)
#' simrowData(sd)
#' paramsFG(sd)
#' paramsGP(sd)
#' topGenes(sd)
#' }
#' @name shinyDesign-accessors
NULL


# -------------------------------
# Reference data accessors
# -------------------------------
#' @rdname shinyDesign-accessors
#' @export
setGeneric("refCounts", function(x) standardGeneric("refCounts"))

#' @rdname shinyDesign-accessors
#' @export
setGeneric("refcolData", function(x) standardGeneric("refcolData"))

#' @rdname shinyDesign-accessors
#' @export
setGeneric("refrowData", function(x) standardGeneric("refrowData"))

# -------------------------------
# Simulated data accessors
# -------------------------------
#' @rdname shinyDesign-accessors
#' @export
setGeneric("simCounts", function(x) standardGeneric("simCounts"))

#' @rdname shinyDesign-accessors
#' @export
setGeneric("simcolData", function(x) standardGeneric("simcolData"))

#' @rdname shinyDesign-accessors
#' @export
setGeneric("simrowData", function(x) standardGeneric("simrowData"))

# -------------------------------
# Model parameter accessors
# -------------------------------
#' @rdname shinyDesign-accessors
#' @export
setGeneric("paramsFG", function(x) standardGeneric("paramsFG"))

#' @rdname shinyDesign-accessors
#' @export
setGeneric("paramsGP", function(x) standardGeneric("paramsGP"))

# -------------------------------
# Top genes
# -------------------------------
#' @rdname shinyDesign-accessors
#' @export
setGeneric("topGenes", function(x) standardGeneric("topGenes"))