#' @title Class \code{shinyDesign}
#' @description The \code{shinyDesign} class is used to store and manage data relevant to power analysis within the \code{shinyDesign} package.
#' @aliases shinyDesign-class
#' @importFrom S4Vectors DataFrame
#' @importFrom methods setClassUnion new
#' @import Matrix
#'
#' @description The \code{shinyDesign} class holds data relevant for performing experimental design with the \code{shinyDesign} package.
#' It includes slots for storing both reference and simulated count data, along with associated spatial coordinates and domain information, and optional parameters.
#' @details
#' The \code{shinyDesign} class includes the following slots:
#' \describe{
#'   \item{\code{refCounts}}{A slot for reference count data}
#'   \item{\code{refcolData}}{A \code{data.frame} containing column-level metadata for reference counts.}
#'   \item{\code{refrowData}}{A \code{data.frame} containing row-level metadata for reference counts.}
#'   \item{\code{simCounts}}{A slot for simulated count data}
#'   \item{\code{simcolData}}{A \code{data.frame} containing column-level metadata for simulated counts.}
#'   \item{\code{simrowData}}{A \code{data.frame} containing row-level metadata for simulated counts.}
#'   \item{\code{paramsGP}}{An optional list of parameters for Gaussian processes, or \code{NULL}.}
#'   \item{\code{paramsFG}}{An optional list of parameters for Fisher Gaussian kernal mixture model, or \code{NULL}.}
#'   \item{\code{topGenes}}{An optional list of top genes, or \code{NULL}.}
#'   \item{\code{NMI}}{A slot for normalized mutual information}
#' }
#'
#' @return The accessor functions \code{refCounts}, \code{\refcolData}, \code{\refrowData}, \cide{\simCounts}, \code{\simcolData}, \code{\simrowData},
#'  \code{paramsGP}, \code{paramsFG}, \code{\topGenes} and \code{NMI} return the corresponding elements of a \code{shinyDesign} object.
#'
#' @export
#'
setClassUnion(name = 'OptionalList', members = c('NULL', 'list'))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

setClass('shinyDesign',
    slots = list(
        refCounts = 'ANY',
        refcolData = 'data.frame',
        refrowData = 'data.frame',
        simCounts = 'ANY',
        simcolData = 'data.frame',
        simrowData = 'data.frame',
        paramsGP = 'OptionalList',
        paramsFG = 'OptionalList',
		topGenes = 'OptionalList',
        NMI = 'ANY'
    ))


