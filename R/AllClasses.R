#' @title Class \code{spaDesign}
#' @description The \code{spaDesign} class is used to store and manage data relevant to power analysis within the \code{spaDesign} package.
#' @aliases spaDesign-class
#' @importFrom S4Vectors DataFrame
#' @importFrom methods setClassUnion new
#' @import Matrix
#'
#' @description The \code{spaDesign} class holds data relevant for performing experimental design with the \code{spaDesign} package.
#' It includes slots for storing both reference and simulated count data, along with associated spatial coordinates and domain information, and optional parameters.
#' @details
#' The \code{spaDesign} class includes the following slots:
#' \describe{
#'   \item{\code{refCounts}}{A slot for reference count data}
#'   \item{\code{refcolData}}{A \code{data.frame} containing column-level metadata for reference counts.}
#'   \item{\code{refrowData}}{A \code{data.frame} containing row-level metadata for reference counts.}
#'   \item{\code{simCounts}}{A slot for simulated count data}
#'   \item{\code{simcolData}}{A \code{data.frame} containing column-level metadata for simulated counts.}
#'   \item{\code{simrowData}}{A \code{data.frame} containing row-level metadata for simulated counts.}
#'   \item{\code{paramsPoissonGP}}{An optional list of parameters for Poisson Gaussian processes, or \code{NULL}.}
#'   \item{\code{paramsFG}}{An optional list of parameters for Fisher Gaussian kernal mixture model, or \code{NULL}.}
#'   \item{\code{topGenes}}{An optional list of top genes, or \code{NULL}.}
#'   \item{\code{GP.post}}{An optional list of posterior samples from Gaussian processes, or \code{NULL}.}
#'   \item{\code{NMI}}{A slot for normalized mutual information}
#' }
#'
#' @return The accessor functions \code{refCounts}, \code{\refcolData}, \code{\refrowData}, \cide{\simCounts}, \code{\simcolData}, \code{\simrowData},
#'  \code{paramsPoissonGP}, \code{paramsFG}, \code{\topGenes} and \code{NMI} return the corresponding elements of a \code{spaDesign} object.
#'
#' @export
#'
setClassUnion(name = 'OptionalList', members = c('NULL', 'list'))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

setClass('spaDesign',
    slots = list(
        refCounts = 'ANY',
        refcolData = 'data.frame',
        refrowData = 'data.frame',
        simCounts = 'ANY',
        simcolData = 'data.frame',
        simrowData = 'data.frame',
        paramsPoissonGP = 'OptionalList',
        paramsFG = 'OptionalList',
		topGenes = 'OptionalList',
		GP.post = 'OptionalList',
        NMI = 'ANY'
    ))


