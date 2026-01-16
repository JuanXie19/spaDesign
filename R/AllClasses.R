#' @title Class \code{spaDesign}
#' 
#' @description 
#' The \code{spaDesign} class is used to store and manage data relevant to power analysis 
#' within the \code{spaDesign} package. It holds both reference and simulated count data, associated meta data,
#' and optional parameters from fitted models.
#' 
#' @slot refCounts A matrix-like object (e.g., \code{Matrix::dgCMatrix}) of reference counts.
#' @slot refcolData A \code{data.frame} containing column-level metadata for reference counts.
#' @slot refrowData A \code{data.frame} containing row-level metadata for reference counts.
#' @slot simCounts A matrix-like object of simulated counts.
#' @slot simcolData A \code{data.frame} containing column-level metadata for simulated counts.
#' @slot simrowData A \code{data.frame} containing row-level metadata for simulated counts.
#' @slot paramsGP An optional list of parameters for Gaussian processes, or \code{NULL}.
#' @slot paramsFG An optional list of parameters for Fisher Gaussian kernel mixture model, or \code{NULL}.
#' @slot selected_M_list_AIC Best model selection results based on AIC.
#' @slot selected_M_list_BIC Best model selection results based on BIC.
#' @slot topGenes An optional list of top genes, or \code{NULL}.
#' @slot NMI Normalized mutual information value (numeric or object).
#' 
#' @details
#' The \code{spaDesign} class includes the following slots:
#' \describe{
#'   \item{\code{refCounts}}{A slot for reference count data}
#'   \item{\code{refcolData}}{A \code{data.frame} containing column-level metadata for reference counts.}
#'   \item{\code{refrowData}}{A \code{data.frame} containing row-level metadata for reference counts.}
#'   \item{\code{simCounts}}{A slot for simulated count data}
#'   \item{\code{simcolData}}{A \code{data.frame} containing column-level metadata for simulated counts.}
#'   \item{\code{simrowData}}{A \code{data.frame} containing row-level metadata for simulated counts.}
#'   \item{\code{paramsGP}}{An optional list of parameters for Gaussian processes, or \code{NULL}.}
#'   \item{\code{paramsFG}}{An optional list of parameters for Fisher Gaussian kernal mixture model, or \code{NULL}.}
#'	 \item{\code{selected_M_list_AIC}}{A slot for selected best number of cluster for FG mixture model based on AIC criteria.}
#'	 \item{\code{selected_M_list_BIC}}{A slot for selected best number of cluster for FG mixture model based on BIC criteria.}
#'   \item{\code{topGenes}}{An optional list of top genes, or \code{NULL}.}
#'   \item{\code{NMI}}{A slot for normalized mutual information}
#' }
#'
#' @details 
#' This class provides a container for storing both observed and simulated 
#' spatial transcriptomics data, along with results from downstream modeling.
#' 
#' @return A \code{spaDesign} object.
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom methods setClassUnion setClass
#' @importFrom Matrix Matrix
#' @export
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
        paramsGP = 'OptionalList',
        paramsFG = 'OptionalList',
		selected_M_list_AIC = 'ANY',
		selected_M_list_BIC = 'ANY',
		topGenes = 'OptionalList',
        NMI = 'ANY'
    ))


