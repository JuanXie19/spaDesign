#' @title Extract the reference count matrix
#' @param x shinyDesign object
#' @export
#' @return Returns a reference count matrix
#' @examples
#' refCounts(shinyDesign)
setGeneric(
    name = 'refCounts',
    signature = 'x',
    def = function(x){
        standardGeneric('refCounts')    
    }
)

#' @title Extract the reference colData
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns the colData of reference data
#' @examples
#' refcolData(shinyDesign)
setGeneric(
    name = 'refcolData',
    signature = 'x',
    def = function(x){
        standardGeneric('refcolData')
    }
)

#' @title Access reference rowData 
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns the rowData of reference data
#' @examples
#' refrowData(shinyDesign)
setGeneric(
    name = 'refrowData',
    signature = 'x',
    def = function(x){
        standardGeneric('refrowData')
    }
)

#' @title Access simulated count matrix
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns the simulated count matrix
#' @examples
#' simCounts(shinyDesign)
setGeneric(
    name = 'simCounts',
    signature = 'x',
    def = function(x){
        standardGeneric('simCounts')
    }
)

#' @title Access colData for simulated data 
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns the colData of simulated data
#' @examples
#' simcolData(shinyDesign)
setGeneric(
    name = 'simcolData',
    signature = 'x',
    def = function(x){
        standardGeneric('simcolData')
    }
)

#' @title Access rowData for simulated data 
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns the rowData of simulated data
#' @examples
#' simrowData(shinyDesign)
setGeneric(
    name = 'simrowData',
    signature = 'x',
    def = function(x){
        standardGeneric('simrowData')
    }
)

#' @title Access estimated parameter from the fitted Fisher-Gaussian kernel model
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns a list of estimated parameters from fitted Fisher-Gaussian kernel mixture model
#' @examples
#' paramsFG(shinyDesign)
setGeneric(
    name = 'paramsFG',
    signature = 'x',
    def = function(x){
        standardGeneric('paramsFG')
    }
)

#' @title Access estimated parameter from the fitted Gaussian process model
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns a list of estimated parameters from fitted Gaussian process model
#' @examples
#' paramsGP(shinyDesign)
setGeneric(
    name = 'paramsGP',
    signature = 'x',
    def = function(x){
        standardGeneric('paramsGP')
    }
)

#' @title Access selected top genes for each domain
#' @param x \code{shinyDesign} object
#' @export
#' @return Returns a list of selected top genes for each domain
#' @examples
#' topGenes(shinyDesign)
setGeneric(
    name = 'topGenes',
    signature = 'x',
    def = function(x){
        standardGeneric('topGenes')
    }
)

