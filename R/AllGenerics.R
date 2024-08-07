#' @title Extract the reference count matrix
#' @param x spaDesign object
#' @export
#' @return Returns a reference count matrix
#' @examples
#' refCounts(spaDesign)
setGeneric(
    name = 'refCounts',
    signature = 'x',
    def = function(x){
        standardGeneric('refCounts')    
    }
)

#' @title Extract the reference colData
#' @param x \code{spaDesign} object
#' @export
#' @return Returns the colData of reference data
#' @examples
#' refcolData(spaDesign)
setGeneric(
    name = 'refcolData',
    signature = 'x',
    def = function(x){
        standardGeneric('refcolData')
    }
)

#' @title Access reference rowData 
#' @param x \code{spaDesign} object
#' @export
#' @return Returns the rowData of reference data
#' @examples
#' refrowData(spaDesign)
setGeneric(
    name = 'refrowData',
    signature = 'x',
    def = function(x){
        standardGeneric('refrowData')
    }
)

#' @title Access simulated count matrix
#' @param x \code{spaDesign} object
#' @export
#' @return Returns the simulated count matrix
#' @examples
#' simCounts(spaDesign)
setGeneric(
    name = 'simCounts',
    signature = 'x',
    def = function(x){
        standardGeneric('simCounts')
    }
)

#' @title Access colData for simulated data 
#' @param x \code{spaDesign} object
#' @export
#' @return Returns the colData of simulated data
#' @examples
#' simcolData(spaDesign)
setGeneric(
    name = 'simcolData',
    signature = 'x',
    def = function(x){
        standardGeneric('simcolData')
    }
)

#' @title Access rowData for simulated data 
#' @param x \code{spaDesign} object
#' @export
#' @return Returns the rowData of simulated data
#' @examples
#' simrowData(spaDesign)
setGeneric(
    name = 'simrowData',
    signature = 'x',
    def = function(x){
        standardGeneric('simrowData')
    }
)

#' @title Access estimated parameter from the fitted Fisher-Gaussian kernel model
#' @param x \code{spaDesign} object
#' @export
#' @return Returns a list of estimated parameters from fitted Fisher-Gaussian kernel mixture model
#' @examples
#' paramsFG(spaDesign)
setGeneric(
    name = 'paramsFG',
    signature = 'x',
    def = function(x){
        standardGeneric('paramsFG')
    }
)

#' @title Access estimated parameter from the fitted Poisson Gaussian process model
#' @param x \code{spaDesign} object
#' @export
#' @return Returns a list of estimated parameters from fitted Poisson Gaussian process model
#' @examples
#' paramsPoissonGP(spaDesign)
setGeneric(
    name = 'paramsPoissonGP',
    signature = 'x',
    def = function(x){
        standardGeneric('paramsPoissonGP')
    }
)

#' @title Access selected top genes for each domain
#' @param x \code{spaDesign} object
#' @export
#' @return Returns a list of selected top genes for each domain
#' @examples
#' topGenes(spaDesign)
setGeneric(
    name = 'topGenes',
    signature = 'x',
    def = function(x){
        standardGeneric('topGenes')
    }
)

