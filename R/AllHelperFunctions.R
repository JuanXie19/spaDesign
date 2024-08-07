#' Access the reference count matrix from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return The reference count matrix
#' @export
setMethod(
    f = 'refCounts',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"refCounts" %in% slotNames(x)) stop("Slot 'refCounts' not found in 'spaDesign' object")
        x@refCounts
    }
)

#' Access the reference column data from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return The reference column data
#' @export
setMethod(
    f = 'refcolData',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"refcolData" %in% slotNames(x)) stop("Slot 'refcolData' not found in 'spaDesign' object")
        x@refcolData    
    }
)

#' Access the reference row data from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return The reference row data
#' @export
setMethod(
    f = 'refrowData',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"refrowData" %in% slotNames(x)) stop("Slot 'refrowData' not found in 'spaDesign' object")
        x@refrowData    
    }
)

#' Access the simulated count matrix from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return The simulated count matrix
#' @export
setMethod(
    f = 'simCounts',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"simCounts" %in% slotNames(x)) stop("Slot 'simCounts' not found in 'spaDesign' object")
        x@simCounts    
    }
)

#' Access the simulated column data from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return The simulated column data
#' @export
setMethod(
    f = 'simcolData',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"simcolData" %in% slotNames(x)) stop("Slot 'simcolData' not found in 'spaDesign' object")
        x@simcolData    
    }
)

#' Access fitted Fisher-Gaussian kernel model parameters from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return A list of estimated parameters from fitted Fisher-Gaussian kernel mixture model
#' @export
setMethod(
    f = 'paramsFG',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"paramsFG" %in% slotNames(x)) stop("Slot 'paramsFG' not found in 'spaDesign' object")
        x@paramsFG    
    }
)

#' Access fitted Poisson-Gaussian model parameters from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return A list of estimated parameters from fitted Poisson Gaussian process model
#' @export
setMethod(
    f = 'paramsPoissonGP',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"paramsPoissonGP" %in% slotNames(x)) stop("Slot 'paramsPoissonGP' not found in 'spaDesign' object")
        x@paramsPoissonGP    
    }
)

#' Access selected top genes for each domain from a spaDesign object
#' @param x A \code{spaDesign} object
#' @return A list of selected top genes
#' @export
setMethod(
    f = 'topGenes',
    signature = 'spaDesign',
    definition = function(x) {
        if (!"topGenes" %in% slotNames(x)) stop("Slot 'topGenes' not found in 'spaDesign' object")
        x@topGenes    
    }
)

#' Show method for spaDesign objects
#' @param object A \code{spaDesign} object
#' @export
setMethod(
  f = 'show',
  signature = 'spaDesign',
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


#' Conditional Sampling from Fitted Poisson Gaussian Process Model
#'
#' This function performs conditional sampling from a fitted Poisson Gaussian process model.
#' 
#' @param SEED Random seed for reproducibility
#' @param coords_norm_sub Normalized coordinates for spots within the domain
#' @param alpha.est Estimated alpha parameter
#' @param rho.est Estimated rho parameter
#' @param mu.est Estimated mean parameter
#' @param s.est Estimated s parameter
#' @return A numeric vector representing the sampled expression levels within the domain
#' @export
mean_in <- function(SEED, coords_norm_sub, alpha.est, rho.est, mu.est, s.est) {
    if (missing(SEED) || !is.numeric(SEED)) stop("SEED must be a numeric value.")
    if (!is.data.frame(coords_norm_sub) || !all(c("x", "y") %in% names(coords_norm_sub))) stop("coords_norm_sub must be a data frame with 'x' and 'y' columns.")
    if (!is.numeric(alpha.est) || !is.numeric(rho.est) || !is.numeric(mu.est) || !is.numeric(s.est)) stop("Parameters must be numeric values.")
    
    set.seed(SEED)
    IDX <- sample(seq_len(nrow(coords_norm_sub)), round(nrow(coords_norm_sub) * 0.70))
    data_train <- coords_norm_sub[IDX, ]
    
    dists_train <- as.matrix(dist(cbind(data_train$x, data_train$y)))
    Sigma <- alpha.est^2 * exp(-dists_train / rho.est) + diag(1e-6, nrow(dists_train))
    
    DXX <- as.matrix(dist(cbind(coords_norm_sub$x, coords_norm_sub$y)))
    SXX <- alpha.est^2 * exp(-DXX / rho.est) + diag(1e-6, nrow(coords_norm_sub))
    
    mat1 <- as.matrix(cbind(coords_norm_sub$x, coords_norm_sub$y))
    mat2 <- as.matrix(cbind(data_train$x, data_train$y))
    DX <- as.matrix(pdist(mat1, mat2))
    
    SX <- alpha.est^2 * exp(-DX / rho.est)
    
    Si <- solve(Sigma)
    mup <- mu.est + SX %*% Si %*% (s.est[IDX] - mu.est)
    Sigmap <- SXX - SX %*% Si %*% t(SX)
    Sigmap <- Sigmap + diag(1e-6, nrow(Sigmap))
    
    set.seed(SEED)
    b.condition <- MASS::mvrnorm(1, mup, Sigmap, tol = 1e-3)
    return(b.condition)
}

#' Calculate the Mean for Outside Domain Spots
#'
#' This function calculates the mean of the lower half of counts outside the domain for a given gene.
#' 
#' @param counts A matrix where rows represent genes and columns represent spots.
#' @param gene A character string representing the gene of interest.
#' @param spot_idx Numeric vector indicating the indices of spots within the domain.
#' @return A numeric value representing the mean of the lower half of counts outside the domain for the given gene.
#' @export
mean_out <- function(counts, gene, spot_idx) {
    if (!inherits(counts, "matrix") && !inherits(counts, "CsparseMatrix")) {
      stop("counts must be a matrix or dgCMatrix.")
    }
    if (!is.character(gene) || length(gene) != 1) stop("gene must be a single character string.")
    if (!is.numeric(spot_idx)) stop("spot_idx must be numeric.")
    
    gene.count <- counts[which(rownames(counts) == gene), ]
    OUT <- gene.count[-spot_idx]
    OUT.lower <- OUT[OUT <= median(OUT)]
    mean.out.low <- mean(OUT.lower)
    return(mean.out.low)
}

