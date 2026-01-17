#' Evaluate clustering performance on simulated spatial transcriptomics data using Seurat
#'
#' This function applies a Seurat-based clustering pipeline to the simulated
#'  count matrix stored in a \code{spaDesign} object
#'  and evaluate clustering accuracy against the ground-truth domain labels. Performance
#'  is quantified using Normalized Mutual Information (NMI).
#'
#' @param spaDesign A \code{spaDesign} object containing:
#'   \itemize{
#'     \item \code{simCounts}: simulated count matrix (genes x spots).
#'     \item \code{simcolData}: spot-level metadata with true domain labels.
#'   }
#' @return The updated \code{spaDesign} object with:
#'   \itemize{
#'     \item \code{simcolData$label_pred}: predicted cluster labels for each spot.
#'     \item \code{NMI}: numeric NMI between true and predicted labels.
#'   }
#'
#' @import aricode
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData RunPCA
#' @importFrom Seurat FindNeighbors FindClusters RenameIdents Idents AggregateEexpression
#' @noRd
#' 
#' @examples
#' \dontrun{
#' # spaDesign must already contain simulated counts and truth labels:
#' # simCounts(spaDesign) : genes x spots matrix
#' # simcolData(spaDesign)$domain : true domain assignment per spot
#'
#' out <- evaluatePowerSeurat(spaDesign)
#' out@NMI
#' head(out@simcolData$label_pred)
#' }
evaluatePowerSeurat <- function(spaDesign) {
   
    # Check if simCounts and simcolData are available
    if (is.null(simCounts(spaDesign))) {
        stop("simCounts is not available in the provided object. Please run data simulation first.")
    }
  
    if (is.null(simcolData(spaDesign))) {
        stop("simcolData is not available in the provided object. Please double check!")
    }
    
    
    counts <- simCounts(spaDesign)
    n.cluster <- length(unique(simcolData(spaDesign)$domain))
	
	  seurat <- CreateSeuratObject(counts = counts, assay = 'RNA')
	  seurat$true_domain <- simcolData(spaDesign)$domain
	
	  seurat <- NormalizeData(seurat)
	  seurat <- FindVariableFeatures(seurat, selection.method = 'vst')
	  seurat <- ScaleData(seurat, features = rownames(seurat))
	  seurat <- RunPCA(seurat, verbose = FALSE)
	  seurat <- FindNeighbors(seurat, dims = 1:20, verbose = FALSE)
	  seurat <- FindClusters(seurat, verbose = FALSE,res = 2 )

	  X <- Seurat::AggregateExpression(seurat, 
	                                 assays=SeuratObject::DefaultAssay(seurat),  
					                          slot= "scale.data", 
					                          group.by = "seurat_clusters")[[1]] 

	  clust2 <- cutree(hclust(dist(t(X))), k = n.cluster)
	  new_labels <- clust2
	  names(new_labels) <- levels(seurat)
	  seurat <- RenameIdents(seurat,new_labels)
	  seurat$label_pred <- Idents(seurat)
	
    nmi <- aricode::NMI(seurat$true_domain, seurat$label_pred, variant = "sqrt")
    spaDesign@simcolData$label_pred <- seurat$label_pred
    spaDesign@NMI <- nmi

    message('Performance evaluation complete. NMI: ', round(nmi,3))
    return(spaDesign)
}
