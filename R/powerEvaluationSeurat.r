#' Evaluate clustering performance on simulated spatial transcriptomics data using Seurat
#'
#' This function applies a Seurat-based clustering pipeline to a \code{spaDesign} object
#' with simulated counts. Clustering performance is evaluated against the true domains
#' using Normalized Mutual Information (NMI).
#'
#' @param spaDesign A \code{spaDesign} object containing:
#'   \itemize{
#'     \item \code{simCounts}: simulated count matrix (genes x spots).
#'     \item \code{simcolData}: spot metadata with domain assignments.
#'   }
#' @return The updated \code{spaDesign} object with:
#'   \itemize{
#'     \item \code{simcolData$label_pred}: predicted cluster labels.
#'     \item \code{NMI}: numeric NMI between true and predicted labels.
#'   }
#'
#' @import aricode
#' @import Seurat
#' @examples
#' # Assuming spaDesign is a valid spaDesign object
#' # result <- evaluatePower(spaDesign)
#' @export
#'
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
