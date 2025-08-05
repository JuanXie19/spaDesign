#' Apply spaGCN to the simulated dataset to obtain predicted labels
#' @param shinyDesign A \code{shinyDesign} object with simulated count matrix
#' @return The updated \code{shinyDesign} object with predicted labels and NMI
#' @import pbapply
#' @import readr
#' @import aricode
#' @import Seurat
#' @examples
#' # Assuming shinyDesign is a valid shinyDesign object
#' # result <- evaluatePower(shinyDesign)
#' @export
#'
evaluatePowerSeurat <- function(shinyDesign) {
    library(pbapply)
    
    # Check if simCounts and simcolData are available
    if (is.null(simCounts(shinyDesign))) {
        stop("simCounts is not available in the provided object. Please run data simulation first.")
    }
  
    if (is.null(simcolData(shinyDesign))) {
        stop("simcolData is not available in the provided object. Please double check!")
    }
    
    
    counts <- simCounts(shinyDesign)
    n.cluster <- length(unique(simcolData(shinyDesign)$domain))
	
	seurat <- CreateSeuratObject(counts = counts, assay = 'RNA')
	seurat$cluster <- shinyDesign@simcolData$domain
	
	seurat <- NormalizeData(seurat)
	seurat <- FindVariableFeatures(seurat, selection.method = 'vst')

	all.genes <- rownames(seurat)
	seurat <- ScaleData(seurat, features = all.genes)
	seurat <- RunPCA(seurat, verbose = FALSE)

	seurat <- FindNeighbors(seurat, dims = 1:20, verbose = FALSE)
	seurat <- FindClusters(seurat, verbose = FALSE,res = 2 )


	X <- Seurat::AggregateExpression(seurat, assays=SeuratObject::DefaultAssay(seurat),  
					slot= "scale.data", group.by = "seurat_clusters")[[1]] # get average scaled expression for each variable gene and cluster

	dist1 <- dist(t(X))
	hclust1 <- hclust(dist1)

	clust2 <- cutree(hclust1, k = n.cluster)
	new_labels <- clust2
	names(new_labels) <- levels(seurat)
	seurat <- RenameIdents(seurat,new_labels)
	seurat$pred <- Idents(seurat)
	
    nmi <- NMI(seurat$cluster, seurat$pred, variant = "sqrt")
    shinyDesign@simcolData$label_pred <- seurat$pred
    shinyDesign@NMI <- nmi

    message('Performance evaluation complete. NMI: ', round(nmi,3))
    return(shinyDesign)
}
