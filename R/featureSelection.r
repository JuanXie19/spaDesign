#' Select Domain-Informative Genes for Each Spatial Domain
#' This function selects genes with large absolute fold change and high within-domain expression
#' for each spatial domain. The selected genes are stored in the `topGenes` slot of the
#' `spaDesign` object.
#'
#' @param spaDesign A \code{spaDesign} object.
#' @param logfc_cutoff Numeric value specifying the threshold for log fold change, 
#' used to select genes with significant differential expression within vs outside domain.
#' @param mean_in_cutoff Numeric value specifying the minimum mean log-transformed expression within the domain.
#' @param max_num_gene Integer specifying the maximum number of genes that can be selected per spatial domain.
#' @return Returns the input \code{spaDesign} object with the `topGenes` slot updated
#'         with the selected genes for each domain.
#'         
#' @details Genes are selected by filtering based on `logfc_cutoff` and `mean_in_cutoff`.
#'          If too few genes meet both criteria, the union of individual criteria is used.
#'          The top genes are then ordered by descending absolute logFC and truncated
#'          to `max_num_gene`.

#' @import dplyr
#' @import pbmcapply
#' @export
#' @examples
#' # Assuming you have a spaDesign object named `my_spadesign`
#' my_spadesign <- featureSelection(my_spadesign, logfc_cutoff = 0.7, mean_in_cutoff = 2, max_num_gene = 10)
#' # Access the selected genes
#' top_genes <- my_spadesign@topGenes

featureSelection <- function(spaDesign, logfc_cutoff, mean_in_cutoff, max_num_gene, n_cores){
	
    if (!is.numeric(logfc_cutoff) || logfc_cutoff <= 0) stop("'logfc_cutoff' must be a positive numeric value.")
    if (!is.numeric(mean_in_cutoff) || mean_in_cutoff <= 0) stop("'mean_in_cutoff' must be a positive numeric value.")


	count_matrix <- refCounts(spaDesign)
	loc_file <- refcolData(spaDesign)[, c('x','y','domain')]
	
	FC_list <- geneSummary(count_matrix, loc_file)
	
	top_genes <- pbmcapply::pbmclapply(FC_list, function(DF) {
        message("Selecting genes with large absolute fold change and large within-domain expression")
        idx <- which(DF$mean_in >= mean_in_cutoff & abs(DF$logFC_low) >= logfc_cutoff)

        if (length(idx) == 0) {
            idx1 <- which(DF$mean_in >= mean_in_cutoff)
            idx2 <- which(DF$logFC_low >= logfc_cutoff)
            idx <- base::union(idx1, idx2)
        }
        
        selected <- DF[idx, ,drop = FALSE]
        if (nrow(selected) > max_num_gene) {
            selected_genes <- DF[idx, ] %>%
                dplyr::arrange(-abs(logFC_low)) %>%
                head(max_num_gene)
        } else {
            selected_genes <- DF[idx, ] %>%
                dplyr::arrange(-abs(logFC_low))
        }
        message("Completed gene selection")
        return(selected_genes)
    }, mc.cores = n_cores)
    names(top_genes) <- names(FC_list)
	  spaDesign@topGenes <- top_genes
	  message("Completed gene selection for all domains.")
	  
    return(spaDesign)
}


#' calculate the fold change for each gene across spatial domains
#'
#' @param count_matrix A gene expression count matrix. 
#' @param loc A dataframe containing the spatial coordiantes as well as domain information for each spot
#' @import pbapply
#' @return A list with fold change for each domain per gene
#' @export
#'
geneSummary <- function(count_matrix, loc){
	
	count_matrix <- as.matrix(count_matrix)
	domains <- sort(unique(loc$domain))
	log_count <- log(count_matrix + 1)
	
	fc_results <- pbmclapply(domains, function(domain) {
        message("Calculating fold change for domain: ", domain)
        rst <- sapply(seq_len(nrow(log_count)), function(gene) {
            spot_idx <- which(loc$domain == domain)
            count_inside <- log_count[gene, spot_idx]
            count_outside <- log_count[gene, -spot_idx]
            out_lower <- count_outside[count_outside <= median(count_outside)]

            mean_in <- mean(count_inside)
            mean_out <- mean(count_outside)
            mean_out_lower <- mean(out_lower)

            fc <- mean_in - mean_out
            fc_lower <- mean_in - mean_out_lower

            c(logFC = fc, logFC_low = fc_lower, mean_in = mean_in, mean_out = mean_out, mean_out_lower = mean_out_lower)
        })
        rst <- as.data.frame(t(rst))
        rownames(rst) <- rownames(log_count)
        colnames(rst) <- c("logFC", "logFC_low", "mean_in", "mean_out", "mean_out_low")
        message("Completed fold change calculation for domain: ", domain)
        rst
    }, mc.cores = 4)
    names(fc_results) <- domains
    message("Completed fold change calculations for all domains")
    return(fc_results)
}

