#' Select Domain-Informative Genes for Each Spatial Domain
#' 
#' Identifies and selects genes that are highly expressed and differential expressed
#' within each spatial domain compared to outside the domain. Selected genes are stored
#' in the \code{topGenes} slot of the \code{spaDesign} object for downstream analysis.
#'
#' @param spaDesign A \code{spaDesign} object containing:
#'   \itemize{
#'     \item Count matrix (accessible via \code{refCounts})
#'     \item Spatial coordinates and domain assignments (accessible via \code{refcolData})
#'   }
#' @param logfc_cutoff Numeric scalar > 0. Minimum absolute log fold change threshold.
#'   Genes must have |logFC| >= this value to be considered domain-informative.
#' @param mean_in_cutoff Numeric scalar > 0. Minimum mean log-transformed expression
#'   required within the domain. Filters out lowly expressed genes.
#' @param max_num_gene Integer. Maximum number of genes to select per spatial domain.
#'   If more genes pass the filters, the top \code{max_num_gene} genes with highest
#'   absolute logFC are retained.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing across domains.
#'
#'#' @return A \code{spaDesign} object with the \code{topGenes} slot updated.
#'   The \code{topGenes} slot contains a named list where each element corresponds
#'   to a spatial domain and contains a data frame of selected genes with their
#'   fold change statistics.

#' @import dplyr
#' @importFrom pbmcapply pbmclapply
#' @export
#'
#' @examples
#' \dontrun{
#' # Select domain-informative genes
#' my_spadesign <- featureSelection(
#'   spaDesign = my_spadesign,
#'   logfc_cutoff = 0.7,      # Require at least 0.7 log FC
#'   mean_in_cutoff = 2,       # Require mean log-expression >= 2 in domain
#'   max_num_gene = 50,        # Select maximum 50 genes per domain
#'   n_cores = 4
#' )
#'
#' # Access the selected genes for each domain
#' top_genes <- my_spadesign@topGenes
#' 
#' # View genes for first domain
#' head(top_genes[[1]])
#' }

featureSelection <- function(spaDesign, logfc_cutoff, mean_in_cutoff, max_num_gene, n_cores){
	
    if (!is.numeric(logfc_cutoff) || logfc_cutoff <= 0) stop("'logfc_cutoff' must be a positive numeric value.")
    if (!is.numeric(mean_in_cutoff) || mean_in_cutoff <= 0) stop("'mean_in_cutoff' must be a positive numeric value.")


	count_matrix <- refCounts(spaDesign)
	loc_file <- refcolData(spaDesign)[, c('x','y','domain')]
	
	FC_list <- geneSummary(count_matrix, loc_file, n_cores)
	
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


#' Calculate Gene-Level Fold Change Statistics Across Spatial Domains
#'
#' Computes log fold change and mean expression statistics for each gene in each spatial domain.
#' Genes are compared between inside-domain and outside-domain expression levels to identify
#' domain-specific markers.
#'
#' @param count_matrix Numeric matrix of gene expression counts (genes × spots).
#'   Rows are genes, columns are spatial spots.
#' @param loc Data frame containing spatial coordinates and domain assignments with columns:
#'   \describe{
#'     \item{x}{Numeric x-coordinate}
#'     \item{y}{Numeric y-coordinate}
#'     \item{domain}{Character or factor indicating domain membership for each spot}
#'   }
#' @param n_cores Integer. Number of CPU cores to use for parallel processing across domains.
#'   Default is 4.
#'
#' @details
#' For each domain and each gene, the function calculates:
#' \itemize{
#'   \item \strong{logFC}: Log fold change comparing mean inside-domain to mean outside-domain expression
#'   \item \strong{logFC_low}: Log fold change comparing mean inside-domain to mean of lower 50% outside-domain expression
#'   \item \strong{mean_in}: Mean log-transformed expression inside the domain
#'   \item \strong{mean_out}: Mean log-transformed expression outside the domain
#'   \item \strong{mean_out_low}: Mean log-transformed expression in the lower 50% of outside-domain spots
#' }
#'
#' Expression values are log-transformed as log(count + 1) to stabilize variance.
#' The \code{logFC_low} metric is more robust to outliers and better identifies genes
#' that are specifically enriched in a domain rather than just avoiding low-expression regions.
#'
#' @return A named list with one element per domain. Each element is a data frame with
#'   rows corresponding to genes and columns:
#'   \describe{
#'     \item{logFC}{Log fold change (inside vs outside)}
#'     \item{logFC_low}{Log fold change (inside vs lower 50% outside)}
#'     \item{mean_in}{Mean log-expression inside domain}
#'     \item{mean_out}{Mean log-expression outside domain}
#'     \item{mean_out_low}{Mean log-expression in lower 50% outside domain}
#'   }
#'   List names correspond to domain identifiers.
#'
#' @importFrom pbmcapply pbmclapply
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Assuming you have a count matrix and location data
#' fc_stats <- geneSummary(
#'   count_matrix = my_counts,
#'   loc = my_locations,
#'   n_cores = 4
#' )
#'
#' # View statistics for first domain
#' head(fc_stats[[1]])
#'} 
geneSummary <- function(count_matrix, loc, n_cores){
	
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
    }, mc.cores = n_cores)
    names(fc_results) <- domains
    message("Completed fold change calculations for all domains")
    return(fc_results)
}

