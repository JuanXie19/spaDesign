#' Data Input Server Module
#'
#' Server logic for the data input tab. Handles file reading, clustering, and processing.
#'
#' @param id A character string specifying the Shiny module namespace.
#' @param reference_data A list of reference datasets (optional).
#'
#' @return A reactive object containing the final processed data.
#' @import shiny
#' @import ggplot2
#' @import data.table
#' @import Seurat
#' @noRd

dataInputServer <- function(id, reference_data_paths) {
  moduleServer(id, function(input, output, session) {
    
    # update choices for reference data once
    updateSelectInput(session, "reference_dataset", choices = names(reference_data_paths))
    
    #================================================================
    # 1. Reactive for Raw Data Input
    #    - This reactive's job is to return a list containing the raw
    #      counts and coordinates, whether from reference or upload.
    #    - It only re-runs if the data source or selected files change.
    #================================================================
    raw_data <- reactive({
      if (input$data_source == "reference") {
        req(input$reference_dataset)
        # read the file to extract counts/coords for preview
        data <- readRDS(reference_data_paths[[input$reference_dataset]])
        
        # Assuming reference data is a list with 'counts' and 'coords'
        return(list(counts = data@refCounts, coords = data@refcolData))
      } else {
        req(input$h5_file, input$positions_file)
        
        # read Space ranger files
        withProgress(message = 'Reading uploaded data...', {
          # Read h5 file
          counts <- Seurat::Read10X_h5(filename = input$h5_file$datapath)
          
          # Read tissue_positions_list.csv
          # Assuming standard SpaceRanger tissue_positions_list.csv format:
          # barcode, in_tissue, array_row, array_col, px_row_in_full_image, px_col_in_full_image
          coords_raw <- data.table::fread(input$positions_file$datapath, header = FALSE)
          colnames(coords_raw) <- c("barcode", "in_tissue", "array_row", "array_col", "px_row", "px_col")
          
          # Filter for spots "in_tissue" and select relevant columns for coords
          coords <- coords_raw[coords_raw$in_tissue == 1, ]
          coords <- data.frame(
            row.names = coords$barcode,
            x = coords$px_row, # Use pixel column as x
            y = coords$px_col  # Use pixel row as y
          )
          
          # Ensure barcodes match between counts and coords
          common_barcodes <- intersect(colnames(counts), rownames(coords))
          if (length(common_barcodes) == 0) {
            stop("No common barcodes found between expression data (H5) and spatial coordinates. Please check your uploaded files.")
          }
          counts <- counts[, common_barcodes]
          coords <- coords[common_barcodes, ]
          
          return(list(counts = counts, coords = coords))
          
        })
      }
    })
    
    
    # reactive to check if annotation file is provided
    anno_provided <- reactive({
      # reactive on the fileInput value, which is a list
      # if the file is not selected, the list will be NULL
      !is.null(input$anno_file)
    })
    
    #================================================================
    # 2. UI Logic for Clustering
    #    - This output controls the visibility of the n_clusters input.
    #    - It's based on the `anno_provided` flag
    #================================================================
    output$need_clusters_ui <- reactive({
      input$data_source == "upload" && !is.null(raw_data()) && !anno_provided()
    })
    outputOptions(output, "need_clusters_ui", suspendWhenHidden = FALSE)
    
    
    #================================================================
    # 3. Reactive for Processed Coordinates
    #    - This reactive has two branches
    #     1. if 'anno_provided()' is TRUE, read the file and join
    #     2. if 'anno_provided()' is FALSE, run Seurat clustering
    #================================================================
    # 3a. Create Base Seurat Object (Only runs once on upload)
    base_seurat_object <- reactive({
      req(input$data_source == 'upload', !anno_provided())
      data <- raw_data()
      req(data)
      
      withProgress(message = "Preprocessing expression data...", {
        s_obj <- Seurat::CreateSeuratObject(counts = data$counts, assay = 'RNA')
        s_obj <- Seurat::NormalizeData(s_obj)
        s_obj <- Seurat::FindVariableFeatures(s_obj, selection.method = 'vst')
        s_obj <- Seurat::ScaleData(s_obj, features = rownames(s_obj))
        s_obj <- Seurat::RunPCA(s_obj, verbose = FALSE)
        # Find neighbors and high-res clusters once
        s_obj <- Seurat::FindNeighbors(s_obj, dims = 1:20)
        s_obj <- Seurat::FindClusters(s_obj, resolution = 2, verbose = FALSE)
        return(s_obj)
      })
    })
    
    
    
    
    processed_coords <- reactive({
      req(raw_data())
      data <- raw_data()
      coords <- data$coords
      
      if (input$data_source =='upload') {
        if(anno_provided()) {
          # branch 1: use provided annotations
          withProgress(message = 'Reading user-provided annotations...', {
            anno_data <- fread(input$anno_file$datapath, header = TRUE)
            req(colnames(anno_data) %in% c('barcode', 'domain'))
            
            # ensure barcodes in anno_data match thosee in coords
            common_barcodes <- intersect(rownames(coords), anno_data$barcode)
            if (length(common_barcodes) == 0) {
              stop("No common barcodes found between spatial coordinates and annotation file. Please check your files.")
            }
            
            # Merge annotations with coordinates
            coords$barcode <- rownames(coords)
            coords <- merge(coords, anno_data, by = "barcode", all.x = TRUE)
            rownames(coords) <- coords$barcode
            coords$barcode <- NULL
            
            # Remove spots without an annotation
            coords <- coords[!is.na(coords$domain), ]
          })
        } else{
          # branch 2: no annotation provided, so run clustering
          req(input$n_clusters, input$n_clusters > 1)
          
          seurat_obj <- base_seurat_object()

          withProgress(message = "Running clustering to predict domains...", {
            
            X <- Seurat::AggregateExpression(seurat_obj, assays=DefaultAssay(seurat),
                                             slot="scale.data", group.by="seurat_clusters")[[1]]
            clust2 <- cutree(hclust(dist(t(X))), k = input$n_clusters)
            # Map back to cells
            current_clusters <- Idents(seurat_obj)
            # Create a named vector for mapping
            cluster_map <- clust2[as.character(current_clusters)]
            names(cluster_map) <- names(current_clusters)
            
            coords$domain <- cluster_map[rownames(coords)]
          })
        }
      }
      return(coords)
    })
    
    #=====================
    #. event Reactive that runs only when button click
    #=====================
    processed_design_object <- eventReactive(input$prepare_data_button, {
      req(input$data_source == 'upload')
      counts <- raw_data()$counts
      coords <- processed_coords()
      req(counts, coords)
      
      logfc_cutoff <- if (isTRUE(input$show_advanced)) input$logfc_cutoff else 0.7
      mean_in_cutoff <- if (isTRUE(input$show_advanced)) input$mean_in_cutoff else 1.8
      max_num_gene <- if (isTRUE(input$show_advanced)) input$max_num_gene else 10
      
      
      tryCatch({
        withProgress(message = "Creating design object...", {
          DATA <- createDesignObject(count_matrix = counts, loc = coords)
        })
        withProgress(message = "Selecting domain informative genes...", {
          DATA <- featureSelection(DATA,
                                   logfc_cutoff = logfc_cutoff,
                                   mean_in_cutoff = mean_in_cutoff,
                                   max_num_gene = max_num_gene, n_cores = 1)
        })
        #withProgress(message = "Learning gene spatial expression patterns...", {
        # DATA <- estimation_NNGP(DATA, n_neighbors = 10, order = 'AMMD')
        #})
        withProgress(message = "Learning domain spatial patterns...", {
          DATA <- estimation_FGEM(DATA, iter_max = 1000, M_candidates = 2:7, tol = 1e-1, n_cores = 1)
        })
        
        return(DATA)
      })
    })
    
    
    #================================================================
    # 4. Final Reactive: The Main Logic Branch
    #    - If reference data is selected, it's returned immediately.
    #    - If data is uploaded, this triggers the full processing pipeline.
    #================================================================
    final_data_obj <- reactive({
      if (input$data_source == "reference") {
        # --- PATH 1: REFERENCE DATA ---
        # Data is pre-processed, just return it.
        req(input$reference_dataset)
        return(reference_data_paths[[input$reference_dataset]])
        
      } else {
        return(processed_design_object())
      }
    })
    
    #================================================================
    # 5. Outputs
    #    - These render functions now depend on the final reactive object,
    #      which works for both reference and uploaded data paths.
    #================================================================
    output$summary_ref <- renderPrint({
      obj <- final_data_obj()
      req(obj)
      
      expr <- obj@refCounts
      coords <- obj@refcolData
      ngenes <- nrow(expr)
      nspots <- ncol(expr)
      ndomains <- length(unique(coords$domain))
      domain_counts <- table(coords$domain)
      
      cat("Genes:  ", ngenes, "\n")
      cat("Spots:  ", nspots, "\n")
      cat("Domains:", ndomains, "\n\n")
      cat("Number of spots per domain:\n")
      print(domain_counts)
    })
    
    output$domain_plot <- renderPlot({
      obj <- final_data_obj()
      req(obj)
      
      df <- obj@refcolData
      ggplot(df, aes(x=x, y=y, color=as.factor(domain))) +
        geom_point(size=2) + theme_classic() + labs(title="Spatial domains", x="X", y="Y")
    })
    
    # Return the final reactive object for use in other modules
    return(final_data_obj)
  })
}


