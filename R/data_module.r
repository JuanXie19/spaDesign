library(shiny)
library(ggplot2)
library(Matrix)
library(data.table)

dataInputUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns("data_source"), "Choose data source:",
                   choices = c("Use reference data" = "reference",
                               "Upload your own data" = "upload"),
                   selected = "reference"),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'reference'", ns("data_source")),
        selectInput(ns("reference_dataset"), "Select reference dataset:", choices = NULL)
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
        h4('Upload SpaceRanger Output Files'),
        fileInput(ns("h5_file"), "Choose filtered_feature_bc_matrix.h5 file", accept = ".h5"),
        fileInput(ns("positions_file"), "Upload tissue_positions_list.csv file", accept = '.csv'),
        fileInput(ns("anno_file"), "Upload annotations (.csv, optional):"),
        
        checkboxInput(ns("show_advanced"), "Show advanced feature selection options", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s']", ns("show_advanced")),
          sliderInput(ns("logfc_cutoff"), "logFC cutoff:", min=0.1, max=2, value=0.7, step=0.1),
          sliderInput(ns("mean_in_cutoff"), "Mean in cutoff:", min=0.5, max=5, value=1.8, step=0.1),
          numericInput(ns("max_num_gene"), "Max number of genes:", value=10, min=1, step=1)
        )
      ),
      # NEW: number of clusters if clustering needed
      conditionalPanel(
        condition = sprintf("output['%s']", ns("need_clusters_ui")),
        numericInput(ns("n_clusters"), "Number of expected clusters:", value = 7, min=2, step=1)
      )
    ),
    mainPanel(
      verbatimTextOutput(ns("summary_ref")),
      plotOutput(ns("domain_plot"))
    )
  )
}

dataInputServer <- function(id, reference_data) {
  moduleServer(id, function(input, output, session) {
  
	# update choices for reference data once
    updateSelectInput(session, "reference_dataset", choices = names(reference_data))
    
	 #================================================================
    # 1. Reactive for Raw Data Input
    #    - This reactive's job is to return a list containing the raw
    #      counts and coordinates, whether from reference or upload.
    #    - It only re-runs if the data source or selected files change.
    #================================================================
	raw_data <- reactive({
		if (input$data_source == "reference") {
        req(input$reference_dataset)
        data <- reference_data[[input$reference_dataset]]
        # Assuming reference data is a list with 'counts' and 'coords'
        return(list(counts = data@refCounts, coords = data@refcolData, needs_clustering = FALSE))
      } else {
        
        req(input$h5_file, input$positions_file)
        # read Space ranger files
        withProgress(message = 'Reading uploaded data...', {
          # Read matrix.mtx
          counts <- Read10X_h5(filename = input$h5_file$datapath)
         
          # Read tissue_positions_list.csv
          # Assuming standard SpaceRanger tissue_positions_list.csv format:
          # barcode, in_tissue, array_row, array_col, px_row_in_full_image, px_col_in_full_image
          coords_raw <- fread(input$positions_file$datapath, header = FALSE)
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
          
          # Check if clustering is needed (i.e., 'domain' column is missing)
          needs_clustering <- !"domain" %in% colnames(coords)
          
          return(list(counts = counts, coords = coords, needs_clustering = needs_clustering))
          
        })
      }
	})
	
	#================================================================
    # 2. UI Logic for Clustering
    #    - This output controls the visibility of the n_clusters input.
    #    - It's based on the `needs_clustering` flag from the raw_data() reactive.
    #================================================================
    output$need_clusters_ui <- reactive({
      raw_data()$needs_clustering
    })
    outputOptions(output, "need_clusters_ui", suspendWhenHidden = FALSE)
	
	#================================================================
    # 3. Reactive for Processed Coordinates
    #    - This reactive takes the raw coordinates and, if necessary,
    #      runs the expensive Seurat clustering to add the 'domain' column.
    #    - This step is now isolated and only runs when the raw data changes
    #      or the number of clusters is adjusted.
    #================================================================
    processed_coords <- reactive({
      data <- raw_data()
      coords <- data$coords
      
      if (data$needs_clustering) {
        req(input$n_clusters, input$n_clusters > 1)
        counts <- data$counts
        
        withProgress(message = "Running clustering to predict domains...", {
          # --- Seurat Clustering Pipeline ---
          # (Your original Seurat code is perfect here)
          library(Seurat)
          library(SeuratObject)
          seurat <- CreateSeuratObject(counts = counts, assay = 'RNA')
          seurat <- NormalizeData(seurat)
          seurat <- FindVariableFeatures(seurat, selection.method = 'vst')
          all.genes <- rownames(seurat)
          seurat <- ScaleData(seurat, features = all.genes)
          seurat <- RunPCA(seurat)
          seurat <- RunUMAP(seurat, dims = 1:30)
          seurat <- FindNeighbors(seurat, dims = 1:30)
          seurat <- FindClusters(seurat, resolution = 2)
          
          X <- Seurat::AggregateExpression(seurat, assays=DefaultAssay(seurat),
                                           slot="scale.data", group.by="seurat_clusters")[[1]]
          dist1 <- dist(t(X))
          hclust1 <- hclust(dist1)
          clust2 <- cutree(hclust1, k = input$n_clusters)
          new_labels <- clust2
          names(new_labels) <- levels(seurat)
          
          clusters <- Idents(seurat)
          coords$domain <- new_labels[as.character(clusters[colnames(counts)])]
        })
      }
      return(coords)
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
        return(reference_data[[input$reference_dataset]])
        
      } else {
        # --- PATH 2: UPLOADED DATA ---
        # Run the full pipeline on the new data.
        counts <- raw_data()$counts
        coords <- processed_coords()
        req(counts, coords)
        
        # Get advanced parameters
        logfc_cutoff <- if (isTRUE(input$show_advanced)) input$logfc_cutoff else 0.7
        mean_in_cutoff <- if (isTRUE(input$show_advanced)) input$mean_in_cutoff else 1.8
        max_num_gene <- if (isTRUE(input$show_advanced)) input$max_num_gene else 10
        
        # create a temporary file to sink the output
        nullfile <- tempfile()
        sink(nullfile)
        
        tryCatch({
        # Now build shinyDesign pipeline
        withProgress(message = "Creating design object...", {
          DATA <- shinyDesign2::createDesignObject(count_matrix = counts, loc = coords)
        })
        withProgress(message = "Running feature selection...", {
          DATA <- shinyDesign2::featureSelection(DATA,
                                                 logfc_cutoff = logfc_cutoff,
                                                 mean_in_cutoff = mean_in_cutoff,
                                                 max_num_gene = max_num_gene)
        })
        withProgress(message = "Fitting BRISC model...", {
          DATA <- estimation_NNGP(DATA, n_neighbors = 10, ORDER = 'AMMD')
        })
        withProgress(message = "Fitting FG model...", {
          DATA <- estimation_FGEM(DATA, iter_max = 1000, M_candidates = 2:7, tol = 1e-1)
        })
        
        return(DATA)
        }, finally = {
        sink(type = 'message')
        sink()
        })
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
