library(shiny)
library(ggplot2)
library(Matrix)
library(data.table)

dataInputUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      
   
      radioButtons(ns("data_source"), "Choose data source:",
                   choices = c("Use reference Data" = "reference",
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
        
        # conditional panel for n_clusters.It shows up only if data source is 'upload' AND annotation file is not provided
        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload' && !input['%s']", ns("data_source"), ns("anno_file")),
          numericInput(ns("n_clusters"), "Number of expected spatial domains:", value = 7, min=2, step=1)
        ),
        
        checkboxInput(ns("show_advanced"), "Show advanced feature selection options", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s']", ns("show_advanced")),
          sliderInput(ns("logfc_cutoff"), "logFC cutoff:", min=0.1, max=2, value=0.7, step=0.1),
          sliderInput(ns("mean_in_cutoff"), "Mean in cutoff:", min=0.5, max=5, value=1.8, step=0.1),
          numericInput(ns("max_num_gene"), "Max number of genes:", value=10, min=1, step=1)
        ),
        
        # action button to triggole the analysis
        hr(),
        actionButton(ns('prepare_data_button'), "Run Data Preparation", class =  "btn-primary")
      )
    ),
    mainPanel(
      plotOutput(ns("domain_plot")),
      verbatimTextOutput(ns("summary_ref"))
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
        return(list(counts = data@refCounts, coords = data@refcolData))
      } else {
        req(input$h5_file, input$positions_file)
        # read Space ranger files
        withProgress(message = 'Reading uploaded data...', {
          # Read h5 file
          library(Seurat)
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
          counts <- data$counts
        
          withProgress(message = "Running clustering to predict domains...", {
          # --- Seurat Clustering Pipeline ---
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
	    
	    nullfile <- tempfile()
	    sink(nullfile)
	    
	    tryCatch({
	      withProgress(message = "Creating design object...", {
	        DATA <- shinyDesign2::createDesignObject(count_matrix = counts, loc = coords)
	      })
	      withProgress(message = "Selecting domain informative genes...", {
	        DATA <- shinyDesign2::featureSelection(DATA,
	                                               logfc_cutoff = logfc_cutoff,
	                                               mean_in_cutoff = mean_in_cutoff,
	                                               max_num_gene = max_num_gene)
	      })
	      withProgress(message = "Learning gene spatial expression patterns...", {
	        DATA <- estimation_NNGP(DATA, n_neighbors = 10, ORDER = 'AMMD')
	      })
	      withProgress(message = "Learning domain spatial patterns...", {
	        DATA <- estimation_FGEM(DATA, iter_max = 1000, M_candidates = 2:7, tol = 1e-1)
	      })
	      
	      return(DATA)
	    }, finally = {
	      sink(type = 'message')
	      sink()
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
        return(reference_data[[input$reference_dataset]])
        
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
      df$domain <- as.factor(df$domain)
      ggplot(df, aes(x=x, y=y, color=domain)) +
        geom_point(size=2) + theme_classic() + labs(title="Spatial domains", x="X", y="Y")
    })
    
    # Return the final reactive object for use in other modules
    return(final_data_obj)
  })
}
