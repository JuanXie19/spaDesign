#' Internal dependency declarations
#'
#' @importFrom Seurat AggregateExpression CreateSeuratObject FindClusters FindNeighbors FindVariableFeatures Idents NormalizeData Read10X_h5 RunPCA ScaleData
#' @importFrom SeuratObject DefaultAssay
#' @importFrom data.table fread
#' @importFrom ggplot2 aes element_text geom_point ggplot labs theme theme_classic
#' @importFrom graphics plot.new text
#' @importFrom shiny moduleServer observeEvent outputOptions reactive reactiveVal renderPlot renderPrint req showNotification updateSelectInput withProgress
#' @importFrom stats cutree dist hclust
#'
#' @noRd


dataInputServer <- function(id, reference_data_paths) {
  moduleServer(id, function(input, output, session) {
    
    updateSelectInput(session, "reference_dataset", choices = names(reference_data_paths))
    
    #================================================================
    # 1. Reactive for Raw Data Input
    #================================================================
    raw_data <- reactive({
      if (input$data_source == "reference") {
        req(input$reference_dataset)
        
        tryCatch({
          data <- readRDS(reference_data_paths[[input$reference_dataset]])
          return(list(counts = data@refCounts, coords = data@refcolData))
        }, error = function(e) {
          showNotification(
            paste("Error reading reference data:", e$message),
            type = "error",
            duration = 10
          )
          return(NULL)
        })
        
      } else {
        req(input$h5_file, input$positions_file)
        
        tryCatch({
          withProgress(message = 'Reading uploaded data...', {
            counts <- Seurat::Read10X_h5(filename = input$h5_file$datapath)
            
            position_names <- c("barcode", "in_tissue", "array_row", "array_col",
                                "px_row", "px_col")
            coords_raw <- data.table::fread(input$positions_file$datapath, header = FALSE,
				colClasses = "character")
            if (ncol(coords_raw) != length(position_names)) {
              stop("Positions CSV must contain exactly six columns.")
            }
            first_row <- unlist(coords_raw[1L], use.names = FALSE)
            if (identical(first_row, c("barcode", "in_tissue", "array_row",
                                      "array_col", "pxl_row_in_fullres",
                                      "pxl_col_in_fullres")) ||
                identical(first_row, position_names)) {
              coords_raw <- coords_raw[-1L, ]
            }
            names(coords_raw) <- position_names
            coords_raw$px_row <- suppressWarnings(as.numeric(coords_raw$px_row))
            coords_raw$px_col <- suppressWarnings(as.numeric(coords_raw$px_col))
            if (anyNA(coords_raw$in_tissue) ||
                !all(coords_raw$in_tissue %in% c("0", "1")) ||
                any(!is.finite(coords_raw$px_row)) ||
                any(!is.finite(coords_raw$px_col))) {
              stop("Positions CSV contains invalid tissue flags or pixel coordinates.")
            }
			
            coords <- coords_raw[coords_raw$in_tissue == 1, ]
            
            if (nrow(coords) == 0) {
              stop("No spots marked as 'in_tissue' found.")
            }
            
            coords <- data.frame(
              row.names = coords$barcode,
              x = coords$px_row,
              y = coords$px_col
            )
            
            common_barcodes <- intersect(colnames(counts), rownames(coords))
            if (length(common_barcodes) == 0) {
              stop("No common barcodes found between H5 and positions file.")
            }
            
            counts <- counts[, common_barcodes, drop = FALSE]
            coords <- coords[common_barcodes, ]
            
            showNotification(
              paste("Successfully loaded", length(common_barcodes), "spots"),
              type = "message",
              duration = 3
            )
            
            return(list(counts = counts, coords = coords))
          })
        }, error = function(e) {
          showNotification(
            paste("Error reading uploaded files:", e$message),
            type = "error",
            duration = 10
          )
          return(NULL)
        })
      }
    })
    
    anno_provided <- reactive({
      !is.null(input$anno_file)
    })
    
    #================================================================
    # 2. UI Logic for Clustering
    #================================================================
    output$need_clusters_ui <- reactive({
      input$data_source == "upload" && !is.null(raw_data()) && !anno_provided()
    })
    outputOptions(output, "need_clusters_ui", suspendWhenHidden = FALSE)
    
    #================================================================
    # 3. Reactive for Processed Coordinates (from working version)
    #================================================================
    processed_coords <- reactive({
      req(raw_data())
      data <- raw_data()
      coords <- data$coords
      
      if (input$data_source == 'upload') {
        if (anno_provided()) {
          withProgress(message = 'Reading user-provided annotations...', {
            anno_data <- fread(input$anno_file$datapath, header = TRUE)
            if (!all(c("barcode", "domain") %in% names(anno_data))) {
              stop("Annotation CSV must contain barcode and domain columns.")
            }
            anno_data <- anno_data[, c("barcode", "domain"), with = FALSE]
            if (anyDuplicated(anno_data$barcode)) {
              stop("Annotation CSV contains duplicate barcodes.")
            }
            
            common_barcodes <- intersect(rownames(coords), anno_data$barcode)
            if (length(common_barcodes) == 0) {
              stop("No common barcodes found between spatial coordinates and annotation file.")
            }
            
            coords$barcode <- rownames(coords)
            coords <- merge(coords, anno_data, by = "barcode", all.x = TRUE)
            rownames(coords) <- coords$barcode
            coords$barcode <- NULL
            coords <- coords[!is.na(coords$domain), ]
          })
        } else {
          req(input$n_clusters, input$n_clusters > 1)
          counts <- data$counts
          
          withProgress(message = "Running clustering to predict domains...", {
            seurat <- Seurat::CreateSeuratObject(counts = counts, assay = 'RNA')
            seurat <- Seurat::NormalizeData(seurat)
            seurat <- Seurat::FindVariableFeatures(seurat, selection.method = 'vst')
            all.genes <- rownames(seurat)
            seurat <- Seurat::ScaleData(seurat, features = all.genes)
            seurat <- Seurat::RunPCA(seurat)
            seurat <- Seurat::FindNeighbors(seurat, dims = 1:30)
            seurat <- Seurat::FindClusters(seurat, resolution = 2)
            
            X <- Seurat::AverageExpression(seurat, assays = DefaultAssay(seurat),
                                             layer = "scale.data", group.by = "seurat_clusters")[[1]]
            clust2 <- cutree(hclust(dist(t(X))), k = input$n_clusters)
            new_labels <- clust2
            names(new_labels) <- levels(seurat)
            
            clusters <- Idents(seurat)
            coords$domain <- new_labels[as.character(clusters[colnames(counts)])]
          })
        }
      }
      return(coords)
    })
    
    #================================================================
    # 4. Event Reactive (from working version)
    #================================================================
    upload_files <- reactive({
      list(input$h5_file, input$positions_file, input$anno_file)
    })
    prepared_upload <- reactiveVal(NULL)

    observeEvent(upload_files(), {
      prepared_upload(NULL)
    }, ignoreInit = TRUE, priority = 100)

    observeEvent(input$prepare_data_button, {
      req(input$data_source == "upload")
      prepared_upload(NULL)
      
      logfc_cutoff <- if (isTRUE(input$show_advanced)) input$logfc_cutoff else 0.7
      mean_in_cutoff <- if (isTRUE(input$show_advanced)) input$mean_in_cutoff else 1.8
      max_num_gene <- if (isTRUE(input$show_advanced)) input$max_num_gene else 10
      n_cores <- if(isTRUE(input$show_advanced)) input$n_cores else 1
   
      if (.Platform$OS.type == "windows") n_cores <- 1L
   
      tryCatch({
        counts <- raw_data()$counts
        coords <- processed_coords()
        req(counts, coords)
        counts <- counts[, rownames(coords), drop = FALSE]
		  
        withProgress(message = "Creating design object...", value = 0.2, {
          DATA <- spaDesign::createDesignObject(count_matrix = counts, loc = coords)
          message("✓ Design object created")
        })
        
        withProgress(message = "Selecting domain informative genes...", value = 0.4, {
          DATA <- spaDesign::featureSelection(
            DATA,
            logfc_cutoff = logfc_cutoff,
            mean_in_cutoff = mean_in_cutoff,
            max_num_gene = max_num_gene
          )
          message("✓ Feature selection completed")
        })
        
        withProgress(message = "Learning spatial expression pattern...", value = 0.4, {
          DATA <- spaDesign::estimation_NNGP(
            DATA,
            n_neighbors =  10,
            order = 'AMMD',
            X = NULL,
            verbose = FALSE
          )
          message("✓ NNGP estimation completed")
        })
        
        
        withProgress(message = "Learning domain spatial patterns...", value = 0.8, {
          message("Starting estimation_FGEM...")
          DATA <- spaDesign::estimation_FGEM(
            DATA,
            iter_max = 500,
            M_candidates = 2:5,
            tol = 1e-1,
            n_cores = n_cores
          )
          message("✓ FGEM estimation completed")
        })
        
        showNotification(
          "Data processing complete!",
          type = "message",
          duration = 5
        )
        
        prepared_upload(DATA)
        
      }, error = function(e) {
        message("ERROR: ", e$message)
        showNotification(
          paste("Error in data processing:", e$message),
          type = "error",
          duration = 10
        )
        return(NULL)
      })
    })
    
    #================================================================
    # 5. Final Reactive
    #================================================================
    final_data_obj <- reactive({
      if (input$data_source == "reference") {
        req(input$reference_dataset)
        
        tryCatch({
          data <- readRDS(reference_data_paths[[input$reference_dataset]])
          return(data)
        }, error = function(e) {
          showNotification(
            paste("Error loading reference data:", e$message),
            type = "error",
            duration = 10
          )
          return(NULL)
        })
        
      } else {
        return(prepared_upload())
      }
    })
    
    #================================================================
    # 6. Outputs
    #================================================================
    output$summary_ref <- renderPrint({
      obj <- final_data_obj()
      req(obj)
      
      tryCatch({
        expr <- obj@refCounts
        coords <- obj@refcolData
        ngenes <- nrow(expr)
        nspots <- ncol(expr)
        ndomains <- length(unique(coords$domain))
        domain_counts <- table(coords$domain)
        
        cat("=== Data Summary ===\n\n")
        cat("Genes:  ", ngenes, "\n")
        cat("Spots:  ", nspots, "\n")
        cat("Domains:", ndomains, "\n\n")
        cat("Number of spots per domain:\n")
        print(domain_counts)
      }, error = function(e) {
        cat("Error displaying summary:", e$message, "\n")
      })
    })
    
    output$domain_plot <- renderPlot({
      obj <- final_data_obj()
      req(obj)
      
      tryCatch({
        df <- obj@refcolData
        
        ggplot(df, aes(x = x, y = y, color = as.factor(domain))) +
          geom_point(size = 2) +
          theme_classic() +
          labs(
            title = "Spatial Domains",
            x = "X Coordinate",
            y = "Y Coordinate",
            color = "Domain"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right"
          )
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error creating plot:", e$message), col = "red")
      })
    })
    
    return(final_data_obj)
  })
}