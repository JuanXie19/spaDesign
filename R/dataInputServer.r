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
            
            coords_raw <- data.table::fread(input$positions_file$datapath, header = FALSE)
            
            if (ncol(coords_raw) < 6) {
              stop("Positions file must have at least 6 columns.")
            }
            
            colnames(coords_raw) <- c("barcode", "in_tissue", "array_row", "array_col", "px_row", "px_col")
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
            
            counts <- counts[, common_barcodes]
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
            req(colnames(anno_data) %in% c('barcode', 'domain'))
            
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
            seurat <- Seurat::RunUMAP(seurat, dims = 1:30)
            seurat <- Seurat::FindNeighbors(seurat, dims = 1:30)
            seurat <- Seurat::FindClusters(seurat, resolution = 2)
            
            X <- Seurat::AggregateExpression(seurat, assays = DefaultAssay(seurat),
                                             slot = "scale.data", group.by = "seurat_clusters")[[1]]
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
    processed_design_object <- eventReactive(input$prepare_data_button, {
      req(input$data_source == 'upload')
      counts <- raw_data()$counts
      coords <- processed_coords()
      req(counts, coords)
      
      logfc_cutoff <- if (isTRUE(input$show_advanced)) input$logfc_cutoff else 0.7
      mean_in_cutoff <- if (isTRUE(input$show_advanced)) input$mean_in_cutoff else 1.8
      max_num_gene <- if (isTRUE(input$show_advanced)) input$max_num_gene else 10
      
      tryCatch({
        withProgress(message = "Creating design object...", value = 0.2, {
          DATA <- spaDesign::createDesignObject(count_matrix = counts, loc = coords)
          message("✓ Design object created")
        })
        
        withProgress(message = "Selecting domain informative genes...", value = 0.4, {
          DATA <- spaDesign::featureSelection(
            DATA,
            logfc_cutoff = logfc_cutoff,
            mean_in_cutoff = mean_in_cutoff,
            max_num_gene = max_num_gene,
            n_cores = 1
          )
          message("✓ Feature selection completed")
        })
        
        withProgress(message = "Learning domain spatial patterns...", value = 0.8, {
          message("Starting estimation_FGEM...")
          DATA <- spaDesign::estimation_FGEM(
            DATA,
            iter_max = 1000,
            M_candidates = 2:7,
            tol = 1e-1,
            n_cores = 1
          )
          message("✓ FGEM estimation completed")
        })
        
        showNotification(
          "Data processing complete!",
          type = "message",
          duration = 5
        )
        
        return(DATA)
        
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
        return(processed_design_object())
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