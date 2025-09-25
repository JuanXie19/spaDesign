library(shiny)
library(ggplot2)
library(shinyDesign2)
library(DT)
library(shinyBS)
library(hdf5r)
library(scam)
library(plotly)
library(dplyr)
library(parallel)
library(pbmcapply)
library(future.apply)
library(BRISC)

options(shiny.maxRequestSize = 300 * 1024^2)

# pilot/reference data paths (from your package)
reference_data_paths <- list(
  "Chicken Heart" = system.file("extdata/ref_chicken_heart.rds", package = "shinyDesign2"),
  "Human Brain"   = system.file("extdata/ref_human_brain.rds", package = "shinyDesign2")
)


server <- function(input, output, session) {
  data_obj <- dataInputServer("data_input", reference_data_paths)
  analysisServer("analysis", data_obj)
}



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
#' @export

dataInputServer <- function(id, reference_data_paths) {
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
        data <- readRDS(reference_data_paths[[input$reference_dataset]])
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
            seurat <- Seurat::CreateSeuratObject(counts = counts, assay = 'RNA')
            seurat <- Seurat::NormalizeData(seurat)
            seurat <- Seurat::FindVariableFeatures(seurat, selection.method = 'vst')
            all.genes <- rownames(seurat)
            seurat <- Seurat::ScaleData(seurat, features = all.genes)
            seurat <- Seurat::RunPCA(seurat)
            seurat <- Seurat::RunUMAP(seurat, dims = 1:30)
            seurat <- Seurat::FindNeighbors(seurat, dims = 1:30)
            seurat <- Seurat::FindClusters(seurat, resolution = 2)
            
            X <- Seurat::AggregateExpression(seurat, assays=DefaultAssay(seurat),
                                             slot="scale.data", group.by="seurat_clusters")[[1]]
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
      ggplot(df, aes(x=x, y=y, color=as.factor(domain))) +
        geom_point(size=2) + theme_classic() + labs(title="Spatial domains", x="X", y="Y")
    })
    
    # Return the final reactive object for use in other modules
    return(final_data_obj)
  })
}


###
detect_saturation <- function(model, threshold = 0.005, consecutive = 3){
  new_data <- data.frame(real_seq_depth = seq(min(model$model$real_seq_depth), 
                                              max(model$model$real_seq_depth), 
                                              length.out = 200))
  
  new_data$NMI_pred <- predict(model, new_data)
  
  # approximate derivative
  new_data$deriv <- c(NA, diff(new_data$NMI_pred)/diff(new_data$real_seq_depth))
  
  # find first index where slope stays below threshold for 'consecutive' points
  below_thresh <- new_data$deriv < threshold
  run_len <- rle(below_thresh)
  consec_idx <- which(run_len$values & run_len$lengths >= consecutive)
  
  if(length(consec_idx) ==0){
    return(NA)
  }
  
  # map back to position in new_data
  start_pos <- sum(run_len$lengths[seq_len(consec_idx[1] - 1)]) + 1
  sat_depth <- new_data$real_seq_depth[start_pos]
  return(sat_depth)
}


#' Analysis Server Module
#'
#' Server logic for simulation and analysis tab.
#'
#' @param id Module namespace.
#' @param data_obj A reactive returning the processed data.
#' @return Reactive simulation results.
#' @export


analysisServer <- function(id, data_obj){
  moduleServer(id, function(input, output, session){
    sim_results <- reactiveVal(list())
    original_seq_depth <- reactiveVal(NULL)
    
    observeEvent(input$run_sim, {
      req(data_obj())
      original_seq_depth_val <- sum(data_obj()@refCounts)/1e6
      original_seq_depth(original_seq_depth_val)
      
      seq_range <- c(0.5, 1:3, 5, 7, 10)
      withProgress(message = "Running simulations...", {
        base_res <- powerAnalysisEffectSize(data_obj(), es_range = 1, 
                                            seq_depth_range = seq_range, n_rep = 5, n_cores = 1)
        base_res$condition <- 'Original'
        results <- list(base = base_res)
        
        if(input$simulation_type == "effect_size"){
          effect_res <- powerAnalysisEffectSize(data_obj(), es_range = input$effect_size,
                                                seq_depth_range = seq_range, n_rep = 5, n_cores = 1)
          effect_res$condition <- paste('Effect Size:', input$effect_size)
          results$effect <- effect_res
          
        } else if(input$simulation_type == "spatial"){
          spatial_res <- powerAnalysisSpatial(data_obj(), SIGMA = input$sigma, prop_range = 0.8,
                                              seq_depth_range = seq_range, n_rep = 5, n_cores = 1)
          spatial_res$condition <- paste('Spatial (Sigma):', input$sigma)
          results$spatial <- spatial_res
        }
        sim_results(results)
      })
    })
    
    combined_raw_data <- reactive({
      results <- sim_results()
      req(length(results) > 0)
      req(original_seq_depth())
      
      dplyr::bind_rows(results) %>%
        mutate(real_seq_depth = seq_depth * original_seq_depth())
    })
    
    processed_plot_data <- reactive({
      data_to_process <- combined_raw_data()
      req(data_to_process)
      
      data_to_process %>%
        group_by(real_seq_depth, condition) %>%
        summarise(mean_NMI = mean(NMI),
                  se_NMI = sd(NMI) / sqrt(n()),
                  .groups = 'drop') %>%
        mutate(mean_NMI = round(mean_NMI, 3),
               se_NMI = round(se_NMI, 3),
               real_seq_depth = round(real_seq_depth, 3))
    })
    
    saturation_points <- reactiveVal(list())
    
    output$sim_plot <- renderPlotly({
      saturation_points(list())
      plot_data_summary <- processed_plot_data()
      
      validate(
        need(nrow(plot_data_summary) > 0, "Press 'Run simulation' to see the results.")
      )
      
      p <- ggplot(plot_data_summary, aes(x = real_seq_depth, y = mean_NMI, color = condition,
                                         text = paste("Condition:", condition, "<br>", "Seq.Depth (million):",
                                                      real_seq_depth, "<br>", "Mean NMI:", mean_NMI))) +
        geom_point() + 
        geom_errorbar(aes(ymin = mean_NMI - se_NMI, ymax = mean_NMI + se_NMI), width = 0.1) + 
        labs(title = 'Mean NMI vs. Sequencing depth', x = 'Sequencing depth (million)', y = 'Mean NMI') + 
        ylim(0,1) + 
        theme_minimal()
      
      plot_data_raw <- combined_raw_data()
      current_saturation_points <- list()
      
      saturation_df <- data.frame()
      
      for (i in seq_along(unique(plot_data_raw$condition))){
        
        cond <- unique(plot_data_raw$condition)[i]
        df_subset <- plot_data_raw %>% filter(condition == cond)
        
        scam_model <- scam(NMI~s(real_seq_depth, bs = 'mpi', k = 6), data = df_subset)
        
        new_data <- data.frame(real_seq_depth = seq(min(df_subset$real_seq_depth), 
                                                    max(df_subset$real_seq_depth),
                                                    length.out = 100))
        new_data$NMI_pred <- predict(scam_model, new_data)
        new_data$condition <- cond
        
        p <- p + geom_line(data = new_data, aes(y = NMI_pred, text = NULL, color = condition), linewidth = 1)
        
        ## saturation detection
        sat_depth <- detect_saturation(scam_model, threshold = 0.005, consecutive = 3)
        current_saturation_points[[cond]] <- sat_depth
        
        if(!is.na(sat_depth)){
         saturation_df <- rbind(saturation_df,
                                data.frame(condition = cond,
                                           sat_depth = sat_depth,
                                           label = paste('Sat:', round(sat_depth, 2))))
        }
      }
      
      if(nrow(saturation_df) > 0){
        saturation_df$label_y_offset <- (seq_along(saturation_df$sat_depth) %% 2) * 0.08
        saturation_df$label_y <- 0.05 + saturation_df$label_y_offset
        
        p <- p + 
          geom_vline(data = saturation_df, aes(xintercept = sat_depth, color = condition), linetype = 'dashed') +
          # Use geom_text instead of annotate for data-driven labels
          geom_text(data = saturation_df, aes(x = sat_depth, y = label_y, label = label, color = condition), 
                    angle = 90, vjust = -0.5, hjust = 0, size = 3, inherit.aes = FALSE)
      }
      
      saturation_points(current_saturation_points)
      ggplotly(p, tooltip = 'text')
    })
    
    output$sim_table <- DT::renderDT({
      table_data <- processed_plot_data()
      
      validate(
        need(nrow(table_data) > 0, "No data available for the table. Run a simulation first.")
      )
      
      table_data <- table_data %>% 
        dplyr::select(c(real_seq_depth, mean_NMI, se_NMI, condition)) %>%
        dplyr::arrange(condition) %>%
        dplyr::rename('Sequencing depth (millions)' = real_seq_depth)
      
      DT::datatable(
        table_data,
        options = list(
          pageLength = 10,
          lengthMenu = c(5, 10, 25, 50),
          scrollX = TRUE
        ),
        rownames = FALSE,
        selection = 'none'
      )
    }, server = FALSE)
    
    output$download_table <- downloadHandler(
      filename = function() {
        paste("simulation-results-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(processed_plot_data(), file, row.names = FALSE)
      }
    )
  })
}

