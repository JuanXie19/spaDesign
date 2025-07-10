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
        fileInput(ns("expr_file"), "Upload expression matrix (.rds):"),
        fileInput(ns("coord_file"), "Upload coordinates (.rds):"),
        fileInput(ns("anno_file"), "Upload annotations (.rds, optional):"),
        
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
        numericInput(ns("n_clusters"), "Number of expected clusters:", value = 5, min=2, step=1)
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
    updateSelectInput(session, "reference_dataset", choices = names(reference_data))
    
    data_obj <- reactiveVal(NULL)
    need_clusters <- reactiveVal(FALSE)
    coords_cache <- reactiveVal(NULL)
    counts_cache <- reactiveVal(NULL)
    
    observe({
      if (input$data_source == "reference") {
        req(input$reference_dataset)
        data_obj(reference_data[[input$reference_dataset]])
        need_clusters(FALSE)
      }
    })
    
    observeEvent({
      input$expr_file
      input$coord_file
    }, {
      req(input$expr_file, input$coord_file)
      counts <- readRDS(input$expr_file$datapath)
      coords <- readRDS(input$coord_file$datapath)
      counts_cache(counts)
      coords_cache(coords)
      
      # Check for 'domain'
      if (!"domain" %in% colnames(coords)) {
        need_clusters(TRUE)
      } else {
        need_clusters(FALSE)
      }
    })
    
    # provide UI output to show n_clusters input only when needed
    output$need_clusters_ui <- reactive({
      need_clusters()
    })
    outputOptions(output, "need_clusters_ui", suspendWhenHidden = FALSE)
    
    observeEvent({
      coords_cache()
      counts_cache()
      input$n_clusters
    }, {
      req(counts_cache(), coords_cache())
      counts <- counts_cache()
      coords <- coords_cache()
      
      # if need_clusters, run clustering
      if (need_clusters()) {
        req(input$n_clusters)
        
        withProgress(message = "Running clustering to predict domains...", {
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
          
          # map to coords
          clusters <- Idents(seurat)
          coords$domain <- new_labels[as.character(clusters[colnames(counts)])]
        })
      }
      
      logfc_cutoff <- if (isTRUE(input$show_advanced)) input$logfc_cutoff else 0.7
      mean_in_cutoff <- if (isTRUE(input$show_advanced)) input$mean_in_cutoff else 1.8
      max_num_gene <- if (isTRUE(input$show_advanced)) input$max_num_gene else 10
      
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
      withProgress(message = "Fitting FG model...", {
        DATA <- estimation_NNGP(DATA, n_neighbors = 10, ORDER = 'AMMD')
      })
      withProgress(message = "Fitting BRISC model...", {
        DATA <- estimation_FGEM(DATA, iter_max = 1000, M_candidates = 2:7, tol = 1e-1)
      })
      data_obj(DATA)
    })
    
    output$summary_ref <- renderPrint({
      req(data_obj())
      obj <- data_obj()
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
      req(data_obj())
      df <- data_obj()@refcolData
      ggplot(df, aes(x=x, y=y, color=as.factor(domain))) +
        geom_point(size=2) + theme_classic() + labs(title="Spatial domains", x="X", y="Y")
    })
    
    return(data_obj)
  })
}
