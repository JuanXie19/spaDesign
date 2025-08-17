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


analysisUI <- function(id){
  ns <- NS(id)
  sidebarLayout(
    # sidebar panel
    sidebarPanel(
      # effect size checkbox
      div(
        class = "form-group shiny-input-container",
        tags$label(
          tags$label(
            tags$input(id = ns("add_effect_size"), type = "checkbox", class = "shiny-input-checkbox"),
            tags$span("Change effect size?"),
            tags$span(id = ns("effect_size_info"), icon ("info-circle"))
          )
        ),
        bsTooltip(
          id = ns("effect_size_info"),
          title = "If checked, simulates data with changed effect size. A value of 1 means original effect size, while 2 means twice the original effect size. Larger values lead to a stronger signal. ",
          placement = "right",
          options = list(container = "body")
        )
      ),
      
      conditionalPanel(
        condition = sprintf("input['%s']", ns("add_effect_size")),
        sliderInput(ns("effect_size"), "Effect size:", min = 1, max = 3, value = 1.5, step = 0.5)
      ),
      
      ## add spatial checkbox
      div(
        class= "form-group shiny-input-container",
        tags$label(
          tags$input(id = ns("add_spatial"), type = "checkbox", class = "shiny-input-checkbox"),
          tags$span("Add spatial perturbation?"),
          tags$span(
            id = ns("spatial_info"),
            icon("info-circle")
          )
        ),
        bsTooltip(
          id = ns("spatial_info"),
          title = "If checked, simulate data with spatial noise. Larger sigma leads to more disturbed spatial pattern.",
          placement = "right",
          options = list(container = "body")
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s']", ns("add_spatial")),
        sliderInput(ns("sigma"), "Sigma(spatial noise):", min = 0.5, max = 3, value = 1, step = 0.5)
      ),
      actionButton(ns("run_sim"), "Run simulation", class = "btn-primary")
    ),
    
    # main panel
    mainPanel(
      plotlyOutput(ns("sim_plot")),
      tags$h4("Detailed simulation results"),
      DT::DTOutput(ns("sim_table"))
    )
  )
}


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
                                            seq_depth_range = seq_range, n_rep = 5)
        base_res$condition <- 'Original'
        results <- list(base = base_res)
        
        if(input$add_effect_size){
          effect_res <- powerAnalysisEffectSize(data_obj(), es_range = input$effect_size,
                                                seq_depth_range = seq_range, n_rep = 5)
          effect_res$condition <- 'EffectSize'
          results$effect <- effect_res
        }
        
        if(input$add_spatial){
          spatial_res <- powerAnalysisSpatial(data_obj(), SIGMA = input$sigma, prop_range = 0.8,
                                              seq_depth_range = seq_range, n_rep = 5)
          spatial_res$condition <- 'Spatial'
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
      
      # get the default colors used by ggplot2
      default_colors <- ggplot_build(p)$data[[1]]$colour
      
      
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
          current_color <- default_colors[[i]]
          p <- p + geom_vline(xintercept = sat_depth, linetype = 'dashed', color = current_color) + 
            annotate("text", x = sat_depth, y = 0.05,
                     label = paste("Sat:", round(sat_depth, 2)),
                     angle = 90, vjust = -0.5, hjust = 0, size = 3, color = current_color)
        }
      }
      saturation_points(current_saturation_points)
      ggplotly(p, tooltip = 'text')
    })
    
    output$sim_table <- DT::renderDT({
      table_data <- processed_plot_data()
      
      validate(
        need(nrow(table_data) > 0, "No data available for the table. Run a simulation first.")
      )
      
      table_data <- table_data %>% dplyr::select(c(real_seq_depth, mean_NMI, se_NMI, condition))
      
      table_data <- table_data %>% dplyr::rename('Sequencing depth (millions)' = real_seq_depth)
      
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
  })
}

