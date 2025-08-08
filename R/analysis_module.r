## add a function to detect the saturation point
detect_saturation <- function(model, threshold = 0.005, consecutive = 3) {
  # Create fine grid for predictions
  new_data <- data.frame(
    real_seq_depth = seq(min(model$model$real_seq_depth),
                         max(model$model$real_seq_depth),
                         length.out = 200)
  )
  
  # Predicted NMI
  new_data$NMI_pred <- predict(model, new_data)
  
  # Approximate derivative
  new_data$deriv <- c(NA, diff(new_data$NMI_pred) / diff(new_data$real_seq_depth))
  
  # Find first index where slope stays below threshold for 'consecutive' points
  below_thresh <- new_data$deriv < threshold
  run_len <- rle(below_thresh)
  consec_idx <- which(run_len$values & run_len$lengths >= consecutive)
  
  if (length(consec_idx) == 0) {
    return(NA)
  }
  
  # Map back to position in new_data
  start_pos <- sum(run_len$lengths[seq_len(consec_idx[1] - 1)]) + 1
  sat_depth <- new_data$real_seq_depth[start_pos]
  return(sat_depth)
}



analysisUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(

	 # effec size checkbox
      div(
        class = "form-group shiny-input-container", # Important for styling
        tags$label(
          # Input element for the checkbox
          tags$input(id = ns("add_effect_size"), type = "checkbox", class = "shiny-input-checkbox"),
          # The checkbox text label
          tags$span("Change effect size?"),
          # The icon right after the text, but within the label
          tags$span(
            id = ns("effect_size_info"),
            icon("info-circle")
          )
        ),
        # You still need the hidden value for the checkbox, Shiny handles this normally
        # For simplicity, we can let shinyBS handle the tooltip here
        bsTooltip(
          id = ns("effect_size_info"),
          title = "If checked, simulates data with an effect size relative to the baseline. A value of 1 means baseline effect, while 2 means twice the baseline effect. Larger values lead to a stronger signal.",
          placement = "right",
          options = list(container = "body")
        )
      ),
      # Existing helpText can be removed or repositioned if you use this tooltip
      # helpText("If checked, simulates data with modified effect size."), # Remove this if tooltip is enough
      
	  conditionalPanel(
        condition = sprintf("input['%s']", ns("add_effect_size")),
        sliderInput(ns("effect_size"), "Effect size:", min = 1, max = 3, value = 1.5, step = 0.5)
      ),
	  
	  ## add spatial checkbox
	  div(
        class = "form-group shiny-input-container", # Important for styling
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
          title = "If checked, simulate data with spatial noise. Larger sigma leads to more disturbed spatial pattern. ",
          placement = "right",
          options = list(container = "body")
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s']", ns("add_spatial")),
        sliderInput(ns("sigma"), "Sigma (spatial variance):", min = 0.5, max = 3, value = 1.5, step = 0.5)
      ),
      actionButton(ns("run_sim"), "Run Simulation", class = "btn-primary")
    ),
    mainPanel(
      plotlyOutput(ns("sim_plot")),
      verbatimTextOutput(ns("sim_summary")),
	  # add a table summary of the performance vs sequencing depth
	  tags$h4("Detailed simulation results"), # Add a clear heading for the table
      DT::DTOutput(ns("sim_table")) # Placeholder for the interactive table
    )
  )
}

analysisServer <- function(id, data_obj) {
  moduleServer(id, function(input, output, session) {
    sim_results <- reactiveVal(list())
    original_seq_depth <- reactiveVal(NULL)
    
    observeEvent(input$run_sim, {
      req(data_obj())
      #seq_range <- seq(from = input$seq_depth_range[1], to = input$seq_depth_range[2], length.out = 8)
      original_seq_depth_val <- sum(data_obj()@refCounts) / 1e6
      original_seq_depth(original_seq_depth_val)
      
      seq_range <- c(0.5, 1:3, 5, 7, 10)
      withProgress(message = "Running simulations...", {
        base_res <- powerAnalysisEffectSize(data_obj(), es_range = 1, seq_depth_range = seq_range, n_rep = 5)
        base_res$condition <- 'Baseline'
		    results <- list(base = base_res)
        
        if (input$add_effect_size) {
          effect_res <- powerAnalysisEffectSize(data_obj(), es_range = input$effect_size,
                                                seq_depth_range = seq_range, n_rep = 5)
		      effect_res$condition <- 'EffectSize'										
          results$effect <- effect_res
        }
        if (input$add_spatial) {
          spatial_res <- powerAnalysisSpatial(data_obj(),SIGMA = input$sigma, prop_range = 0.7,
                                              seq_depth_range = seq_range, n_rep = 5)
          spatial_res$condition <- 'Spatial'
		      results$spatial <- spatial_res
        }
        sim_results(results)
      })
    })
	  
    # reactive to combine the raw data from all conditions
    combined_raw_data <- reactive({
      results <- sim_results()
      req(length(results) > 0)
      
      req(original_seq_depth())
      bind_rows(results)   %>% 
        mutate(real_seq_depth = seq_depth * original_seq_depth())  
    })
      
	
	  processed_plot_data <- reactive({
	    data_to_process <- combined_raw_data()
	    req(data_to_process)
	    
	    data_to_process %>%
	      group_by(real_seq_depth, condition) %>%
	      summarise(mean_NMI = mean(NMI),
	                se_NMI = sd(NMI) / sqrt(n()),
	                .group = 'drop') %>%
	      mutate(mean_NMI = round(mean_NMI, 3),
	             se_NMI = round(se_NMI, 3),
	             real_seq_depth = round(real_seq_depth, 3))
	  })
	  
	  saturation_points <- list()
	  
	  output$sim_plot <- renderPlotly({
	    plot_data_summary <- processed_plot_data()
	    
	    validate(
	      need(nrow(plot_data_summary) > 0 , "Press 'Run Simulation' to see the results.")
	    )
	    
	    p <- ggplot(plot_data_summary, aes(x = real_seq_depth, y = mean_NMI, color = condition, 
	                                       # Add custom tooltip text
	                                       text = paste("Condition:", condition, "<br>",
	                                                    "Seq. Depth (million):", real_seq_depth, "<br>",
	                                                    "Mean NMI:", mean_NMI))) +
	      geom_point() +
	      geom_errorbar(aes(ymin = mean_NMI - se_NMI, ymax = mean_NMI + se_NMI), 
	                    width = 0.1) +
	      labs(title = 'Mean NMI vs. Sequencing Depth', x = 'Sequencing depth (million)', y = 'Mean NMI') +
	      ylim(0, 1) +
	      theme_minimal()
	    
	    plot_data_raw <- combined_raw_data()
	    for (cond in unique(plot_data_raw$condition)){
	      df_subset <- plot_data_raw %>% filter (condition == cond)
	      
	      # fit the scam model
	      scam_model <- scam(NMI ~ s(real_seq_depth, bs = 'mpi', k = 6), data = df_subset)
	      
	      # Create new data for a smooth prediction line
	      new_data <- data.frame(real_seq_depth = seq(min(df_subset$real_seq_depth), max(df_subset$real_seq_depth), length.out = 100))
	      new_data$NMI_pred <- predict(scam_model, new_data)
	      new_data$condition <- cond
	      p <- p + geom_line(data = new_data, aes(y = NMI_pred, text = NULL), linewidth = 1)
	      
	      ## saturation detection
	      sat_depth <- detect_saturation(scam_model, threshold = 0.005, consecutive = 3)
	      saturation_points[[cond]] <- sat_depth
	      
	      # Add dashed vertical line if found
	      if (!is.na(sat_depth)) {
	        p <- p + geom_vline(xintercept = sat_depth, linetype = "dashed", color = "grey40") +
	          annotate("text", x = sat_depth, y = 0.05, 
	                   label = paste("Sat:", round(sat_depth, 2)),
	                   angle = 90, vjust = -0.5, hjust = 0, size = 3)
	      }
	    }
	    ggplotly(p, tooltip = "text")
	    
	  })

    output$sim_summary <- renderPrint({
      results <- sim_results()
      validate(
        need(length(results) > 0, "No simulations run yet.")
      )
      cat("Simulation summary:\n\n")
      if (!is.null(results$base)) cat("- Base curve (es=1) done\n")
      if (!is.null(results$effect)) cat(sprintf("- Effect curve (es=%.1f) done\n", input$effect_size))
      if (!is.null(results$spatial)) cat(sprintf("- Spatial curve (sigma=%.1f) done\n", input$sigma))
    })
	
	## add data table
	output$sim_table <- DT:: renderDT({
		table_data <- processed_plot_data()
		
		validate(
		  need(nrow(table_data) > 0, "No data available for the table. Run a simulation first.")
		)
		
		table_data <- table_data %>% dplyr::select(c(real_seq_depth, mean_NMI, se_NMI,condition))
		table_data <- table_data %>% dplyr::rename('Sequencing depth (millons)' = real_seq_depth)
		
		DT::datatable(
			table_data,
			options = list(
			pageLength = 10, # Number of rows per page
			lengthMenu = c(5, 10, 25, 50), # Options for rows per page
			scrollX = TRUE
			),
			rownames = FALSE, # Don't show R's default row names
			selection = 'none' # Disable row selection if not needed
		)
    }, server = FALSE) # server = FALSE sends all data to browser, good for small/medium tables.
                      # Set to TRUE for very large datasets for performance.
	
	})
}
