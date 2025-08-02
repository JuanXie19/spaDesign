analysisUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(

      #sliderInput(ns("seq_depth_range"), "Sequencing depth range:",
      #            min = 0.5, max = 2, value = c(0.5, 1.5), step = 0.1),
	 
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
      plotOutput(ns("sim_plot")),
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
    
    observeEvent(input$run_sim, {
      req(data_obj())
      #seq_range <- seq(from = input$seq_depth_range[1], to = input$seq_depth_range[2], length.out = 8)

      seq_range <- c(0.5, 1:3, 5,7,10)
      withProgress(message = "Running simulations...", {
        base_res <- powerAnalysisEffectSize(data_obj(), es_range = 1, seq_depth_range = seq_range, n_rep = 10)
        base_res$condition <- 'Baseline'
		    results <- list(base = base_res)
        
        if (input$add_effect_size) {
          effect_res <- powerAnalysisEffectSize(data_obj(), es_range = input$effect_size,
                                                seq_depth_range = seq_range, n_rep = 10)
		      effect_res$condition <- 'EffectSize'										
          results$effect <- effect_res
        }
        if (input$add_spatial) {
          spatial_res <- powerAnalysisSpatial(data_obj(),SIGMA = input$sigma, prop_range = 0.7,
                                              seq_depth_range = seq_range, n_rep = 10)
          spatial_res$condition <- 'Spatial'
		      results$spatial <- spatial_res
        }
        sim_results(results)
      })
    })
	  
    # reactive to combine the raw data from all conditions
    combined_raw_data <- reactive({
      req(sim_results())
      results <- sim_results()
      bind_rows(results)    
    })
      
	
	  processed_plot_data <- reactive({
      req(combined_raw_data()) # Ensure simulations have run
      
	    combined_raw_data() %>%
	      group_by(seq_depth, condition) %>%
	      summarise(mean_NMI = mean(NMI),
	                se_NMI = sd(NMI) / sqrt(n()),
	                .group = 'droup') %>%
	      mutate(mean_NMI = round(mean_NMI, 3),
	             se_NMI = round(se_NMI, 3))
	  })
	  
	  output$sim_plot <- renderPlot({
	    req(combined_raw_data())
	    plot_data_raw <- combined_raw_data()
	    
	    print("Plotting with data:")
	    print(head(plot_data_raw))
	    print(paste("Number of rows:", nrow(plot_data_raw)))
	    print(paste("Unique conditions:", paste(unique(plot_data_raw$condition), collapse = ", ")))
	    
	    # 
	    p <- ggplot(plot_data_raw, aes(x = seq_depth, y = NMI, color = condition)) + 
	      geom_point(alpha = 0.5) + 
	      labs(title = '', x = 'Sequencing depth', y = 'NMI') + 
	      ylim(0, 1) +
	      theme_minimal()
	    
	    for (cond in unique(plot_data_raw$condition)){
	      df_subset <- plot_data_raw %>% filter (condition == cond)
	      
	      # fit the scam model
	      scam_model <- scam(NMI ~ s(seq_depth, bs = 'mpi', k = 6), data = df_subset)
	      
	      # Create new data for a smooth prediction line
	      new_data <- data.frame(seq_depth = seq(min(df_subset$seq_depth), max(df_subset$seq_depth), length.out = 100))
	      new_data$NMI_pred <- predict(scam_model, new_data)
	      new_data$condition <- cond
	      
	      # Add the non-decreasing curve to the plot
	      p <- p + geom_line(data = new_data, aes(y = NMI_pred), linewidth = 1)
	    }
	    p  
	  })


    output$sim_summary <- renderPrint({
      req(sim_results())
      cat("Simulation summary:\n\n")
      if (!is.null(sim_results()$base)) cat("- Base curve (es=1) done\n")
      if (!is.null(sim_results()$effect)) cat(sprintf("- Effect curve (es=%.1f) done\n", input$effect_size))
      if (!is.null(sim_results()$spatial)) cat(sprintf("- Spatial curve (sigma=%.1f) done\n", input$sigma))
    })
	
	## add data table
	output$sim_table <- DT:: renderDT({
		req(processed_plot_data())
		table_data <- processed_plot_data()
		
		DT::datatable(
			table_data,
			options = list(
			pageLength = 10, # Number of rows per page
			lengthMenu = c(5, 10, 25, 50), # Options for rows per page
			scrollX = TRUE # Enable horizontal scrolling if table is wide
			),
			rownames = FALSE, # Don't show R's default row names
			selection = 'none' # Disable row selection if not needed
		)
    }, server = FALSE) # server = FALSE sends all data to browser, good for small/medium tables.
                      # Set to TRUE for very large datasets for performance.
	
	})
  }

