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
          tags$span("Add Effect Size?"),
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
        sliderInput(ns("effect_size"), "Effect size:", min = 1, max = 3, value = 1.5, step = 0.1)
      ),
	  
	  ## add spatial checkbox
	  div(
        class = "form-group shiny-input-container", # Important for styling
        tags$label(
          tags$input(id = ns("add_spatial"), type = "checkbox", class = "shiny-input-checkbox"),
          tags$span("Add Spatial Perturbation?"),
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
        sliderInput(ns("sigma"), "Sigma (spatial variance):", min = 0.5, max = 3, value = 1.5, step = 0.1)
      ),
      numericInput(ns("n_rep"), "Number of replicates:", value = 5, min = 1, step = 1),
      actionButton(ns("run_sim"), "Run Simulation", class = "btn-primary")
    ),
    mainPanel(
      plotOutput(ns("sim_plot")),
      verbatimTextOutput(ns("sim_summary")),
	  
	  # add a table summary of the performance vs sequencing depth
	  tags$h4("Detailed Simulation Results"), # Add a clear heading for the table
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
      seq_range <- c(1:3, 5,7,10)
	  
      withProgress(message = "Running simulations...", {
        base_res <- powerAnalysisEffectSize(data_obj(), es_range = 1, seq_depth_range = seq_range, n_rep = input$n_rep)
        base_res$condition <- 'Base'
		results <- list(base = base_res)
        
        if (input$add_effect_size) {
          effect_res <- powerAnalysisEffectSize(data_obj(), es_range = input$effect_size,
                                                seq_depth_range = seq_range, n_rep = input$n_rep)
		  effect_res$condition <- 'Effect size'										
          results$effect <- effect_res
        }
        if (input$add_spatial) {
          spatial_res <- powerAnalysisSpatial(data_obj(),sigma = input$sigma, prop_range = 0.7,
                                              seq_depth_range = seq_range, n_rep = input$n_rep)
          spatial_res$condition <- 'Spatial'
		  results$spatial <- spatial_res
        }
        sim_results(results)
      })
    })
	
	
	processed_plot_data <- reactive({
      req(sim_results()) # Ensure simulations have run
      results <- sim_results()

      ## helper to compute summary + fit
      process_curve <- function(df, label){
        df_sum <- df %>%
          group_by(seq_depth) %>%
          summarise(mean_NMI = mean(NMI),
                    se_NMI = sd(NMI) / sqrt (n()),
                    .groups = 'drop') %>%
          mutate(condition = label,
				mean_NMI = round(mean_NMI,3),
				se_NMI = round(se_NMI, 3))
      }

      df_base <- process_curve(results$base, "Base")
      df_base <- as.data.frame(df_base)
      list_curves <- list(df_base)

      if (!is.null(results$effect)) {
        df_effect <- process_curve(results$effect, "Effect size")
        df_effect <- as.data.frame(df_effect)
        list_curves <- append(list_curves, list(df_effect))
      }
      if (!is.null(results$spatial)) {
        df_spatial <- process_curve(results$spatial, "Spatial")
        df_spatial <- as.data.frame(df_spatial)
        list_curves <- append(list_curves, list(df_spatial))
      }

      # Combine for ggplot and table
      bind_rows(list_curves)
    })
    
    output$sim_plot <- renderPlot({
      req(processed_plot_data())
	  plot_data <- processed_plot_data()
	  
	  library(ggplot2)
	  
	  ggplot(plot_data, aes(x = seq_depth, y = mean_NMI, color = condition)) + 
		geom_point() + 
		geom_smooth(method = 'loess', span = 0.7, se = FALSE) + ylim(c(0, 1))+
		theme_minimal()
	
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

