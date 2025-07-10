analysisUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      sliderInput(ns("seq_depth_range"), "Sequencing depth range:",
                  min = 0.5, max = 2, value = c(0.5, 1.5), step = 0.1),
      checkboxInput(ns("add_effect_size"), "Add Effect Size?", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s']", ns("add_effect_size")),
        sliderInput(ns("effect_size"), "Effect size:", min = 1, max = 3, value = 1.5, step = 0.1)
      ),
      checkboxInput(ns("add_spatial"), "Add Spatial Perturbation?", value = FALSE),
      conditionalPanel(
        condition = sprintf("input['%s']", ns("add_spatial")),
        sliderInput(ns("sigma"), "Sigma (spatial variance):", min = 0.5, max = 3, value = 1.5, step = 0.1)
      ),
      numericInput(ns("n_rep"), "Number of replicates:", value = 5, min = 1, step = 1),
      actionButton(ns("run_sim"), "Run Simulation", class = "btn-primary")
    ),
    mainPanel(
      plotOutput(ns("sim_plot")),
      verbatimTextOutput(ns("sim_summary"))
    )
  )
}

analysisServer <- function(id, data_obj) {
  moduleServer(id, function(input, output, session) {
    sim_results <- reactiveVal(list())
    
    observeEvent(input$run_sim, {
      req(data_obj())
      seq_range <- seq(from = input$seq_depth_range[1], to = input$seq_depth_range[2], length.out = 8)
      
      withProgress(message = "Running simulations...", {
        base_res <- powerAnalysisEffectSize(data_obj(), es_range = 1, seq_depth_range = seq_range, n_rep = input$n_rep)
        results <- list(base = base_res)
        
        if (input$add_effect_size) {
          effect_res <- powerAnalysisEffectSize(data_obj(), es_range = input$effect_size,
                                                seq_depth_range = seq_range, n_rep = input$n_rep)
          results$effect <- effect_res
        }
        if (input$add_spatial) {
          spatial_res <- powerAnalysisSpatial(data_obj(), prop_range = 0.5,
                                              seq_depth_range = seq_range, n_rep = input$n_rep)
          results$spatial <- spatial_res
        }
        sim_results(results)
      })
    })
    
    output$sim_plot <- renderPlot({
      req(sim_results())
      results <- sim_results()
      p <- ggplot(results$base, aes(x=seq_depth, y=NMI)) +
        geom_point(color="blue") + geom_line(color="blue") +
        labs(title="Performance vs Sequencing Depth", x="Sequencing Depth", y="NMI") + theme_minimal()
      if (!is.null(results$effect)) {
        p <- p + geom_point(data=results$effect, aes(x=seq_depth, y=NMI), color="red") +
          geom_line(data=results$effect, aes(x=seq_depth, y=NMI), color="red")
      }
      if (!is.null(results$spatial)) {
        p <- p + geom_point(data=results$spatial, aes(x=seq_depth, y=NMI), color="green") +
          geom_line(data=results$spatial, aes(x=seq_depth, y=NMI), color="green")
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
  })
}
