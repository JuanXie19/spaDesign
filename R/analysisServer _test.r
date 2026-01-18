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
    
    # New reactive for saturation detection results
    saturation_results <- reactive({
      raw_data <- combined_raw_data()
      req(raw_data)
      
      # Run saturation detection with your new function
      saturationDetection(
        data = raw_data,
        metric_col = "NMI",
        depth_col = "seq_depth",
        group_cols = "condition",  # stratify by condition (Original, Effect Size, Spatial)
        pilot_depth = original_seq_depth(),
        slope_threshold = 0.005,
        required_metric_percentage = 0.8,
        k = 7,
        grid_size = 200,
        aggregate_reps = TRUE
      )
    })
    
    output$sim_plot <- renderPlotly({
      sat_results <- saturation_results()
      plot_data_summary <- processed_plot_data()
      
      validate(
        need(nrow(plot_data_summary) > 0, "Press 'Run simulation' to see the results.")
      )
      
      # Create base ggplot with observed data
      p <- ggplot(plot_data_summary, aes(x = real_seq_depth, y = mean_NMI, color = condition,
                                         text = paste("Condition:", condition, "<br>",
                                                      "Seq.Depth (million):", real_seq_depth, "<br>",
                                                      "Mean NMI:", mean_NMI))) +
        geom_point() + 
        geom_errorbar(aes(ymin = mean_NMI - se_NMI, ymax = mean_NMI + se_NMI), width = 0.1) +
        theme_minimal() +
        labs(title = 'Mean NMI vs. Sequencing depth', 
             x = 'Sequencing depth (million)', 
             y = 'Mean NMI') +
        ylim(0, 1)
      
      # Add fitted curves from saturation detection
      pred_df <- dplyr::bind_rows(lapply(names(sat_results$predictions), function(g) {
        x <- sat_results$predictions[[g]]
        if (is.null(x)) return(NULL)
        x$condition <- g
        x
      }))
      
      if (nrow(pred_df) > 0) {
        p <- p + geom_line(
          data = pred_df,
          aes(x = absolute_depth, y = metric_pred, color = condition, text = NULL),
          linewidth = 1
        )
      }
      
      # Add saturation vertical lines
      sat_summary <- sat_results$summary
      sat_summary <- sat_summary[!is.na(sat_summary$saturation_absolute_depth), ]
      
      if (nrow(sat_summary) > 0) {
        # Rename group to condition for consistency
        sat_summary$condition <- sat_summary$group
        
        # Add vertical dashed lines
        p <- p + geom_vline(
          data = sat_summary,
          aes(xintercept = saturation_absolute_depth, color = condition),
          linetype = 'dashed',
          show.legend = FALSE
        )
        
        # Add saturation labels
        # Stagger labels vertically to avoid overlap
        sat_summary$label_y_offset <- (seq_len(nrow(sat_summary)) %% 2) * 0.08
        sat_summary$label_y <- 0.05 + sat_summary$label_y_offset
        sat_summary$label <- paste('Sat:', round(sat_summary$saturation_absolute_depth, 2))
        
        p <- p + geom_text(
          data = sat_summary,
          aes(x = saturation_absolute_depth, y = label_y, label = label, color = condition),
          angle = 90, vjust = -0.5, hjust = 0, size = 3,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
      }
      
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