#' Detect Saturation Point in Sequencing Depth Analysis
#'
#' This function fits a monotonic spline model to analyze the relationship between
#' sequencing depth and a performance metric (e.g., NMI), then identifies the saturation
#' point where further sequencing depth yields diminishing returns.
#'
#' @param data A data frame containing sequencing depth and performance metric information.
#' @param nmi_col Character string specifying the column name for the performance 
#'   metric (e.g., NMI, ARI). Default is "NMI".
#' @param depth_col Character string specifying the column name for sequencing depth. 
#'   Default is "seq_depth".
#' @param pilot_depth Numeric value for pilot sequencing depth. If provided, results 
#'   will show absolute sequencing depth (seq_depth * pilot_depth). If NULL (default), 
#'   results show relative sequencing depth.
#' @param slope_threshold Numeric value defining the slope threshold below which a 
#'   point is considered saturated. Default is 0.05.
#' @param required_NMI_percentage Numeric value between 0 and 1 indicating the 
#'   percentage of maximum NMI that must be achieved before considering saturation. 
#'   Default is 0.8 (80%).
#' @param plot Logical indicating whether to generate a plot. Default is TRUE.
#'
#' @return A list with the following components:
#'   \item{saturation_point}{A data frame containing the detected saturation point, 
#'     including sequencing depth and predicted NMI value.}
#'   \item{predictions}{A data frame with model predictions across the full range 
#'     of sequencing depths.}
#'   \item{model}{The fitted SCAM (Shape Constrained Additive Model) object.}
#'   \item{data}{The input data used for model fitting.}
#'   \item{plot}{A ggplot2 object (if plot = TRUE) showing the data points, fitted 
#'     curve, and saturation point.}
#'
#' @details
#' The function uses Shape Constrained Additive Models (SCAM) to fit a monotonically
#' increasing spline to the relationship between sequencing depth and performance
#' metric. Saturation is detected when:
#' \enumerate{
#'   \item The predicted metric reaches the specified percentage of its maximum value
#'   \item The slope of the curve falls below the specified threshold
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with relative depth
#' results <- detect_saturation(
#'   data = my_results,
#'   slope_threshold = 0.05,
#'   required_NMI_percentage = 0.8
#' )
#'
#' # With absolute depth conversion
#' results <- detect_saturation(
#'   data = my_results,
#'   pilot_depth = 19.3,
#'   slope_threshold = 0.05,
#'   required_NMI_percentage = 0.8
#' )
#'
#' # Access saturation point
#' print(results$saturation_point)
#'
#' # View the plot
#' results$plot
#' }
#'
#' @import scam
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @export
detect_saturation <- function(data,
                              nmi_col = "NMI",
                              depth_col = "seq_depth",
                              pilot_depth = NULL,
                              slope_threshold = 0.05,
                              required_NMI_percentage = 0.8,
                              plot = TRUE) {
  
  # Load required packages
  require(scam)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  
  # Add absolute depth if pilot_depth provided
  if (!is.null(pilot_depth)) {
    data$absolute_depth <- data[[depth_col]] * pilot_depth
  }
  
  # Fit SCAM model (simple monotonic spline, no factor interaction)
  formula_str <- paste0(nmi_col, " ~ s(", depth_col, ", bs = 'mpi', k = 8)")
  scam_model <- scam(as.formula(formula_str), data = data)
  
  # Create prediction grid
  new_seq_depth <- seq(0, max(data[[depth_col]]), length.out = 200)
  pred_grid <- data.frame(seq_depth = new_seq_depth)
  names(pred_grid)[1] <- depth_col
  
  # Make predictions
  pred_grid$NMI_pred <- predict(scam_model, newdata = pred_grid)
  
  # Add absolute depth to predictions if pilot_depth provided
  if (!is.null(pilot_depth)) {
    pred_grid$absolute_depth <- pred_grid[[depth_col]] * pilot_depth
  }
  
  # Detect saturation point
  saturation_point <- pred_grid %>%
    arrange(!!sym(depth_col)) %>%
    mutate(
      max_NMI_pred = max(NMI_pred),
      diff_NMI = NMI_pred - lag(NMI_pred),
      diff_depth = !!sym(depth_col) - lag(!!sym(depth_col)),
      slope = diff_NMI / diff_depth,
      NMI_threshold = max_NMI_pred * required_NMI_percentage
    ) %>%
    filter(!is.na(slope)) %>%
    filter(!!sym(depth_col) >= min(data[[depth_col]])) %>%
    filter(NMI_pred >= NMI_threshold) %>%
    filter(slope < slope_threshold) %>%
    slice(1) %>%
    select(!!sym(depth_col), NMI_pred)
  
  # Handle case where no saturation point is found
  if (nrow(saturation_point) == 0) {
    saturation_point <- pred_grid %>%
      filter(!!sym(depth_col) == max(!!sym(depth_col))) %>%
      select(!!sym(depth_col), NMI_pred)
    warning("No saturation point found meeting the criteria. Returning maximum depth.")
  }
  
  # Add absolute depth to saturation point if pilot_depth provided
  if (!is.null(pilot_depth)) {
    saturation_point$absolute_depth <- saturation_point[[depth_col]] * pilot_depth
  }
  
  # Prepare output
  results <- list(
    saturation_point = saturation_point,
    predictions = pred_grid,
    model = scam_model,
    data = data
  )
  
  # Generate plot if requested
  if (plot) {
    # Determine which depth column to use for plotting
    x_label <- if (!is.null(pilot_depth)) "Absolute Sequencing Depth" else "Sequencing Depth"
    
    # Prepare data for plotting
    plot_data <- data
    plot_pred <- pred_grid
    plot_sat <- saturation_point
    
    if (is.null(pilot_depth)) {
      plot_data$plot_depth <- plot_data[[depth_col]]
      plot_pred$plot_depth <- plot_pred[[depth_col]]
      plot_sat$plot_depth <- plot_sat[[depth_col]]
    } else {
      plot_data$plot_depth <- plot_data$absolute_depth
      plot_pred$plot_depth <- plot_pred$absolute_depth
      plot_sat$plot_depth <- plot_sat$absolute_depth
    }
    
    p <- ggplot() +
      geom_point(data = plot_data, 
                 aes(x = plot_depth, y = !!sym(nmi_col)), 
                 alpha = 0.5, color = "steelblue") +
      geom_line(data = plot_pred, 
                aes(x = plot_depth, y = NMI_pred), 
                linewidth = 1, color = "steelblue") +
      geom_vline(data = plot_sat, 
                 aes(xintercept = plot_depth), 
                 linetype = "dashed", color = "red") +
      geom_text_repel(data = plot_sat, 
                      aes(x = plot_depth, y = 0.2, 
                          label = round(plot_depth, 2)),
                      color = "red",
                      direction = "x",
                      nudge_x = 0.5,
                      segment.color = NA,
                      show.legend = FALSE) +
      coord_cartesian(ylim = c(0, 1.05)) +
      theme_classic() +
      labs(x = x_label, y = "NMI")
    
    results$plot <- p
    print(p)
  }
  
  return(invisible(results))
}

