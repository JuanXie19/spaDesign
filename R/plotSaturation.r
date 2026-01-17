#' Plot Saturation Curves from Sequencing Depth Analysis
#'
#' Visualizes saturation curves produced by \code{\link{saturationDetection}}.
#' The plot overlays observed metric values, the fitted monotone curve, and the detected
#' saturation depth for each group. Optionally adds text labels to display the saturation
#' sequencing depth directly on the plot.
#'
#' @param saturation_results Output list returned by \code{\link{saturationDetection}}.
#' @param metric_col Character string specifying the metric name to display on the y-axis.
#'   Default is \code{"NMI"}.
#' @param use_absolute Logical. If TRUE and absolute depth is available in the results,
#'   the x-axis uses absolute sequencing depth. Otherwise, relative depth is used.
#'   Default is TRUE.
#' @param facet Logical. If TRUE, plots one curve per group using faceting.
#'   If FALSE, overlays all curves on one panel. Default is TRUE.
#' @param color_by_group Logical. If TRUE and multiple groups exist, curves are colored by group.
#'   Default is TRUE.
#' @param show_points Logical. If TRUE, observed data points are plotted. Default is TRUE.
#' @param alpha Numeric scalar in (0,1]. Transparency for points. Default is 0.5.
#' @param label_saturation Logical. If TRUE, add text labels showing saturation depth.
#'   Default is TRUE.
#' @param label_digits Integer. Number of digits used when rounding saturation depth.
#'   Default is 2.
#' @param label_prefix Character string prefix for the label text. Default is \code{"Sat = "}.
#'
#' @return A \code{ggplot2} object.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' res <- saturationDetection(rst, group_cols = "effect_size", aggregate_reps = TRUE)
#'
#' # Faceted plot with saturation labels
#' p1 <- plotSaturation(res, facet = TRUE, label_saturation = TRUE)
#' print(p1)
#'
#' # Overlay plot with saturation labels
#' p2 <- plotSaturation(res, facet = FALSE, color_by_group = TRUE, label_saturation = TRUE)
#' print(p2)
#' }
plotSaturation <- function(saturation_results,
                            metric_col = "NMI",
                            use_absolute = TRUE,
                            facet = TRUE,
                            color_by_group = TRUE,
                            show_points = TRUE,
                            alpha = 0.5,
                            label_saturation = TRUE,
                            label_digits = 2,
                            label_prefix = "Sat = ") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Package 'rlang' is required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
  
  summary_df <- saturation_results$summary
  pred_list  <- saturation_results$predictions
  data_list  <- saturation_results$data
  
  pred_df <- dplyr::bind_rows(lapply(names(pred_list), function(g) {
    x <- pred_list[[g]]
    if (is.null(x)) return(NULL)
    x$group <- g
    x
  }))
  
  obs_df <- dplyr::bind_rows(lapply(names(data_list), function(g) {
    x <- data_list[[g]]
    if (is.null(x)) return(NULL)
    x$group <- g
    x
  }))
  
  if (nrow(pred_df) == 0) stop("No prediction curves found in saturation_results$predictions.")
  
  # ---- choose x-axis ----
  has_abs <- "absolute_depth" %in% colnames(pred_df)
  if (use_absolute && has_abs) {
    x_pred <- "absolute_depth"
    x_obs  <- if ("absolute_depth" %in% colnames(obs_df)) "absolute_depth" else NULL
    xlab <- "Absolute Sequencing Depth"
  } else {
    x_pred <- if ("seq_depth" %in% colnames(pred_df)) "seq_depth" else {
      candidate_cols <- setdiff(colnames(pred_df), c("metric_pred", "absolute_depth", "group"))
      num_cols <- candidate_cols[sapply(pred_df[candidate_cols], is.numeric)]
      if (length(num_cols) == 0) stop("Could not infer depth column in prediction data.")
      num_cols[1]
    }
    x_obs <- x_pred
    xlab <- "Sequencing Depth"
  }
  
  # ---- saturation x column ----
  sat_df <- summary_df
  sat_df$group <- sat_df$group
  
  if (use_absolute && "saturation_absolute_depth" %in% colnames(sat_df)) {
    sat_x_col <- "saturation_absolute_depth"
  } else {
    sat_x_col <- "saturation_depth"
  }
  
  sat_df <- sat_df[!is.na(sat_df[[sat_x_col]]), , drop = FALSE]
  
  n_groups <- length(unique(pred_df$group))
  
  # ---- aesthetics: separate curves properly ----
  if (n_groups > 1) {
    line_aes <- if (color_by_group) {
      ggplot2::aes(x = .data[[x_pred]], y = .data[["metric_pred"]], group = group, color = group)
    } else {
      ggplot2::aes(x = .data[[x_pred]], y = .data[["metric_pred"]], group = group)
    }
    
    point_aes <- if (color_by_group) {
      ggplot2::aes(x = .data[[x_obs]], y = .data[[metric_col]], group = group, color = group)
    } else {
      ggplot2::aes(x = .data[[x_obs]], y = .data[[metric_col]], group = group)
    }
  } else {
    line_aes  <- ggplot2::aes(x = .data[[x_pred]], y = .data[["metric_pred"]])
    point_aes <- ggplot2::aes(x = .data[[x_obs]], y = .data[[metric_col]])
  }
  
  p <- ggplot2::ggplot()
  
  if (show_points && !is.null(x_obs)) {
    p <- p + ggplot2::geom_point(data = obs_df, mapping = point_aes, alpha = alpha)
  }
  
  p <- p + ggplot2::geom_line(data = pred_df, mapping = line_aes, linewidth = 1)
  
  # ---- dashed vertical line ----
  if (nrow(sat_df) > 0) {
    if (n_groups > 1 && color_by_group) {
      p <- p + ggplot2::geom_vline(
        data = sat_df,
        ggplot2::aes(xintercept = .data[[sat_x_col]], color = group),
        linetype = "dashed",
        linewidth = 0.7,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_vline(
        data = sat_df,
        ggplot2::aes(xintercept = .data[[sat_x_col]]),
        linetype = "dashed",
        linewidth = 0.7
      )
    }
  }
  
  # ---- add text label for saturation depth ----
  if (label_saturation && nrow(sat_df) > 0) {
    
    # place the label near the top of the curve (group-specific if possible)
    y_top <- pred_df |>
      dplyr::group_by(group) |>
      dplyr::summarise(y = max(metric_pred, na.rm = TRUE), .groups = "drop")
    
    label_df <- dplyr::left_join(sat_df, y_top, by = "group")
    label_df$label <- paste0(label_prefix, round(label_df[[sat_x_col]], label_digits))
    
    # small y-offset so label doesn't sit on the curve
    label_df$y <- pmin(label_df$y + 0.02, 1)
    
    if (n_groups > 1 && color_by_group) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = .data[[sat_x_col]], y = y, label = label, color = group),
        size = 3.5,
        direction = "x",
        segment.color = NA,
        show.legend = FALSE
      )
    } else {
      p <- p + ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = .data[[sat_x_col]], y = y, label = label),
        size = 3.5,
        direction = "x",
        segment.color = NA,
        show.legend = FALSE
      )
    }
  }
  
  # ---- facet ----
  if (facet && n_groups > 1) {
    p <- p + ggplot2::facet_wrap(~group, scales = "free_x")
  }
  
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::labs(x = xlab, y = metric_col)
  
  if (n_groups > 1 && color_by_group) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Group"))
  }
  
  p <- p + ggplot2::coord_cartesian(ylim = c(0, 1.05))
  return(p)
}

