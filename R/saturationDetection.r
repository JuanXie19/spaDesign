#' Detect Saturation Point in Sequencing Depth Analysis
#'
#' Fits a monotone increasing spline model to the relationship between sequencing depth
#' and a performance metric (e.g., NMI), then identifies the saturation point where
#' additional sequencing depth yields diminishing improvement.
#'
#' The function supports two common scenarios:
#' \itemize{
#'   \item \strong{Single curve:} only one scenario is present (e.g., one \code{prop} value)
#'   \item \strong{Stratified curves:} multiple scenarios exist (e.g., multiple \code{effect_size} levels),
#'         and saturation is detected separately for each level using \code{group_cols}
#' }
#'
#' @param data A data frame containing sequencing depth and metric values.
#' @param metric_col Character string specifying the column name for the performance metric.
#'   Default is \code{"NMI"}.
#' @param depth_col Character string specifying the column name for sequencing depth.
#'   Default is \code{"seq_depth"}.
#' @param group_cols Character vector of column names used to stratify the analysis.
#'   If \code{NULL} (default), a single curve is fitted using all rows in \code{data}.
#'   If provided, one curve is fitted per unique group.
#' @param pilot_depth Numeric value for pilot sequencing depth. If provided, absolute depth is
#'   computed as \code{depth * pilot_depth} and stored in the output predictions and summary.
#'   Default is \code{NULL}.
#' @param slope_threshold Numeric scalar > 0. Slope threshold below which a point is considered
#'   saturated. Default is 0.05.
#' @param required_metric_percentage Numeric scalar in (0, 1]. Saturation is only considered after
#'   the predicted metric reaches at least this fraction of the maximum predicted value.
#'   Default is 0.8.
#' @param k Integer. Maximum spline basis dimension for the monotone spline. The function will
#'   automatically reduce \code{k} within each group if the dataset is small. Default is 8.
#' @param grid_size Integer. Number of depth points used for prediction and slope calculation.
#'   Default is 200.
#' @param aggregate_reps Logical. If TRUE, replicate rows at the same depth (within each group)
#'   are averaged before fitting. This can improve stability for small \code{n_rep}. Default is FALSE.
#' @param plot Logical. If TRUE, plots are generated (optional). Default is TRUE.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{summary}: data frame of saturation results (one row per group)
#'     \item \code{predictions}: list of prediction data frames (one per group)
#'     \item \code{models}: list of fitted SCAM models (one per group; may contain NULL)
#'     \item \code{data}: list of input data split by group
#'   }
#'
#' @importFrom scam scam
#' @importFrom dplyr arrange mutate filter bind_rows group_by summarise
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: single curve (no stratification)
#' res <- detect_saturation(
#'   data = rst2,
#'   metric_col = "NMI",
#'   depth_col = "seq_depth",
#'   group_cols = NULL,
#'   aggregate_reps = TRUE
#' )
#' res$summary
#'
#' # Example 2: stratified curves (by effect size)
#' res <- detect_saturation(
#'   data = rst,
#'   metric_col = "NMI",
#'   depth_col = "seq_depth",
#'   group_cols = "effect_size",
#'   aggregate_reps = TRUE
#' )
#' res$summary
#' }
saturationDetection <- function(data,
                              metric_col = "NMI",
                              depth_col = "seq_depth",
                              group_cols = NULL,
                              pilot_depth = NULL,
                              slope_threshold = 0.05,
                              required_metric_percentage = 0.8,
                              k = 8,
                              grid_size = 200,
                              aggregate_reps = FALSE,
                              plot = TRUE) {
  # ...
}

detect_saturation <- function(data,
                              metric_col = "NMI",
                              depth_col = "seq_depth",
                              group_cols = NULL,
                              pilot_depth = NULL,
                              slope_threshold = 0.05,
                              required_metric_percentage = 0.8,
                              k = 8,
                              grid_size = 200,
                              aggregate_reps = FALSE,
                              plot = TRUE) {
  
  if (!requireNamespace("scam", quietly = TRUE)) stop("Package 'scam' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Package 'rlang' is required.")
  
  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!metric_col %in% names(data)) stop("metric_col not found in data.")
  if (!depth_col %in% names(data)) stop("depth_col not found in data.")
  
  if (!is.null(group_cols) && length(group_cols) > 0) {
    missing_cols <- setdiff(group_cols, names(data))
    if (length(missing_cols) > 0) stop("group_cols not found in data: ", paste(missing_cols, collapse = ", "))
  }
  
  # add absolute depth if requested
  if (!is.null(pilot_depth)) {
    data$absolute_depth <- data[[depth_col]] * pilot_depth
  }
  
  # grouping setup
  if (is.null(group_cols) || length(group_cols) == 0) {
    data$.group_id <- "ALL"
  } else {
    data$.group_id <- apply(data[, group_cols, drop = FALSE], 1, paste, collapse = " | ")
  }
  groups <- split(data, data$.group_id)
  
  # helper: fit + saturation for one group
  fit_one_group <- function(df) {
    
    # optional: aggregate replicates within each depth (and within group_cols)
    if (aggregate_reps) {
      df <- df |>
        dplyr::group_by(.data[[depth_col]]) |>
        dplyr::summarise(
          !!metric_col := mean(.data[[metric_col]], na.rm = TRUE),
          .groups = "drop"
        )
      # after aggregation, df has only unique depth rows
    }
    
    n_unique_depth <- length(unique(df[[depth_col]]))
    n_obs <- nrow(df)
    
    # need at least 3 unique depths to meaningfully detect saturation
    if (n_unique_depth < 3) {
      return(list(
        saturation_point = NA_real_,
        saturation_metric = NA_real_,
        model = NULL,
        predictions = NULL,
        warning = "Too few unique depth values to fit saturation curve."
      ))
    }
    
    # --- KEY FIX: shrink k automatically ---
    # scam spline complexity must be small enough for the data
    # safe and conservative rule:
    k_eff <- min(k, n_unique_depth)
    k_eff <- max(3, k_eff)  # allow at least 3 if possible
    k_eff <- min(k_eff, n_obs - 1)  # avoid more coefficients than data
    
    # If still too small, fallback to linear model later
    if (k_eff < 3) {
      return(list(
        saturation_point = NA_real_,
        saturation_metric = NA_real_,
        model = NULL,
        predictions = NULL,
        warning = "Not enough data to fit SCAM. Consider aggregate_reps=TRUE or more depth points."
      ))
    }
    
    formula_str <- paste0(metric_col, " ~ s(", depth_col, ", bs = 'mpi', k = ", k_eff, ")")
    
    model <- scam::scam(stats::as.formula(formula_str), data = df)
    
    # prediction grid
    new_depth <- seq(min(df[[depth_col]]), max(df[[depth_col]]), length.out = grid_size)
    pred_grid <- data.frame(tmp = new_depth)
    names(pred_grid)[1] <- depth_col
    pred_grid$metric_pred <- stats::predict(model, newdata = pred_grid)
    
    if (!is.null(pilot_depth)) {
      pred_grid$absolute_depth <- pred_grid[[depth_col]] * pilot_depth
    }
    
    # saturation detection
    pred2 <- pred_grid |>
      dplyr::arrange(.data[[depth_col]]) |>
      dplyr::mutate(
        max_metric = max(metric_pred, na.rm = TRUE),
        metric_threshold = max_metric * required_metric_percentage,
        diff_metric = metric_pred - dplyr::lag(metric_pred),
        diff_depth = .data[[depth_col]] - dplyr::lag(.data[[depth_col]]),
        slope = diff_metric / diff_depth
      ) |>
      dplyr::filter(!is.na(slope)) |>
      dplyr::filter(metric_pred >= metric_threshold) |>
      dplyr::filter(slope < slope_threshold)
    
    if (nrow(pred2) == 0) {
      sat_depth <- max(pred_grid[[depth_col]])
      sat_metric <- pred_grid$metric_pred[which.max(pred_grid[[depth_col]])]
      warn_msg <- "No saturation point found meeting criteria. Returning maximum depth."
    } else {
      sat_depth <- pred2[[depth_col]][1]
      sat_metric <- pred2$metric_pred[1]
      warn_msg <- ""
    }
    
    list(
      saturation_point = sat_depth,
      saturation_metric = sat_metric,
      model = model,
      predictions = pred_grid,
      warning = warn_msg
    )
  }
  
  out <- lapply(groups, fit_one_group)
  
  summary_df <- data.frame(
    group = names(out),
    saturation_depth = sapply(out, `[[`, "saturation_point"),
    saturation_metric = sapply(out, `[[`, "saturation_metric"),
    warning = sapply(out, `[[`, "warning"),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(pilot_depth)) {
    summary_df$saturation_absolute_depth <- summary_df$saturation_depth * pilot_depth
  }
  
  results <- list(
    summary = summary_df,
    predictions = lapply(out, `[[`, "predictions"),
    models = lapply(out, `[[`, "model"),
    data = groups
  )
  
  return(results)
}
