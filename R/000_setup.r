#' @keywords internal
#'
#' @importFrom shiny navbarPage tabPanel
#' @importFrom utils globalVariables
"_PACKAGE"

# Define the global variable in the package namespace
reference_data_paths <- NULL

utils::globalVariables(c(
  ":=", "absolute_depth", "condition", "diff_depth", "diff_metric", "domain",
  "group", "label", "label_y", "logFC_low", "max_metric", "mean_NMI",
  "metric_pred", "metric_threshold", "NMI", "real_seq_depth",
  "saturation_absolute_depth", "se_NMI", "seq_depth", "slope", "x", "y"
))

.onLoad <- function(libname, pkgname) {
  
  # Set Shiny upload limit
  options(shiny.maxRequestSize = 300 * 1024^2)
  
  # Store file paths inside the package namespace
  reference_data_paths <<- list(
    "Chicken Heart" = system.file("extdata/ref_chickenHeart.rds", package = pkgname),
    "Human Brain"   = system.file("extdata/ref_humanBrain.rds", package = pkgname)
  )
}


