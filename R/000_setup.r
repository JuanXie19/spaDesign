#' @keywords internal
#' @importFrom graphics plot.new text
#' @import hdf5r
#' @importFrom dplyr %>% n
#' @importFrom ggplot2 aes element_text geom_errorbar geom_line geom_point geom_text geom_vline ggplot labs theme theme_classic theme_minimal ylim
#' @importFrom plotly ggplotly
#' @importFrom plotly renderPlotly
#' @importFrom rlang .data
#' @importFrom shiny NS actionButton br checkboxInput conditionalPanel downloadButton downloadHandler eventReactive fileInput h4 h5 helpText hr mainPanel moduleServer need numericInput observeEvent outputOptions plotOutput radioButtons reactive reactiveVal renderPlot renderPrint req selectInput showNotification sidebarLayout sidebarPanel sliderInput tags updateSelectInput validate verbatimTextOutput withProgress
#' @importFrom stats cutree dist hclust kmeans median rpois sd
#' @importFrom utils head tail write.csv
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


