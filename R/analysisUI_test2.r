#' Analysis UI Module
#'
#' UI for simulation and analysis tab. Allows effect size and spatial perturbation simulations.
#'
#' @param id Module namespace.
#' @importFrom plotly plotlyOutput
#' @importFrom shinyBS bsTooltip
#' @return Shiny UI element.
#' @export
analysisUI <- function(id){
  ns <- NS(id)
  sidebarLayout(
    # sidebar panel
    sidebarPanel(
      radioButtons(ns("simulation_type"), "Select simulation modification:",
                   choices = c("None (run original only)" = "original",
                               "Change effect size" = "effect_size",
                               "Add spatial perturbation" = "spatial"),
                   selected = "original"),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'effect_size'", ns("simulation_type")),
        sliderInput(ns("effect_size"), "Effect size:", min = 1, max = 3, value = 1.5, step = 0.5),
        shinyBS::bsTooltip(id = ns("effect_size"), 
                           title = "A value of 1 means original effect size, while 2 means twice the original effect size. Larger values lead to a stronger signal.", 
                           placement = "right")
      ),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'spatial'", ns("simulation_type")),
        sliderInput(ns("sigma"), "Sigma (spatial noise):", min = 0.5, max = 3, value = 1, step = 0.5),
        shinyBS::bsTooltip(id = ns("sigma"), 
                           title = "Larger sigma leads to more disturbed spatial pattern.", 
                           placement = "right")
      ),
      
      # Add advanced saturation detection options
      checkboxInput(ns("show_saturation_advanced"), "Advanced saturation detection options", value = FALSE),
      
      conditionalPanel(
        condition = sprintf("input['%s']", ns("show_saturation_advanced")),
        sliderInput(ns("slope_threshold"), 
                    "Slope threshold:", 
                    min = 0.001, max = 0.1, value = 0.005, step = 0.001),
        shinyBS::bsTooltip(id = ns("slope_threshold"), 
                           title = "Slope below this value indicates saturation. Lower values = stricter saturation criteria.", 
                           placement = "right"),
        
        sliderInput(ns("metric_percentage"), 
                    "Required metric percentage:", 
                    min = 0.5, max = 1.0, value = 0.8, step = 0.05),
        shinyBS::bsTooltip(id = ns("metric_percentage"), 
                           title = "Saturation is only considered after metric reaches this fraction of maximum. Higher values = more conservative.", 
                           placement = "right")
      ),
      
      actionButton(ns("run_sim"), "Run simulation", class = "btn-primary")
    ),
    
    # main panel
    mainPanel(
      plotly::plotlyOutput(ns("sim_plot")),
      tags$h4("Detailed simulation results"),
      DT::DTOutput(ns("sim_table")),
      hr(),
      downloadButton(ns("download_table"), "Download Table", class = "btn-sm")
    )
  )
}