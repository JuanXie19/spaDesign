##
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
                   choices = c("None (run baseline only)" = "baseline",
                               "Change effect size" = "effect_size",
                               "Add spatial perturbation" = "spatial"),
                   selected = "baseline"),
      
      conditionalPanel(
        # Condition now depends on the 'effect_size' radio button choice
        condition = sprintf("input['%s'] == 'effect_size'", ns("simulation_type")),
        sliderInput(ns("effect_size"), "Effect size:", min = 0.5, max = 5, value = 1, step = 0.5),
        # You can add the info icon and tooltip here if you like
        shinyBS::bsTooltip(id = ns("effect_size"), title = "A value of 1 means original effect size, while 2 means twice the original effect size. Larger values lead to a stronger signal.", placement = "right")
      ),
      
      # Conditional panel for the spatial noise slider
      conditionalPanel(
        # Condition now depends on the 'spatial' radio button choice
        condition = sprintf("input['%s'] == 'spatial'", ns("simulation_type")),
        sliderInput(ns("sigma"), "Sigma(spatial noise):", min = 0.5, max = 3, value = 1, step = 0.5),
        shinyBS::bsTooltip(id = ns("sigma"), title = "Larger sigma leads to more disturbed spatial pattern.", placement = "right")
      ),
      actionButton(ns("run_sim"), "Run simulation", class = "btn-primary")
    ),
    
    # main panel
    mainPanel(
      plotly::plotlyOutput(ns("sim_plot")),
      
      br(),
      # ---- Advanced saturation detection controls moved here ----
      checkboxInput(
        ns("show_saturation_advanced"),
        "Advanced saturation detection options",
        value = FALSE
      ),
      
      conditionalPanel(
        condition = sprintf("input['%s']", ns("show_saturation_advanced")),
        
        sliderInput(
          ns("slope_threshold"),
          "Slope threshold:",
          min = 0.01, max = 0.1, value = 0.05, step = 0.01
        ),
        shinyBS::bsTooltip(
          id = ns("slope_threshold"),
          title = "Slope below this value indicates saturation. Lower values = stricter saturation criteria.",
          placement = "right"
        ),
        
        sliderInput(
          ns("metric_percentage"),
          "Required metric percentage:",
          min = 0.5, max = 1.0, value = 0.8, step = 0.1
        ),
        shinyBS::bsTooltip(
          id = ns("metric_percentage"),
          title = "Saturation is only considered after metric reaches this fraction of maximum. Higher values = more conservative.",
          placement = "right"
        )
      ),
      tags$h4("Detailed simulation results"),
      DT::DTOutput(ns("sim_table")),
      hr(),
      downloadButton(ns("download_table"), "Download Table", class = "btn-sm")
    )
  )
}


