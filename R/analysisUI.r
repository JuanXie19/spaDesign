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
                   choices = c("None (run original only)" = "original",
                               "Change effect size" = "effect_size",
                               "Add spatial perturbation" = "spatial"),
                   selected = "original"),
      
      conditionalPanel(
        # Condition now depends on the 'effect_size' radio button choice
        condition = sprintf("input['%s'] == 'effect_size'", ns("simulation_type")),
        sliderInput(ns("effect_size"), "Effect size:", min = 1, max = 3, value = 1.5, step = 0.5),
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
      tags$h4("Detailed simulation results"),
      DT::DTOutput(ns("sim_table")),
      hr(),
      downloadButton(ns("download_table"), "Download Table", class = "btn-sm")
    )
  )
}


