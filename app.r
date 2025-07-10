library(shiny)
library(ggplot2)
library(shinyDesign2)

source("modules/data_module.R")
source("modules/analysis_module.R")

# Load reference data ONCE
reference_data <- list(
  "Chicken Heart" = readRDS('C:/Users/xie15/OneDrive/Documents/GitHub/shinyDesign2/inst/extdata/ref_chicken_heart.rds'),
  "Human Brain" = readRDS('C:/Users/xie15/OneDrive/Documents/GitHub/shinyDesign2/inst/extdata/ref_human_brain.rds')
)

ui <- navbarPage("Spatial Sequencing Analysis",
  tabPanel("Data Input",
    dataInputUI("data_input")
  ),
  tabPanel("Analysis",
    analysisUI("analysis")
  )
)

server <- function(input, output, session) {
  data_obj <- dataInputServer("data_input", reference_data)
  analysisServer("analysis", data_obj)
}

shinyApp(ui, server)
