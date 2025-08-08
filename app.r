library(shiny)
library(ggplot2)
library(shinyDesign2)
library(DT)
library(shinyBS)
library(hdf5r)
library(scam)
library(plotly)
library(dplyr)
library(parallel)
library(future.apply)


# Load reference data ONCE
reference_data <- list(
  "Chicken Heart" = readRDS('/Users/xie.1097/Documents/GitHub/shinyDesign2/inst/extdata//ref_chicken_heart.rds'),
  "Human Brain" = readRDS('/Users/xie.1097/Documents/GitHub/shinyDesign2/inst/extdata/ref_human_brain.rds')
)

options(shiny.maxRequestSize = 300 * 1024^2)

ui <- navbarPage("Spatial Transcriptomics Experimental Design",
  tabPanel("Data Preparation",
    dataInputUI("data_input")
  ),
  tabPanel("Sequencing Depth Estimation",
    analysisUI("analysis")
  )
)

server <- function(input, output, session) {
  data_obj <- dataInputServer("data_input", reference_data)
  analysisServer("analysis", data_obj)
}

shinyApp(ui, server)
