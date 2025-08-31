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
library(pbmcapply)
library(future.apply)
library(BRISC)


# Load reference data ONCE
reference_data <- list(
  "Chicken Heart" = readRDS('/Users/juanxie/Documents/GitHub/shinyDesign2/inst/extdata/ref_chicken_heart.rds'),
  "Human Brain" = readRDS('/Users/juanxie/Documents/GitHub/shinyDesign2/inst/extdata/ref_human_brain.rds')
)

reference_data_paths <- list(
  "Chicken Heart" = system.file("extdata/ref_chicken_heart.rds", package = "shinyDesign2"),
  "Human Brain"   = system.file("extdata/ref_human_brain.rds", package = "shinyDesign2")
)



options(shiny.maxRequestSize = 300 * 1024^2)

ui <- navbarPage("spaDesign: Spatial Transcriptomics Experimental Design",
  tabPanel("Step1: Data Preparation",
    dataInputUI("data_input")
  ),
  tabPanel("Step2: Sequencing Depth Estimation",
    analysisUI("analysis")
  )
)

server <- function(input, output, session) {
  data_obj <- dataInputServer("data_input", reference_data_paths)
  analysisServer("analysis", data_obj)
}

shinyApp(ui, server)
