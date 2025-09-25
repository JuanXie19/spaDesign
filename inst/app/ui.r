
ui <- navbarPage("spaDesign: Spatial Transcriptomics Experimental Design",
                 tabPanel("Step1: Data Preparation",
                          dataInputUI("data_input")
                 ),
                 tabPanel("Step2: Sequencing Depth Estimation",
                          analysisUI("analysis")
                 )
)


#' Data Input UI Module
#'
#' UI for the data input tab. Allows users to choose reference data or upload their own SpaceRanger output.
#'
#' @param id A character string specifying the Shiny module namespace.
#'
#' @return A Shiny UI element (sidebarLayout) for data input.
#' @export
dataInputUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns("data_source"), "Choose data source:",
                   choices = c("Use reference data" = "reference",
                               "Upload your own data" = "upload"),
                   selected = "reference"),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'reference'", ns("data_source")),
        selectInput(ns("reference_dataset"), "Select reference dataset:", choices = NULL)
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
        h4('Upload SpaceRanger Output Files'),
        fileInput(ns("h5_file"), "Choose filtered_feature_bc_matrix.h5 file", accept = ".h5"),
        fileInput(ns("positions_file"), "Upload tissue_positions_list.csv file", accept = '.csv'),
        fileInput(ns("anno_file"), "Upload annotations (.csv, optional):"),
        
        # conditional panel for n_clusters.It shows up only if data source is 'upload' AND annotation file is not provided
        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload' && !input['%s']", ns("data_source"), ns("anno_file")),
          numericInput(ns("n_clusters"), "Number of expected clusters:", value = 7, min=2, step=1)
        ),
        
        checkboxInput(ns("show_advanced"), "Show advanced feature selection options", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s']", ns("show_advanced")),
          sliderInput(ns("logfc_cutoff"), "logFC cutoff:", min=0.1, max=2, value=0.7, step=0.1),
          sliderInput(ns("mean_in_cutoff"), "Mean in cutoff:", min=0.5, max=5, value=1.8, step=0.1),
          numericInput(ns("max_num_gene"), "Max number of genes:", value=10, min=1, step=1)
        ),
        
        # action button to triggole the analysis
        hr(),
        actionButton(ns('prepare_data_button'), "Run Data Preparation", class =  "btn-primary")
      )
    ),
    mainPanel(
      plotOutput(ns("domain_plot")),
      verbatimTextOutput(ns("summary_ref"))
    )
  )
}



##
#' Analysis UI Module
#'
#' UI for simulation and analysis tab. Allows effect size and spatial perturbation simulations.
#'
#' @param id Module namespace.
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
        bsTooltip(id = ns("effect_size"), title = "A value of 1 means original effect size, while 2 means twice the original effect size. Larger values lead to a stronger signal.", placement = "right")
      ),
      
      # Conditional panel for the spatial noise slider
      conditionalPanel(
        # Condition now depends on the 'spatial' radio button choice
        condition = sprintf("input['%s'] == 'spatial'", ns("simulation_type")),
        sliderInput(ns("sigma"), "Sigma(spatial noise):", min = 0.5, max = 3, value = 1, step = 0.5),
        bsTooltip(id = ns("sigma"), title = "Larger sigma leads to more disturbed spatial pattern.", placement = "right")
      ),
      actionButton(ns("run_sim"), "Run simulation", class = "btn-primary")
    ),
    
    # main panel
    mainPanel(
      plotlyOutput(ns("sim_plot")),
      tags$h4("Detailed simulation results"),
      DT::DTOutput(ns("sim_table")),
      hr(),
      downloadButton(ns("download_table"), "Download Table", class = "btn-sm")
    )
  )
}



