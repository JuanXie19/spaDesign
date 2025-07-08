library(shiny)
library(ggplot2)


reference_data <- list(
  "Mouse Brain" = readRDS('C:/Users/xie15/OneDrive/Documents/GitHub/shinyDesign2/inst/extdata/ref_chicken_heart.rds'),
  "Human Brain" = readRDS('C:/Users/xie15/OneDrive/Documents/GitHub/shinyDesign2/inst/extdata/ref_human_brain.rds')
)



  ui <- navbarPage("Spatial Sequencing Analysis",
    tabPanel("Data Input",
      sidebarLayout(
        sidebarPanel(
          radioButtons("data_source", "Choose data source:",
                       choices = c("Use reference data" = "reference",
                                   "Upload your own data" = "upload"),
                       selected = "reference"),
          conditionalPanel(
            condition = "input.data_source == 'reference'",
            selectInput("reference_dataset", "Select reference dataset:",
                        choices = names(reference_data))
          ),
          conditionalPanel(
            condition = "input.data_source == 'upload'",
            fileInput("expr_file", "Upload expression matrix (.rds):"),
            fileInput("coord_file", "Upload coordinates (.rds):"),
            fileInput("anno_file", "Upload annotations (.rds, optional):")
          )
        ),
        mainPanel(
          conditionalPanel(
            condition = "input.data_source == 'reference'",
            verbatimTextOutput("summary_ref"),
            plotOutput("domain_plot")
          ),
          conditionalPanel(
            condition = "input.data_source == 'upload'",
            h4("Your uploaded data summary will appear here.")
          )
        )
      )
    ),
    tabPanel("Analysis",
      h4("This is where you'll add simulation & clustering next.")
    )
  )
  
  server <- function(input, output, session) {
    
    data_obj <- reactive({
      if (input$data_source == "reference") {
        req(input$reference_dataset)
        reference_data[[input$reference_dataset]]
      } else {
        # Placeholder for uploaded data
        validate(need(FALSE, "Please implement uploaded data creation"))
      }
    })
    
    output$summary_ref <- renderPrint({
      req(input$data_source == "reference")
  obj <- req(data_obj())
  
  expr <- obj@refCounts
  coords <- obj@refcolData
  
  ngenes <- nrow(expr)
  nspots <- ncol(expr)
  ndomains <- length(unique(coords$domain))
  domain_counts <- table(coords$domain)
  
  cat("Summary of reference dataset:\n")
  cat("-----------------------------\n")
  cat("Genes:  ", ngenes, "\n")
  cat("Spots:  ", nspots, "\n")
  cat("Domains:", ndomains, "\n\n")
  
  cat("Number of spots per domain:\n")
  print(domain_counts)
    })
    
    output$domain_plot <- renderPlot({
      req(input$data_source == "reference")
      obj <- req(data_obj())
      
      df <- obj@refcolData
      
      
      ggplot(df, aes(x=x, y=y, color=domain)) +
        geom_point(size=2) +
        theme_classic() +
        labs(title = paste("Spatial domains in", input$reference_dataset),
             x = "X coordinate", y = "Y coordinate")
    })
  }
  
  shinyApp(ui, server)