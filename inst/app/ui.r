
ui <- navbarPage("spaDesign: Spatial Transcriptomics Experimental Design",
                 tabPanel("Step1: Data Preparation",
                          dataInputUI("data_input")
                 ),
                 tabPanel("Step2: Sequencing Depth Estimation",
                          analysisUI("analysis")
                 )
)
