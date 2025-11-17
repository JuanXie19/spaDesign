
server <- function(input, output, session) {
  data_obj <- dataInputServer("data_input", reference_data_paths)
  analysisServer("analysis", data_obj)
}



