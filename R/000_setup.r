#' @keywords internal
"_PACKAGE"

# Define the global variable in the package namespace
reference_data_paths <- NULL

.onLoad <- function(libname, pkgname) {
  
  # Set Shiny upload limit
  options(shiny.maxRequestSize = 300 * 1024^2)
  
  # Store file paths inside the package namespace
  reference_data_paths <<- list(
    "Chicken Heart" = system.file("extdata/ref_chicken_heart.rds", package = pkgname),
    "Human Brain"   = system.file("extdata/ref_human_brain.rds", package = pkgname)
  )
}



