.onLoad <- function(libname, pkgname) {
  options(shiny.maxRequestSize = 300 * 1024^2)
}

# pilot/reference data paths (from your package)
reference_data_paths <- list(
  "Chicken Heart" = system.file("extdata/ref_chicken_heart.rds", package = pkgname),
  "Human Brain"   = system.file("extdata/ref_human_brain.rds", package = pkgname)
)
