#' Launch the spaDesign Shiny app
#'
#' Starts the interactive Shiny application bundled with \pkg{spaDesign}. The app
#' supports built-in reference datasets and user-uploaded SpaceRanger output.
#'
#' @param ... Additional arguments passed to \code{shiny::runApp()}.
#' @import plotly
#' @import shiny
#' @import shinyBS
#' @return This function is called for its side effect of launching a Shiny app.
#' @export
#'
#' @examples
#' \dontrun{
#' runSpaDesignApp()
#' }
runSpaDesignApp <- function(...) {
  app_dir <- system.file("app", package = "spaDesign")
  if (!nzchar(app_dir)) {
    stop("Could not find the bundled spaDesign Shiny app.", call. = FALSE)
  }

  app_env <- new.env(parent = asNamespace("spaDesign"))
  app_env$reference_data_paths <- list(
    "Chicken Heart" = system.file("extdata/ref_chickenHeart.rds", package = "spaDesign"),
    "Human Brain" = system.file("extdata/ref_humanBrain.rds", package = "spaDesign")
  )

  sys.source(file.path(app_dir, "ui.r"), envir = app_env)
  sys.source(file.path(app_dir, "server.r"), envir = app_env)

  shiny::runApp(
    shiny::shinyApp(ui = app_env$ui, server = app_env$server),
    ...
  )
}
