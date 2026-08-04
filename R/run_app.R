#' Run the local Shiny application.
#' @param ... Passed to [shiny::runApp()].
#' @export
run_app <- function(...) {
  shiny::runApp(shiny::shinyApp(ui = app_ui(), server = app_server), ...)
}
