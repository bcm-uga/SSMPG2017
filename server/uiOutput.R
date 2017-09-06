output$ui <- renderUI({
  
  fluidPage(
    theme = shinythemes::shinytheme("flatly"),
    if (input$tab == "challenge" || input$tab == "challenge_2") {
      source(file.path("ui", "tab-challenge.R"), local = TRUE)$value
    } else if (input$tab == "vignette") {
      source(file.path("ui", "tab-vignette.R"), local = TRUE)$value
    } else if (input$tab == "download") {
      source(file.path("ui", "tab-download.R"), local = TRUE)$value    
    })
  
})