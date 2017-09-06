output$ui <- renderUI({
  
  if (input$tab == "challenge" || input$tab == "challenge_2") {
    fluidPage(
      title = " ",
      source(file.path("ui", "tab-challenge.R"), local = TRUE)$value
    )
  } else if (input$tab == "vignette") {
    fluidPage(
      title = " ",
      style = "background-color:#ffffff;",
      source(file.path("ui", "tab-vignette.R"), local = TRUE)$value
    )
  } else if (input$tab == "download") {
    fluidPage(
      title = " ",
      source(file.path("ui", "tab-download.R"), local = TRUE)$value    
    )
  }
  
})