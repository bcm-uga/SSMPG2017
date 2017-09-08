server <- function(input, output) {
  n.challenge <- reactive({
    n <- 0
    if (input$tab == "challenge"){
      n <- 1
    } else if (input$tab == "challenge_2"){
      n <- 2
    }
    return(n)
  })
  
  gt <- reactive({
    n <- 1
    if (input$tab == "challenge"){
      n <- 1
    } else if (input$tab == "challenge_2"){
      n <- 2
    }
    return(all.gt[[n]])
  })
  
  
  source(file.path("server", "uiOutput.R"), local = TRUE)$value
  source(file.path("server", "infobox.R"), local = TRUE)$value
  source(file.path("server", "users.R"), local = TRUE)$value
  source(file.path("server", "submission.R"), local = TRUE)$value
  source(file.path("server", "tab-barchart.R"), local = TRUE)$value
  source(file.path("server", "tab-table.R"), local = TRUE)$value
  source(file.path("server", "tab-barchart-region.R"), local = TRUE)$value
  source(file.path("server", "tab-table-region.R"), local = TRUE)$value
}
